	subroutine RKPR_version (version, thermo_name)

	character(20)  ::thermo_name, version

	version = "0.3.4"
	thermo_name = "RKPR-EOS"
	return
	endsubroutine RKPR_version
!---------------------------------------------------------------------------------------
!
!
!
!	
	SUBROUTINE HelmRKPR (NDE, NTD, NC, rn, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)
	
	IMPLICIT none
	
	integer, parameter                       :: NCM = 30 
	real(8), parameter                       :: RGAS = 82.05D0 !0.08314472d0
	                                        
	integer, intent(in)                      :: NC, NDE, NTD
	                                       
	real(8), intent(in)                      :: V, T	
	real(8), dimension(NCM), intent(in)      :: rn(NCM)
	
	real(8), intent(out)                     :: Ar, ArTV, ArT, ArTT, ArV, ArV2
	real(8), dimension(NCM), intent(out)     :: Arn, ArVn, ArTn
	real(8), dimension(NCM,NCM), intent(out) :: Arn2
!
!	Internal variables:

	integer                                  :: i, j, NComb

	real(8)                                  :: AUX, auxD2, Bmix, D, D1, D11, D12, D2, dDdT, dDdT2, f, fB, &
	                                            fBD1, fD1, fD1D1, FFB, FFBB, FFBV, fv, fV2, fvD1, g, gV, gV2, &
	                                            TOTN
	real(8), dimension(NCM)                  :: dBi, dDi, dD1i, dDiT
	real(8), dimension(NCM,NCM)              :: aij, daijdT, daijdT2, dBij, dDij, dD1ij

	COMMON /rule/                               ncomb
	
	TOTN = sum(rn)
	call DELTAnder (nc, rn, D1, dD1i, dD1ij)
	D2 = (1D0 - D1)/(1D0 + D1)

	if (ncomb < 2) then
	  
	  call Bnder (NC, rn, Bmix, dBi, dBij)
	  call DandTnder (NTD, NC, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
	  
! 	else
! 	  
! 	  call Bcubicnder (NC, rn, Bmix, dBi, dBij)
! 	  call DCubicandTnder (NTD, NC, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
	  
	end if
!	The f's and g's used here are for Ar, not F (reduced Ar)					
!	This requires to multiply by R all g, f and its derivatives as defined by Mollerup ****
	D12 = D1 - D2
	f = log((V + D1*Bmix)/(V + D2*Bmix))/Bmix/D12
	g = RGAS*log(1D0 - Bmix/V)
	fv = -1D0/((V + D1*Bmix)*(V + D2*Bmix))
	fB = -(f + V*fv)/Bmix
	gv = RGAS*Bmix/(V*(V - Bmix))
	fv2 = (-1D0/(V + D1*Bmix)**2 + 1D0/(V + D2*Bmix)**2)/Bmix/D12
	gv2=RGAS*(1D0/V**2 - 1D0/(V-Bmix)**2)
!	Derivatives of f with respect to delta1
	D11 = 1D0 + D1
	auxD2 = (1D0 + 2D0/D11**2)
	fD1 = (1D0/(V + D1*Bmix) + 2D0/(V + D2*Bmix)/D11**2) - f*auxD2
	fD1 = fD1/D12
	
	fBD1 = -(fB*auxD2 + D1/(V + D1*Bmix)**2 + 2*D2/(V + D2*Bmix)**2/D11**2)
	fBD1 = fBD1/D12
	
	fVD1 = -(fV*auxD2 + 1D0/(V + D1*Bmix)**2 + 2D0/(V + D2*Bmix)**2/D11**2)/D12
	
	fD1D1 = 4*(f - 1D0/(V + D2*Bmix))/D11**3 + Bmix*(-1D0/(V + D1*Bmix)**2             &
	        + 4D0/(V + D2*Bmix)**2/D11**4) - 2*fD1*(1D0 + 2D0/D11**2)
	fD1D1=fD1D1/D12
	
!	Reduced Helmholtz Energy and derivatives
	Ar = -TOTN*g*T - D*f
	ArV = -TOTN*gv*T - D*fv
	ArV2 = -TOTN*gv2*T - D*fv2
!
	AUX = RGAS*T/(V - Bmix)
	FFB = TOTN*AUX - D*fB
	FFBV = -TOTN*AUX/(V - Bmix) + D*(2*fv + V*fv2)/Bmix
	FFBB = TOTN*AUX/(V - Bmix) - D*(2*f + 4*V*fv + V**2*fv2)/Bmix**2
	do i = 1, NC
	  
	  Arn(i) = -g*T + FFB*dBi(i) - f*dDi(i) - D*fD1*dD1i(i)
	  ArVn(i) = -gv*T+FFBV*dBi(i) - fv*dDi(i) - D*fVD1*dD1i(i)
	  IF (NDE == 2) THEN
	    do j = 1, i
	      Arn2(i,j) = AUX*(dBi(i) + dBi(j)) - fB*(dBi(i)*dDi(j) + dBi(j)*dDi(i))      &
	                  + FFB*dBij(i,j) + FFBB*dBi(i)*dBi(j) - f*dDij(i,j)
	      Arn2(i,j) = Arn2(i,j) - D*fBD1*(dBi(i)*dD1i(j) + dBi(j)*dD1i(i))            &
	                  - fD1*(dDi(i)*dD1i(j) + dDi(j)*dD1i(i)) - D*fD1*dD1ij(i,j)      & 
	                  - D*fD1D1*dD1i(i)*dD1i(j)
	      Arn2(j,i)=Arn2(i,j)
	    end do
	  END IF
	end do
!	Temperature derivatives
	IF (NTD == 1) THEN
	
	  ArT=-TOTN*g-dDdT*f
	  ArTV=-TOTN*gv-dDdT*fV
	  ArTT=-dDdT2*f
	  do i = 1, NC
	    
	    ArTn(i) = -g + (TOTN*AUX/T - dDdT*fB)*dBi(i) - f*dDiT(i) - dDdT*fD1*dD1i(i)
     
	  end do
	END IF
	endsubroutine HelmRKPR
!---------------------------------------------------------------------------------------
!
!
!
!
	subroutine readcomp (NC, nin, nout)
	
	implicit DOUBLE PRECISION (A-H, O-Z)
	
	integer, PARAMETER              :: NCM = 30
	real(8), parameter              :: RGAS = 0.08314472d0, Rgas2 = 82.05D0, third = 1.0d0/3
	                               
	integer, intent(in)             :: NC, NIn, NOut

	real(8), dimension(NCM)         :: ac, b, bb1, DCeos, del1, diam, HHA, HHB, HHC, HHD, HHE, &
	                                   HHF, HHG, om, Peneloux, PM, Pc, rk, Tc, vc, vcEOS
	real(8), dimension(NCM,NCM)     :: bij, Kijinf, Kijp, Tstar, lij
	real(8), dimension(NCM,NCM,NCM) :: bijk, Kinf1, Kinf2, K01, K02, lijk
	
	real(8), dimension(2)           :: Apar, Bpar, Cpar
	real(8), dimension(6)           :: Dpar

	CHARACTER(10)                   :: fluid(NCM)
	character(20)                   :: versRKPR

	COMMON /name/                      fluid

	COMMON /COMPONENTS/                ac, b, del1, rk, Kijinf, NTDEP
                                         
	COMMON /rule/                      ncomb
	COMMON /Tdep/                      Kijp, Tstar
	COMMON /COVOL/                     common_b
	COMMON /bcross/                    bij
! 	COMMON /Kcubic/                    Kinf1,Kinf2,K01,K02,Tstar1,Tstar2,C1,C2
! 	COMMON /bcrosscub/                 bijk
      COMMON /CRIT/                      TC, PC, om, DCeos
	common /GCPROM/                    PM, Peneloux, HHA, HHB, HHC, HHD, HHE, HHF, HHG      
	common /versSUB/                   versRKPR
	common /ENTALP/                    iEntalp	

	data Apar /0.0017D0, -2.4407D0/
	data Bpar /1.9681D0, 7.4513D0/
	data Cpar /-2.7238D0, 12.504D0/
	data Dpar / 0.428363D0, 18.496215D0, 0.338426D0, 0.66D0, 789.723105D0, 2.512392D0/

	versRKPR = "0.3.5"
	thermoName = "RKPR EOS"
	Kijinf(:NC,:NC) = 0
	kijp(:NC,:NC) = 0
	Tstar(:NC,:NC) = 1
	lij(:NC,:NC) = 0
	read(NIN,*) ncomb, NTDEP
	do i = 1, NC
	  
	  READ (NIN, *) fluid(i), j
	  if (j == 0) then
	    
! 	    Input are critical properties
	    
	    READ (NIN, *) Tc(i), Pc(i), OM(i), Vc(i), Zrat
	    RT = RGAS*Tc(i)
	    vcEOS(i) = Zrat*vc(i)
	    Zc = Pc(i)*Vceos(i)/RT
	    del1(i) =  Dpar(1) + Dpar(2)*(Dpar(3) - Zc)**Dpar(4) + Dpar(5)*(Dpar(3) - Zc)**Dpar(6)
	    d1 = del1(i) + 1
	    y = 1 + (2*d1)**third + (4D0/d1)**third
	    d1 = (1 + del1(i)*del1(i))/(1 + del1(i))
	    OmegaB = 1D0/(3*y + d1 - 1)
	    OmegaA = (3*y*y + 3*y*d1 + d1*d1 + d1 - 1)*OmegaB**2
	    ac(i) = OmegaA*RT*RT/Pc(i)
	    b(i) = OmegaB*RT/Pc(i)
	    rk(i) = (Apar(2)*Zc + Apar(1))*om(i)**2 + (Bpar(2)*Zc + Bpar(1))*om(i) + (Cpar(2)*Zc + Cpar(1))
	    Zcin = Zc/Zrat
	    Vc(i) = Vceos(i)/Zrat
	    dceos(i) = 1D0/Vceos(i)
	    
	  elseif (j == 1) then
	    
! 	    Input are compound parameters. Critical properties are calculated thereof.
	    
	    READ (NIN, *) ac(i), b(i), del1(i), rk(i)
	    d1 = del1(i) + 1
	    y = 1 + (2*d1)**third + (4D0/d1)**third
	    d1 = (1 + del1(i)*del1(i))/(1 + del1(i))
	    OmegaB = 1D0/(3*y + d1 - 1)
	    OmegaA = (3*y*y + 3*y*d1 + d1*d1 + d1 - 1)*OmegaB**2
	    Zc = y*OmegaB
	    Tc(i) = ac(i)*OmegaB/b(i)/OmegaA/Rgas
	    Pc(i) = OmegaB*Rgas*Tc(i)/b(i)
	    AA = Apar(2)*Zc + Apar(1); BB = Bpar(2)*Zc + Bpar(1); CC = Cpar(2)*Zc + Cpar(1) - rk(i)
	    om(i) = (-BB + dsqrt(Bb*BB - 4*AA*CC))/2/AA
	    
	  endif
	  read (NIn, *) PM(I), Peneloux(i), HHA(I), HHB(I), HHC(I), HHD(I), HHE(i), HHF(i), HHG(i)
	  
	enddo

	write (NOut, '("  ",100("-"))')
	WRITE (NOut, 11)
11	FORMAT(X,'Pure component properties',//,3X,'Compound',6X,'Tc(K)   Pc(bar)   omega  ac(bar L2/mol2)   b(L/mol)     k-value    d1     M(g/mol)  c(cm3/mol)'/)

	DO I = 1, NC
	
	  WRITE (NOut, 12) fluid(I), TC(I), PC(I), om(i), ac(i), b(i), rk(I), del1(i), PM(i), Peneloux(i)
	  
	enddo
12	FORMAT (3X, A10, X, 2(F8.2,2X), F6.4, 5X, G12.7, G12.4, 2X, 2(F7.4, 4X), 2(F7.2,4X))
	iEntalp = NC
	do i = 1, NC
	  if ((HHA(i) == 0) .AND. (HHB(i) == 0) .AND. (HHC(i) == 0) .AND. (HHD(i) == 0) .AND. (HHE(i) == 0) .AND. (HHF(i) == 0)) then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	if (iEntalp /= 0) then

	  write (NOut, '(/," Passut-Danner constant for ideal gas enthalpy calculation:",/,"   Compound   A(BTU/lb)  B(BTU/(lb 튣))  ", &
     &             "C(BTU/(lb 튣2)�1E3  D(BTU/(lb 튣3))�1E6  E(BTU/(lb 튣4))�1E10  F(BTU/(lb 튣5))�1E14     G(BTU/(lb 튣))",/,135("-"))')
	  do i = 1,NC

	    write (NOut, 14) fluid(i), HHA(i), HHB(i), 1.D3*HHC(i), 1.D6*HHD(i), 1D10*HHE(i), 1D14*HHF(i), HHG(i)

	  enddo

	endif
14	format (3X, A10, X, F9.6, 3X, F9.6, 10X, F9.6, 12X, F9.6, 13X, F9.6, 13X, F9.6, 12X, F9.6)
	
! 	Units inside GC package are (K, atm, cm3, mol)
	Pc(:NC) = Pc(:NC)*Rgas2/Rgas/1e3
	b(:NC) = b(:NC)*1D3
	ac(:NC) = ac(:NC)*Rgas2/Rgas*1D3
	
! 	Binary interaction parameters
	read (Nin, *) NownIJ
	do k = 1, NownIJ
	  
	  if (NTDep == 0) then
	    
	    read (NIn, *) i, j, omegab, d1
	    kijinf(i,j) = d1
	    kijinf(j,i) = d1
	    lij(i,j) = omegab
	    lij(j,i) = omegab
	    
	  else
	    
	    read (NIn, *) i, j, omegab, d1, omegaa, y
	    kijinf(i,j) = d1
	    kijinf(j,i) = d1
	    lij(i,j) = omegab
	    lij(j,i) = omegab
	    kijp(i,j) = omegaa
	    kijp(j,i) = omegaa
	    Tstar(i,j) = y
	    Tstar(j,i) = y
	    
	  endif
	  
	enddo
	  
! 	do i = 1, NC

! 	  IF(i.gt.1)then
! 	    if(ncomb.lt.2)then
! 	      READ(NIN,*) (Kijinf(j,i),j=1,i-1)
! 	      if(NTDEP >= 1)READ(NIN,*) (Kijp(j,i),j=1,i-1)
! 	      if(NTDEP >= 1)READ(NIN,*)(Tstar(j,i),j=1,i-1)		
! 	      READ(NIN,*) (lij(j,i),j=1,i-1)
! ! 	    else
! ! 	      READ (NIN, *) K01,K02
! ! 	      if (NTDEP >= 1) READ(NIN,*) Kinf1,Kinf2
! ! 	      if (NTDEP == 1) READ(NIN,*)Tstar1,Tstar2
! ! 	      if (NTDEP == 2) READ(NIN,*)C1,C2
! ! 	      READ(NIN,*) Lijk(1,1,2),Lijk(1,2,2)
! 	    end if
! 	  ENDIF
! 	end do



	  
	write (NOut, 115) fluid(:NC)
115	FORMAT (/,1X,'Binary interaction parameter matrix at {T -> infinity}, [k^inf]', //, 18X, <NC>A10)
116	FORMAT (3X, A10, 2X, '|',<NC>(F8.4, 2X),'|')
	do i = 1, NC

	  write (NOut, 116) fluid(i), kijinf(I,:I)

	enddo
	write (NOut, 117) fluid(:NC)
117	FORMAT (/,1X,'Binary interaction parameter matrix at {T = 0 K}, [k^0]', //, 18X, <NC>A10)
	do i = 1, NC

	  write (NOut, 116) fluid(i), kijp(I,:I)

	enddo 
	write (NOut, 118) fluid(:NC)
118	FORMAT (/,1X,'Reference temperature matrix, [T*] (K)', //, 18X, <NC>A10)
	do i = 1, NC

	  write (NOut, 116) fluid(i), Tstar(I,:I)

	enddo
	write (NOut, 119) fluid(:NC)
119	FORMAT (/,1X,'Repulsive interaction parameter, [l]', //, 18X, <NC>A10)
	do i = 1, NC

	  write (NOut, 116) fluid(i), lij(I,:I)

	enddo 	
	  
! 	if (ncomb < 2) then
! 	  write (NOUT, *)'  kij matrix'
! 	  if(NTDEP.EQ.0)then
! 	    write(NOUT,*)'    K12 = ',Kijinf(1,2)
! 	    write(NOUT,*)
! 	  else
! 	    write(NOUT,*)' Kinf12 = ',Kijinf(1,2)
! 	    write(NOUT,*)
! 	    write(NOUT,*)'     K` = ',Kijp(1,2)
! 	    write(NOUT,*)
! 	    write(NOUT,*)'     T* = ',Tstar(1,2)
! 	    write(NOUT,*)
! 	  end if		
! ! 		DO I=1,NC
! ! 		write(NOUT,6)FLUID(I),(Kij(j,i),j=1,i-1)
! ! 		END DO
! 	  write(NOUT,*)
! 	  write(NOUT,*)'  LIJ MATRIX'
! 	  DO I=1,NC
! 	    write(NOUT,6)FLUID(I),(Lij(j,i),j=1,i-1)
! 	  END DO
	  
! 	else
! 	  
! 	  if(NTDEP.EQ.0)then
! 	    
! 	    write(NOUT,*)' Kijk:     112      122'
! 	    write(NOUT,7)K01,K02
! 	    write(NOUT,*)
! 	    
! 	  else
! 	    
! 	    write(NOUT,*)' K0ijk:    112      122'
! 	    write(NOUT,7)K01,K02
! 	    write(NOUT,*)
! 	    write(NOUT,*)'Kinfijk:   112      122'
! 	    write(NOUT,7)Kinf1,Kinf2
! 	    write(NOUT,*)
! 	    
! 	  end if
! 	  if(NTDEP.EQ.2)then
! 	    
! 	    write(NOUT,*)' Cijk:     112      122'
! 	    write(NOUT,7)C1,C2
! 	    write(NOUT,*)
! 	    
! 	  end if
! 	  write(NOUT,*)' Lijk:     112      122'
! 	  write(NOUT,7)Lijk(1,1,2),Lijk(1,2,2)
! 	  write(NOUT,*)
	  
! 	end if
	write (NOUT, '(/, " Combining rules:")')
	if (ncomb == 0) then
	write(NOUT,*)' 0: Classical or van der Waals '
	  do i = 1, NC
	    do j = i, NC
		
	      bij(i,j) = (1 - lij(i,j))*(b(i) + b(j))/2
	      bij(j,i) = bij(i,j)
		
	    end do
	  end do
! 	elseif (ncomb == 3) then
! 	  
! 	  do i = 1, NC
! 	    
! 	    bijk(i,i,i)=b(i)
! 	    do j=i + 1,NC
! 		
! 	      bijk(i,i,j)=(1-lijk(i,i,j))*(2*b(i)+b(j))/3
! 	      bijk(i,j,i)=bijk(i,i,j)
! 	      bijk(j,i,i)=bijk(i,i,j)
! 	      bijk(i,j,j)=(1-lijk(i,j,j))*(b(i)+2*b(j))/3
! 	      bijk(j,i,j)=bijk(i,j,j)
! 	      bijk(j,j,i)=bijk(i,j,j)
! 	      do k=j+1,nc	! only possible with three or more components
! 		  
! 	        bijk(i,j,k)=(1-lijk(i,j,k))*(b(i)+b(j)+b(k))/3
! 	        bijk(j,i,k)=bijk(i,j,k)
! 	        bijk(i,k,j)=bijk(i,j,k)
! 	        bijk(j,k,i)=bijk(i,j,k)
! 	        bijk(k,i,j)=bijk(i,j,k)
! 	        bijk(k,j,i)=bijk(i,j,k)
! 		  
! 	      end do
! 	      
! 	    end do
! 	    
! 	  end do
	  
	else
	  
	  write (NOUT, *) ' 1: Lorentz-Berthelot'
	  do i = 1, NC
	    
	    diam(i) = b(i)**third
	    
	  enddo
	  do i = 1, NC
	    do j = i, NC
		
	      bij(i,j) = ((1D0 - lij(i,j))*(diam(i) + diam(j))/2)**3
	      bij(j,i) = bij(i,j)
		
	    end do
	  end do
	  
	end if
 6	FORMAT(A10,4F8.5)
 7	FORMAT(9x,F7.4,2x,F7.4)
	endsubroutine readcomp
!-----------------------------------------------------------------------------------------
!
!	Then Kij values will be called indicating the lower index first, e.g. Kij(1,3)
!
!
	subroutine aTder (ac, Tc, rk, T, a, dadT, dadT2)
!       implicit DOUBLE PRECISION (A-H,O-Z)
	implicit none
	
	real(8), intent(in) :: ac, rk, T, Tc
	real(8), intent(out) :: a, dadT, dadT2
	
	integer              :: NModel
	real(8)              :: Tr
!	Given ac,Tc and the k parameter of the RKPR correlation, as well as the actual T,
!	this subroutine calculates a(T) and its first and second derivatives with T.
!
	COMMON /MODEL/ NMODEL
	
	Tr = T/Tc
	IF (NMODEL < 3) THEN
	
	  a = ac*(1 + rk*(1 - dsqrt(Tr)))**2
	  dadT = ac*rk*(rk - (rk + 1)/dsqrt(Tr))/Tc
	  dadT2 = ac*rk*(rk + 1)/(2*Tc**2*Tr**1.5D0)

	ELSEif (nModel == 4) then
	
	  a = ac*(3.D0/(2 + Tr))**rk
	  dadT = -rk*a/Tc/(2 + Tr)
	  dadT2 = -(rk + 1)*dadT/Tc/(2 + Tr)

	ELSEif (nModel == 6) then
	
	  a = ac
	  dadT = 0
	  dadT2 = 0
	  
	END IF
	endsubroutine aTder
!-------------------------------------------------------------------------------
!
!
!
	subroutine aijTder (NTD, NC, T, aij, daijdT, daijdT2)
	
	implicit none !DOUBLE PRECISION (A-H,O-Z)
	integer, PARAMETER                       :: NCM = 30
                                               
	integer, intent(in)                      :: NTD, NC
                                               
	real(8), intent(in)                      :: T
      
	real(8), dimension(NCM,NCM), intent(out) :: aij, daijdT, daijdT2
      
	integer                                  :: i, j
	real(8)                                  :: aux, barrgij, ratK
      
	integer                                  :: NTDep, NComb
	real(8), dimension(NCM)                  :: ac, ai, daidT, daidT2, b, d1, rhoc, omega, Pc, rk, Tc
	real(8), dimension(NCM,NCM)              :: bij, kij, kijinf, kijp, lij, Tstar
      
! 	DOUBLE PRECISION Kij0(nco,nco),Kij(nco,nco)
      COMMON /CRIT/                               TC, PC, omega, rhoc
	COMMON /COMPONENTS/                         ac, b, d1, rk, Kijinf, NTDep
	COMMON /bcross/                             bij
	COMMON /rule/                               ncomb
	COMMON /Tdep/                               kijp, Tstar
	
	Kij(:NC,:NC) = 0.0D0
	IF (NTDEP >= 1) THEN
	
	  do i = 1, NC - 1
	    do j = i+1, NC
		
	      Kij(j,i) = Kijinf(j,i) + kijp(j,i)*dexp(-T/Tstar(j,i))
	      Kij(i,j) = Kij(j,i)
		
	    enddo
	  enddo

	ELSE
	
	  Kij = Kijinf
	  
	END IF
	DO i = 1, NC
	
	  call aTder (ac(i), Tc(i), rk(i), T, ai(i), daidT(i), daidT2(i))
	  aij(i,i) = ai(i)
	  daijdT(i,i) = daidT(i)
	  daijdT2(i,i) = daidT2(i)
	  IF (i > 1) THEN
	    do j = 1, i - 1
		
	      aij(j,i) = dsqrt(ai(i)*ai(j))*(1-Kij(j,i))
	      aij(i,j) = aij(j,i)
	      if (NTD == 1) then
		  
	        daijdT(j,i) = (1D0 - Kij(j,i))*( dsqrt(ai(i)/ai(j))*daidT(j)               &
	                                        + dsqrt(ai(j)/ai(i))*daidT(i) )/2
	        daijdT(i,j) = daijdT(j,i)
	        daijdT2(j,i) = (1D0 - Kij(j,i))*( daidT(j)*daidT(i)/dsqrt(ai(i)*ai(j))     &
	                           + dsqrt(ai(i)/ai(j))*(daidT2(j) - daidT(j)**2/(2*ai(j)))&
	                           + dsqrt(ai(j)/ai(i))*(daidT2(i)-daidT(i)**2/(2*ai(i))) )/2
	        daijdT2(i,j)=daijdT2(j,i)
		  
	      end if
	      
	    end do
	  END IF
	  
	END DO
	if (ncomb == 1) then
	  DO i = 1, NC - 1
	    DO j = i + 1, NC
	    
	      barrgij = bij(i,j)/dsqrt(b(i)*b(j))
	      aij(i,j) = barrgij*aij(i,j)
	      aij(j,i) = aij(i,j)
	      daijdT(i,j) = barrgij*daijdT(i,j)
	      daijdT(j,i) = daijdT(i,j)
	      daijdT2(i,j) = barrgij*daijdT2(i,j)
	      daijdT2(j,i) = daijdT2(i,j)
		
	    END DO
	  END DO
	end if
	IF(NTDEP >= 1.and.NTD == 1) THEN
	  do i = 1, NC
	    do j = i+1, NC
	      aux = daijdT(j,i)
! 	      ratK=Kij(1,2)/(1-Kij(1,2))/Tstar
	      ratK = (Kij(j,i) - Kijinf(j,i))/(1 - Kij(j,i))/Tstar(j,i)
	      daijdT(j,i) = aux + aij(j,i)*ratK
	      daijdT(i,j) = daijdT(j,i)
	      daijdT2(j,i) = daijdT2(j,i) + (2*aux - aij(j,i)/Tstar(j,i))*ratK  ! 2* was missing (before aux)
!	    										since implementation in 2005	(06-03-2008)
	      daijdT2(i,j) = daijdT2(j,i)
	    enddo
	  enddo
	END IF
	return
	
	endsubroutine aijTder
!-------------------------------------------------------------------------------
!
!
!
!
!
!
!
! 	subroutine aijkTder(NTD,NC,T,aijk,daijkdT,daijkdT2)
! 	implicit DOUBLE PRECISION (A-H,O-Z)
! 	PARAMETER (nco=2)
! 	DOUBLE PRECISION Kinf1,Kinf2,K01,K02,Kijk(nco,nco,nco),Kij(nco,nco)
! 	dimension ai(nco),daidT(nco),daidT2(nco)
! 	dimension aijk(nco,nco,nco),daijkdT(nco,nco,nco),daijkdT2(nco,nco,nco)
!       COMMON/CRIT/TC(nco),PC(nco),DC(nco)
! 	COMMON /COMPONENTS/ ac(nco),b(nco),d1(nco),rk(nco),Kij,NTDEP
! !	COMMON /Tdep/ Kinf, Tstar
! 	COMMON /Kcubic/Kinf1,Kinf2,K01,K02,Tstar1,Tstar2,C1,C2
! 	COMMON /bcross/bij(nco,nco)
! 	Kijk=0.0D0
! 	IF(NTDEP.EQ.1)THEN
! 	  Kijk(1,1,2)=Kinf1+K01*dexp(-T/Tstar1)
! 	  Kijk(1,2,1)=Kijk(1,1,2)
! 	  Kijk(2,1,1)=Kijk(1,1,2)
! 	  Kijk(2,2,1)=Kinf2+K02*dexp(-T/Tstar2)
! 	  Kijk(2,1,2)=Kijk(2,2,1)
! 	  Kijk(1,2,2)=Kijk(2,2,1)
! 	ELSE IF(NTDEP.EQ.2)THEN
! 	  A1=K01
! 	  B1=Kinf1
! 	  SUM1=(C1+T)
! 	  rT1=T/SUM1
! 	  Kijk(1,1,2)=A1+B1*rT1
! 	  Kijk(1,2,1)=Kijk(1,1,2)
! 	  Kijk(2,1,1)=Kijk(1,1,2)
! 	  A2=K02
! 	  B2=Kinf2
! 	  SUM2=(C2+T)
! 	  rT2=T/SUM2
! 	  Kijk(2,2,1)=A2+B2*rT2
! 	  Kijk(2,1,2)=Kijk(2,2,1)
! 	  Kijk(1,2,2)=Kijk(2,2,1)
! 	ELSE
! 	  Kijk(1,1,2)=K01
! 	  Kijk(1,2,1)=K01
! 	  Kijk(2,1,1)=K01
! 	  Kijk(2,2,1)=K02
! 	  Kijk(2,1,2)=K02
! 	  Kijk(1,2,2)=K02
! 	END IF
! 	third=1.0d0/3
! 	DO i=1,nc
! 	  call aTder(ac(i),Tc(i),rk(i),T,ai(i),daidT(i),daidT2(i))
! 	  aijk(i,i,i)=ai(i)
! 	  daijkdT(i,i,i)=daidT(i)
! 	  daijkdT2(i,i,i)=daidT2(i)
! 	  IF (i.gt.1) THEN
! 	    do j=1,i-1
! 	      aijk(j,j,i)=(1-Kijk(j,j,i))*(ai(i)*ai(j)*ai(j))**third
! 	      aijk(j,i,j)=aijk(j,j,i)
! 	      aijk(i,j,j)=aijk(j,j,i)
! 	      aijk(j,i,i)=(1-Kijk(j,i,i))*(ai(i)*ai(i)*ai(j))**third
! 	      aijk(i,i,j)=aijk(j,i,i)
! 	      aijk(i,j,i)=aijk(j,i,i)
! 	      if(NTD.EQ.1)then
! 	        daijkdT(j,j,i)=(2*daidT(j)/ai(j)+
!      &	                               daidT(i)/ai(i))*aijk(j,j,i)/3
! 	        daijkdT(j,i,j)=daijkdT(j,j,i)
! 	        daijkdT(i,j,j)=daijkdT(j,j,i)
! 	        daijkdT(j,i,i)=(daidT(j)/ai(j)+
!      &	                             2*daidT(i)/ai(i))*aijk(j,i,i)/3
! 	        daijkdT(i,i,j)=daijkdT(j,i,i)
! 	        daijkdT(i,j,i)=daijkdT(j,i,i)
! 	        daijkdT2(j,j,i)=daijkdT(j,j,i)**2 / aijk(j,j,i) +
!      1                               (2*daidT2(j)/ai(j)+daidT2(i)/ai(i)-
!      1        (2*(daidT(j)/ai(j))**2+(daidT(i)/ai(i))**2))*aijk(j,j,i)/3
! 	        daijkdT2(j,i,j)=daijkdT2(j,j,i)
! 	        daijkdT2(i,j,j)=daijkdT2(j,j,i)
! 	        daijkdT2(j,i,i)=daijkdT(j,i,i)**2 / aijk(j,i,i) +
!      1                               (daidT2(j)/ai(j)+2*daidT2(i)/ai(i)-
!      1        ((daidT(j)/ai(j))**2+2*(daidT(i)/ai(i))**2))*aijk(j,i,i)/3
! 	        daijkdT2(i,i,j)=daijkdT2(j,i,i)
! 	        daijkdT2(i,j,i)=daijkdT2(j,i,i)
! 	      end if
! 	      IF (j.gt.1) THEN	! only possible with three or more components
! 	        do k=1,j-1
! 	          aijk(k,j,i)=(1-Kijk(k,j,i))*(ai(i)*ai(j)*ai(k))**third
! 	          aijk(k,i,j)=aijk(k,j,i)
! 	          aijk(j,i,k)=aijk(k,j,i)
! 	          aijk(j,k,i)=aijk(k,j,i)
! 	          aijk(i,j,k)=aijk(k,j,i)
! 	          aijk(i,k,j)=aijk(k,j,i)
! 	          if(NTD.EQ.1)then
! 	        ! add here daijkdT(k,j,i) & daijkdT2(k,j,i)
! 	          end if
! 	        end do
! 	      END IF
! 	    end do
! 	  END IF
! 	END DO
! 	IF(NTDEP.ge.1.and.NTD.EQ.1)THEN
! 	  aux1=daijkdT(1,1,2)
! 	  aux2=daijkdT(1,2,2)
! 	  IF(NTDEP.EQ.1)THEN
! 	    ratK1=(Kijk(1,1,2)-Kinf1)/(1-Kijk(1,1,2))/Tstar1
! 	    ratK2=(Kijk(1,2,2)-Kinf2)/(1-Kijk(1,2,2))/Tstar2
! c	first
! 	    daijkdT(1,1,2)=aux1+aijk(1,1,2)*ratK1
! 	    daijkdT(1,2,2)=aux2+aijk(1,2,2)*ratK2
! c	second
! 	    daijkdT2(1,1,2)=daijkdT2(1,1,2)+(2*aux1-aijk(1,1,2)/Tstar1)
!      &                                                            *ratK1
! 	    daijkdT2(1,2,2)=daijkdT2(1,2,2)+(2*aux2-aijk(1,2,2)/Tstar2)
!      &                                                            *ratK1	    
! 	  ELSE IF(NTDEP.EQ.2)THEN
! 	    rB1=B1/SUM1
! 	    dK1dT=rB1*(1-rT1)
! 	    d2K1dT2=-2*dK1dT/SUM1
! 	    rB2=B2/SUM2
! 	    dK2dT=rB2*(1-rT2)
! 	    d2K2dT2=-2*dK2dT/SUM2
! c	first
! 	    daijkdT(1,1,2)=aux1-dK1dT*aijk(1,1,2)/(1-Kijk(1,1,2))
! 	    daijkdT(1,2,2)=aux2-dK2dT*aijk(1,2,2)/(1-Kijk(1,2,2))
! c	second
! 	    daijkdT2(1,1,2)=daijkdT2(1,1,2)
!      1	         -(2*aux1*dK1dT+aijk(1,1,2)*d2K1dT2)/(1-Kijk(1,1,2))
! 	    daijkdT2(1,2,2)=daijkdT2(1,2,2)
!      1               -(2*aux2*dK2dT+aijk(1,2,2)*d2K2dT2)/(1-Kijk(1,2,2))
! 	  END IF
! 	  daijkdT(1,2,1)=daijkdT(1,1,2)
! 	  daijkdT(2,1,1)=daijkdT(1,1,2)
! 	  daijkdT(2,2,1)=daijkdT(1,2,2)
! 	  daijkdT(2,1,2)=daijkdT(1,2,2)
! 	  daijkdT2(1,2,1)=daijkdT2(1,1,2)
! 	  daijkdT2(2,1,1)=daijkdT2(1,1,2)
! 	  daijkdT2(2,2,1)=daijkdT2(1,2,2)
! 	  daijkdT2(2,1,2)=daijkdT2(1,2,2)
! 	END IF
! 	end subroutine aijkTder
!-----------------------------------------------------------------------
!
	subroutine DandTnder (NTD, NC, T, rn, D, dDi, dDiT, dDij, dDdT, dDdT2)
	
      implicit DOUBLE PRECISION (A-H,O-Z)
	
	integer, PARAMETER                       :: NCM = 30
	                                         
	integer, intent(in)                      :: NTD, NC
	real(8), intent(in)                      :: rn(NCM), T
	                                         
	real(8), intent(out)                     :: D, dDdT, dDdT2
	real(8), dimension(NCM), intent(out)     :: dDi, dDiT
	real(8), dimension(NCM,NCM), intent(out) :: dDij
	
	real(8), dimension(NCM,NCM)              :: aij, daijdT, daijdT2
	
	call aijTder (NTD, NC, T, aij, daijdT, daijdT2)
	D = 0D0
	dDdT = 0D0
	dDdT2 = 0.0D0
	DO i = 1, NC
	
	  aux = 0.0D0
	  aux2 = 0.0D0
	  dDi(i) = 0.0D0
	  dDiT(i) = 0.0D0
	  do j = 1, nc
	    
	    dDi(i)=dDi(i) + 2*rn(j)*aij(i,j)
	    if (NTD == 1) then
	
	      dDiT(i) = dDiT(i) + 2*rn(j)*daijdT(i,j)
	      aux2 = aux2 + rn(j)*daijdT2(i,j)
	
	    end if
	    dDij(i,j) = 2*aij(i,j)
	    aux = aux + rn(j)*aij(i,j)
	    
	  end do
	  D = D + rn(i)*aux
	  if (NTD == 1) then
	    
	    dDdT = dDdT + rn(i)*dDiT(i)/2
	    dDdT2 = dDdT2 + rn(i)*aux2
	    
	  end if
	  
	END DO
	endsubroutine DandTnder
!----------------------------------------------------------------------

! 	subroutine DcubicandTnder (NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
!       implicit DOUBLE PRECISION (A-H,O-Z)
!       PARAMETER (nco=2)
! 	dimension rn(nco),dDiT(nco)
! 	dimension dDi(nco),dDij(nco,nco),aux(nco),auxij(nco,nco)
! 	dimension auxT(nco),auxTij(nco,nco),auxT2(nco),auxT2ij(nco,nco)
! 	dimension aijk(nco,nco,nco),daijkdT(nco,nco,nco),
!      &          daijkdT2(nco,nco,nco)
!      
! 	call aijkTder(NTD,nc,T,aijk,daijkdT,daijkdT2)
! 	TOTN = sum(rn)
! 	D=0.0D0
! 	dDdT=0.0D0
! 	dDdT2=0.0D0
! 	aux=0.0D0
! 	auxij=0.0D0
! 	auxT=0.0D0
! 	auxTij=0.0D0
! 	DO i=1,nc
! 	  do j=1,nc
! 	    do k=1,nc
! 	      auxij(i,j)=auxij(i,j)+rn(k)*aijk(i,j,k)
! 	      if(NTD.EQ.1)then
! 	        auxTij(i,j)=auxTij(i,j)+rn(k)*daijkdT(i,j,k)
! 	        auxT2ij(i,j)=auxT2ij(i,j)+rn(k)*daijkdT2(i,j,k)
! 	      end if
! 	    end do
! 	    aux(i)=aux(i)+rn(j)*auxij(i,j)
! 	    if(NTD.EQ.1)then
! 	      auxT(i)=auxT(i)+rn(j)*auxTij(i,j)
! 	      auxT2(i)=auxT2(i)+rn(j)*auxT2ij(i,j)
! 	    end if
! 	  end do
! 	  D=D+rn(i)*aux(i)
! 	  if(NTD.EQ.1)then
! 	    dDdT=dDdT+rn(i)*auxT(i)
! 	    dDdT2=dDdT2+rn(i)*auxT2(i)
! 	  end if
! 	END DO
! 	D=D/TOTN
! 	if(NTD.EQ.1)then
! 	  dDdT=dDdT/TOTN
! 	  dDdT2=dDdT2/TOTN
! 	end if
! 	DO i=1,nc
! 	  dDi(i)=(3*aux(i)-D)/TOTN
! 	  if(NTD.EQ.1)dDiT(i)=(3*auxT(i)-dDdT)/TOTN
! 	  do j=1,i
! 	    
! 	    dDij(i,j)=(6*auxij(i,j)-dDi(i)-dDi(j))/TOTN
! 	    dDij(j,i)=dDij(i,j)
! 	    
! 	  end do
! 	END DO
! 	endsubroutine DcubicandTnder
!------------------------------------------------------------------------------------- 
	subroutine DELTAnder (NC, rn, D1m, dD1i, dD1ij)
	
	implicit none !DOUBLE PRECISION (A-H,O-Z)
	integer, PARAMETER          :: NCM = 30
	                            
	integer, intent(in)         :: NC
	real(8), intent(in)         :: rn(NCM)
	                            
	real(8), intent(out)        :: D1m, dD1i(NCM), dD1ij(NCM,NCM)
	                            
	integer                     :: i, j
	real(8)                     :: totn
	
	integer                     :: NTDep
	real(8), dimension(NCM)     :: ac, b, d1, rk
	real(8), dimension(NCM,NCM) :: kij

	COMMON /COMPONENTS/            ac, b, d1, rk, Kij, NTDEP
	
	TOTN = sum(rn)
	D1m = dot_product(rn,d1)/totn
	do i = 1, NC
	  
	  dD1i(i) = (d1(i) - D1m)/totn
	  do j = 1, NC
	    
	    dD1ij(i,j) = (2*D1m-d1(i)-d1(j))/totn**2
	    
	  end do
	  
	end do
	endsubroutine DELTAnder
!-------------------------------------------------------------------------------------
!
!
	subroutine Bnder (nc, rn, Bmix, dBi, dBij)
	implicit DOUBLE PRECISION (A-H,O-Z)
      integer, PARAMETER   :: NCM = 30
                           
      integer, intent(in)  :: NC
      real(8), intent(in)  :: rn(NCM)
	
	real(8), intent(out) :: Bmix, dBi(NCM), dBij(NCM,NCM)
	
	real(8)              :: aux(NCM), bij(NCM,NCM)
	COMMON /bcross/         bij
	
	TOTN = sum(rn)
	Bmix = 0.0D0
	aux = 0.0D0
	DO i = 1, NC
	  do j = 1, NC
	    aux(i) = aux(i) + rn(j)*bij(i,j)
	  end do
	  Bmix = Bmix + rn(i)*aux(i)
	END DO
	Bmix = Bmix/totn
	DO i = 1, NC
	  dBi(i) = (2*aux(i) - Bmix)/totn
	  do j = 1, i
	    dBij(i,j) = (2*bij(i,j) - dBi(i) - dBi(j))/totn
	    dBij(j,i) = dBij(i,j)
	  end do
	END DO
	endsubroutine Bnder
! 
! 	subroutine Bcubicnder(nc,rn,Bmix,dBi,dBij)
!       implicit DOUBLE PRECISION (A-H,O-Z)
!       PARAMETER (nco=2)
! 	dimension rn(nco),dBi(nco),dBij(nco,nco),aux(nco),auxij(nco,nco)
! 	COMMON /bcrosscub/bijk(nco,nco,nco)
! 	TOTN = sum(rn)
! 	sqn=TOTN*TOTN
! 	Bmix=0.0D0
! 	aux=0.0D0
! 	auxij=0.0D0
! 	DO i=1,nc
! 	  do j=1,nc
! 	    do k=1,nc
! 	      auxij(i,j)=auxij(i,j)+rn(k)*bijk(i,j,k)
! 	    end do
! 	    aux(i)=aux(i)+rn(j)*auxij(i,j)
! 	  end do
! 	  Bmix=Bmix+rn(i)*aux(i)
! 	END DO
! 	Bmix=Bmix/sqn
! 	DO i=1,nc
! 	  dBi(i)=(3*aux(i)-2*totn*Bmix)/sqn
! 	  do j=1,i
! 	    dBij(i,j)=(6*auxij(i,j)-2*(Bmix+totn*dBi(i)+totn*dBi(j)))/sqn
! 	    dBij(j,i)=dBij(i,j)
! 	  end do
! 	END DO
! 	endsubroutine Bcubicnder
!-----------------------------------------------------------------------------------



	subroutine SL_version (version, thermo_name)

	character(20)  ::thermo_name, version

	version = "0.2.0"
	thermo_name = "SL-EOS"
	return
	endsubroutine SL_version
!---------------------------------------------------------------------------------------
!
!
!(NDER, NTD, NC, rn(:NC), V, T, Ar, ArV, ArTV, ArV2, Arn(:NC), ArVn(:NC), ArTn(:NC), Arn2(:NC,:NC), ArT, ArTT)
!	
subroutine Helmhotz_SL (n_deriv, t_deriv, NC, n, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArT2)
  implicit none
  integer, parameter  :: NCM = 30
  real(8), parameter  :: Rgas = 82.05 !0.08314472d0
  
  integer, intent(in) :: n_deriv, NC, t_deriv
  real(8), intent(in) :: n(NCM), V, T
  
  real(8), intent(out) :: Ar, ArTV, ArV, ArV2 , ArT, ArT2
  real(8), intent(out), dimension(NCM) :: Arn, ArTn, ArVn
  real(8), intent(out), dimension(NCM,NCM) :: Arn2

  integer :: i, j, NComb
  
  real(8) :: amix, bmix, da_dT, d2a_dT2, f, fB, fBR, fBV, fR, fV, fV2, g, gB, gBB, gBR, gBV, gR, gRR, gV, gV2, gVR, ntot, rmix
  real(8), dimension(NCM) :: da_dn, d2a_dTdn, db_dn, rseg
  real(8), dimension(NCM,NCM) :: d2a_dn2, d2b_dn2

  COMMON /rule/   NComb
  common /SL_par/ rseg  
  
  ntot = sum(n)
  rmix = dot_product(n(:NC), rseg(:NC))
  if (NComb < 2) then
    call Bnder (NC, n, Bmix, dB_dn, d2B_dn2)
    call DandTnder (t_deriv, NC, T, n, amix, da_dn, d2a_dTdn, d2a_dn2, da_dT, d2a_dT2)
!   else
!     call Bcubicnder(nc,rn,Bmix,dBi,dBij)
!     call DCubicandTnder(NTD,nc,T,rn,D,dDi,dDiT,dDij,dDdT,dDdT2)
  end if
  
  f = 1D0/V
  g = Rgas * rmix*(1 + (V/bmix - 1)*dlog(1 - bmix/V)) !units of "mol"
  
  fV = -f/V
  fV2 = -2*fV/V
  fB = 0

  gV = Rgas*rmix*(dlog(1 - bmix/V)/bmix + (V/bmix - 1)*bmix/(V**2*(1 - bmix/V)))
  gV2 = Rgas*rmix*(2D0/(V**2*(1 - bmix/V)) - (2*(V/bmix - 1))*bmix/(V**3*(1 - bmix/V)) - (V/bmix - 1)*bmix**2/(V**4*(1 - bmix/V)**2))
  gB = Rgas*rmix*( - V*dlog(1 - bmix/V)/bmix**2 - (V/bmix - 1)/(V*(1 - bmix/V)))
  gBB = Rgas*rmix*(2*v*dlog(1-bmix/v)/bmix**3+2/(bmix**2*(1-bmix/v))-(v/bmix-1)/(v**2*(1-bmix/v)**2))
  gBV = Rgas*rmix*( - 2/(V*(1 - bmix/V)*bmix) - dlog(1 - bmix/V)/bmix**2 + (V/bmix - 1)/(V**2*(1 - bmix/V)) + (V/bmix - 1)*bmix/(V**3*(1 - bmix/V)**2))
  gR = g/rmix
  gRR = 0
  gVR = gV/rmix
  gBR = gB/rmix
  
! Reduced Helmholtz Energy and derivatives (units of "mol"
  Ar = g*T - amix*f
  ArV = gv*T - amix*fv
  ArV2 = gv2*T - amix*fv2  
  
  do i = 1, NC
  
    Arn(i) = T*(gB*dB_dn(i) + gR*rseg(i)) - f*da_dn(i)
    ArVn(i) = T*(gBV*dB_dn(i) + gVR*rseg(i)) - fV*da_dn(i)
    if (n_deriv == 2) then
      do j = 1, i
        Arn2(j,i) = T*((gBB*dB_dn(j) + gBR*rseg(j) )*dB_dn(i) + gB*d2B_dn2(j,i) + ( gBR*dB_dn(j) + gRR*rseg(j) )*rseg(i) + gR*0) - f*d2a_dn2(j,i)
        Arn2(i,j) = Arn2(j,i)
      end do
    end if
  end do
  
! Temperature derivatives
  IF (t_deriv == 1) THEN
    ArT = ntot*g - da_dT*f !da/dT = 0; this could be simplified but I'm keeping it only if it's modified in the near future.
    ArTV = ntot*gV - da_dT*fV
    ArT2 = -d2a_dT2*f
    do i = 1, NC
    
      ArTn(i) = (gB*dB_dn(i) + gR*rseg(i)) - f*d2a_dTdn(i)
	
    end do
  END IF  
  return
endsubroutine Helmhotz_SL
!--------------------------------------------------------------------
!
!
!
subroutine read_SL_parameters (inputFile, outputFile, NC)

  implicit none
  integer, parameter          :: NCM = 30
  real(8), parameter          :: Zc1 = 2*dlog(2D0) - 1, Rgas = 0.08314472d0, Rgas2 = 82.05D0, third = 1D0/3
  character, parameter        :: TAB = achar(9)
                              
  integer, intent(in)         :: inputFile, outputFile
                              
  integer                     :: i, j, k, compound, NC, NownIJ, NTDEP = 0, ncomb = 0
   
  real(8)                     :: epsilon_R, r1, OmegaA, OmegaB, Zc
  real(8), dimension(3)       :: OA, OB, OZ
  real(8), dimension(NCM)     :: a, b, Bgpec, dc, del1, m, MW, omega, P_str, Pc, rho_str, rseg, T_str, Tc
  real(8), dimension(NCM,NCM) :: bij, kijinf, kijp, lij = 0, Tstar 
  
  character(10)               :: compound_name(NCM)
  COMMON /name/ compound_name
  COMMON /COVOL/ Bgpec
  COMMON /bcross/ bij
  COMMON /COMPONENTS/ a, b, del1, m, Kijinf, NTDEP  
  COMMON /Tdep/ kijp, Tstar  
  common /SL_par/ rseg
  COMMON /CRIT/ Tc, Pc, omega, Dc  

  DATA OA /0.45D0, -0.0758921D0, -0.06794D0/
  DATA OB /0.465D0, 0.0214636D0, 0.0442549D0/
  DATA OZ /1D-2,   -0.0073054D0, -0.019315D0/

  write (outputFile, '( "Model: Sanchez-Lacombe EOS (1974/1991?)", //, "Input model values:")')
  write (outputFile, '("Fluid", A, "T*(K)", A, "P*(bar)", A, "rho*(kg/m3)", A, "MW(kg/kmol)", A, "v*(m3/kmol)")')  (TAB, i = 1, 5)
  do compound = 1, NC
    
    read (inputFile, *) compound_name(compound)
    read (inputFile, *) T_str(compound), P_str(compound), rho_str(compound), MW(compound)
    write (outputFile, '(A10, 5(A, G12.4))') compound_name(compound), TAB, T_str(compound), TAB, P_str(compound), TAB, rho_str(compound), TAB, MW(compound), TAB, Rgas*T_str(compound)/P_str(compound)
 
 !   Conversion from "literature values" to "typical model values"
    epsilon_R = T_str(compound)
    b(compound) = MW(compound)/rho_str(compound)
    Bgpec(compound) = b(compound)
    rseg(compound) = Mw(compound)*P_str(compound)/Rgas/T_str(compound)/rho_str(compound)
    a(compound) = rseg(compound)*epsilon_R*Rgas*b(compound)
    r1 = (rseg(compound)-1)/rseg(compound)
    OmegaA = rseg(compound)**OA(1) * (Zc1 + OA(2)*r1 + OA(3)*r1*r1)
    OmegaB = (Zc1/2 + OB(2)*r1 + OB(3)*r1*r1) / rseg(compound)**OB(1)
    Zc = rseg(compound)**OZ(1) * (Zc1 + OZ(2)*r1 + OZ(3)*r1*r1)
    Tc(compound) = a(compound)*OmegaB/b(compound)/OmegaA/Rgas
    Pc(compound) = OmegaB*Rgas*Tc(compound)/b(compound)
    Dc(compound) = Pc(compound)/Zc/Rgas/Tc(compound)

  enddo
  write (outputFile, '(/, "Specific model values:", /, "Fluid", A, "r", A, "b(L/mol)", A, "epsilon/kB(K)", A, "a(bar L2/mol2)")') (TAB, i = 1, 4)
  write (outputFile, '(<NC>( A10, 4(A, G12.6), /) )') (compound_name(i), TAB, rseg(i), TAB, b(i), TAB, T_str(i), TAB, a(i), i = 1, NC)

! Units inside GC package are (K, atm, cm3, mol)
  Pc(:NC) = Pc(:NC)/1.01325
  b(:NC) = b(:NC)*1D3
  a(:NC) = a(:NC)*Rgas2/Rgas*1D3
  
! Binary interaction parameters:  
!  write (outputFile, '(/, "Attractive binary interaction parameters [kinf]:")')
!  kijinf = 0
!  do compound = 1, NC - 1
!    read (inputFile, *) kijinf(compound+1:NC,compound)
!    write (outputFile, '(A10, <NC>(A, G12.4))') compound_name(compound), (TAB,  kijinf(i,compound), i = compound + 1, NC)
!    kijinf(compound,compound+1:NC) =  kijinf(compound+1:NC,compound)
!  enddo
!
!  write (outputFile, '(/, "Repulsive/free-volume binary interaction parameters: [l]")')
!  do compound = 1, NC - 1
!    read (inputFile, *) l(compound+1:NC,compound)
!    write (outputFile, '(A10, <NC>(A, G12.4))') compound_name(compound), (TAB,  l(i,compound), i = compound + 1, NC)
!    l(compound,compound+1:NC) =  kijinf(compound+1:NC,compound)
!  enddo
! 	Binary interaction parameters
  read (inputFile, *) NownIJ
  do k = 1, NownIJ
    
    if (NTDep == 0) then
      
      read (inputFile, *) i, j, lij(i,j), kijinf(i,j)
      kijinf(j,i) = kijinf(i,j)
      lij(j,i) = lij(i,j)
      
    else
      
      read (inputFile, *) i, j, lij(i,j), kijinf(i,j), kijp(i,j), Tstar(i,j)

      kijinf(j,i) = kijinf(i,j)
      lij(j,i) = lij(i,j)
      kijp(j,i) = kijp(i,j)
      Tstar(j,i) = Tstar(i,j)
      
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



	  
	write (outputFile, 115) compound_name(:NC)
115	FORMAT (/,1X,'Binary interaction parameter matrix at {T -> infinity}, [k^inf]', //, 18X, <NC>A10)
116	FORMAT (3X, A10, 2X, '|',<NC>(F8.4, 2X),'|')
	do compound = 1, NC

	  write (outputFile, 116) compound_name(compound), kijinf(compound,:compound)

	enddo
	write (outputFile, 117) compound_name(:NC)
117	FORMAT (/,1X,'Binary interaction parameter matrix at {T = 0 K}, [k^0]', //, 18X, <NC>A10)
	do compound = 1, NC

	  write (outputFile, 116) compound_name(compound), kijp(compound,:compound)

	enddo 
	write (outputFile, 118) compound_name(:NC)
118	FORMAT (/,1X,'Reference temperature matrix, [T*] (K)', //, 18X, <NC>A10)
	do compound = 1, NC

	  write (outputFile, 116) compound_name(compound), Tstar(compound,:compound)

	enddo
	write (outputFile, 119) compound_name(:NC)
119	FORMAT (/,1X,'Repulsive interaction parameter, [l]', //, 18X, <NC>A10)
	do compound = 1, NC

	  write (outputFile, 116) compound_name(compound), lij(compound,:compound)

	enddo 	
	  
! 	if (ncomb < 2) then
! 	  write (outputFile, *)'  kij matrix'
! 	  if(NTDEP.EQ.0)then
! 	    write(outputFile,*)'    K12 = ',Kijinf(1,2)
! 	    write(outputFile,*)
! 	  else
! 	    write(outputFile,*)' Kinf12 = ',Kijinf(1,2)
! 	    write(outputFile,*)
! 	    write(outputFile,*)'     K` = ',Kijp(1,2)
! 	    write(outputFile,*)
! 	    write(outputFile,*)'     T* = ',Tstar(1,2)
! 	    write(outputFile,*)
! 	  end if		
! ! 		DO I=1,NC
! ! 		write(outputFile,6)compound_name(i),(Kij(j,i),j=1,i-1)
! ! 		END DO
! 	  write(outputFile,*)
! 	  write(outputFile,*)'  LIJ MATRIX'
! 	  DO I=1,NC
! 	    write(outputFile,6)compound_name(i),(Lij(j,i),j=1,i-1)
! 	  END DO
	  
! 	else
! 	  
! 	  if(NTDEP.EQ.0)then
! 	    
! 	    write(outputFile,*)' Kijk:     112      122'
! 	    write(outputFile,7)K01,K02
! 	    write(outputFile,*)
! 	    
! 	  else
! 	    
! 	    write(outputFile,*)' K0ijk:    112      122'
! 	    write(outputFile,7)K01,K02
! 	    write(outputFile,*)
! 	    write(outputFile,*)'Kinfijk:   112      122'
! 	    write(outputFile,7)Kinf1,Kinf2
! 	    write(outputFile,*)
! 	    
! 	  end if
! 	  if(NTDEP.EQ.2)then
! 	    
! 	    write(outputFile,*)' Cijk:     112      122'
! 	    write(outputFile,7)C1,C2
! 	    write(outputFile,*)
! 	    
! 	  end if
! 	  write(outputFile,*)' Lijk:     112      122'
! 	  write(outputFile,7)Lijk(1,1,2),Lijk(1,2,2)
! 	  write(outputFile,*)
	  
! 	end if
	write (outputFile, '(/, " Combining rules:")')
	if (ncomb == 0) then
	  write(outputFile,*)' 0: Classical or van der Waals '
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
  	  
    write (outputFile, *) ' 1: Lorentz-Berthelot'
    do i = 1, NC
    	    
      del1(i) = b(i)**third !diameters is stored now in del1
    	    
    enddo
    do i = 1, NC
      do j = i, NC
      		
        bij(i,j) = ((1D0 - lij(i,j))*(del1(i) + del1(j))/2)**3
        bij(j,i) = bij(i,j)
      		
      end do
    end do
  	  
  end if



! Segment number "r" is stored in the corresponding "delta_1" parameter (RK-PR) space:
  del1(:NC) = rseg(:NC)
!! Calculation of cross covolumes
!  do i = 1, NC
!    do j = i, NC
!	
!      bij(i,j) = (1 - l(i,j))*(b(i) + b(j))/2
!      bij(j,i) = bij(i,j)
!	
!    end do
!  end do  

  return
endsubroutine read_SL_parameters

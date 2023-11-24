	subroutine PCSAFT_version (submodel, version, thermo_name)

	implicit none
	character(20)  ::thermo_name, version
	integer        ::submodel

	version = "0.5.0"
	if (submodel == 0) then
	  
	  thermo_name = "CS-SAFT"
	  
	elseif (submodel == 1) then
	  
	  thermo_name = "sPC-SAFT"
	  
	elseif (submodel == 2) then
	  
	  thermo_name = "PC-SAFT"
	  
	endif
	return
	endsubroutine PCSAFT_version
!---------------------------------------------------------------------------------------
!	
	subroutine readPCSAFT (nin, nout, NC)
	
	implicit DOUBLE PRECISION (A-H,O-Z)
	
	integer, parameter :: NCM = 30
	real(8), parameter :: NAVO = 6.026D23, PI6 =3.14159265D0/6.D0, RGAS = 82.05D0 !0.08314472d0
	
	
!	SUBMODEL:   0 = CARNAHAN/STARLING, 1 = CS G-FUNCTION, 2 = ORG. PC-SAFT
	iNTEGER SUBMOD
      
	integer, intent(in) :: NC, NIn, NOut
      
	DOUBLE PRECISION Kij, lij(NCM,NCM)
	CHARACTER(10) fluid(NCM)

	real(8), dimension(:,:), allocatable :: sigav
	
	COMMON /PCSAFTDATA/ RM(NCM), KIJ(NCM,NCM), ETA(NCM), SIGMA(NCM), SIGMAIJ(NCM,NCM), SUBMOD
	COMMON /name/ fluid
	COMMON /CRIT/ TC(NCM),PC(NCM),om(NCM)
	COMMON /covol/ B(NCM)
	common /assoc/ NST
	common /GCPROM/ PMM(NCM),Pen(NCM),HHA(ncm),HHB(ncm),hhc(ncm),HHD(ncm),HHE(NCM),HHF(NCM),HHG(NCM)	
	

	allocate(sigav(NC,NC))

!	read(NIN,*) SUBMOD
	SUBMOD = 2
	Zrat = 1
	AVOTRD = (NAVO*PI6/RGAS)**(1.D0/3)
	write (NOut, '(//, " Pure component properties", //, &
&	               "  Compound", 6x, "Tc(K)", 2X, "Pc(atm)", 3X, "omega", 3X, "epsilon/k(K)  sigma(A)", 4X, "m", 6X, "M(g/mol)   c(cm3/mol)", /)')
	do i = 1, NC
	  
	  READ (NIN,'(A)') fluid(i)
	  write (*, *) fluid(i)
	  READ (NIN,*) Tc(i), Pc(i), OM(i)
	  RT = RGAS*Tc(i)
	  
	  READ (NIN, *) ETA(i), SIGMA(i), RM(i)
	  
	  read (NIn, *) PMM(I), Pen(i), HHA(I), HHB(I), HHC(I), HHD(I), HHE(i), HHF(i), HHG(i)	  

	  write (NOut, '(4X, A10, X, 2(F7.1, X), 2X, F6.4, 3X, F8.2, 6X, 2(F6.4, 3X), X,F6.2, 6X, F6.2)') &
	                               fluid(i), Tc(i), Pc(i), om(i), eta(i), sigma(i), rm(i), PMM(i), Pen(i)

! 	  Complex unit change: [A] -> [dm/mol^(1/3)] -> [(K/atm)^(1/3)]
! 	  SIGMA(I) = 1D-9*SIGMA(I)*AVOTRD 

! 	  Complex unit change: [A] -> [cm/mol^(1/3)] -> [(K/atm)^(1/3)]
	  SIGMA(I) = 1D-8*SIGMA(I)*AVOTRD
	  Kij(i,:NC) = 0.0D0
! 	  Compound co-volume:
	  B(i) = RM(i)*SIGMA(i)**3
	  
	end do
!
!	Constantes de Passut-Danner: impresiónn sólo si al menos se le ingresa valores a un compuesto:
	iEntalp = NC
	do i = 1,NC
	  if ((HHA(i) == 0) .AND. (HHB(i) == 0) .AND. (HHC(i) == 0) .AND. (HHD(i) == 0) .AND. (HHE(i) == 0) .AND. (HHF(i) == 0)) then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	if (iEntalp /= 0) then

	  write (NOut, '(/," Passut-Danner constant for ideal gas enthalpy calculation:",/,"   Compound   A(BTU/lb)  B(BTU/(lb ºR))  ", &
     &             "C(BTU/(lb ºR2)·1E3  D(BTU/(lb ºR3))·1E6  E(BTU/(lb ºR4))·1E10  F(BTU/(lb ºR5))·1E14     G(BTU/(lb ºR))",/,148("-"))')
	  do i = 1,NC

	    write (NOut, 14) fluid(i),HHA(i),HHB(i),1.D3*HHC(i),1.D6*HHD(i),1D10*HHE(i),1D14*HHF(i),HHG(i)

	  enddo

	endif
	
! 	Binary interaction parameters
	read (NIn, *) NownIJ
	do k = 1, NownIJ
	  
	  read (NIn, *) i, j, omegab, d1
	  kij(i,j) = d1
	  kij(j,i) = d1
	  lij(i,j) = omegab
	  lij(j,i) = omegab
	    		  
	enddo
	write (NOut, 115) fluid(:NC)
115	FORMAT (/,1X,'Binary interaction parameter matrix, [k]', //, 18X, <NC>A10)
116	FORMAT (3X, A10, 2X, '|',<NC>(F8.4, 2X),'|')
	do i = 1, NC

	  write (NOut, 116) fluid(i), kij(I,:I)

	enddo
	write (NOut, 119) fluid(:NC)
119	FORMAT (/,1X,'Repulsive interaction parameter, [l]', //, 18X, <NC>A10)
	do i = 1, NC

	  write (NOut, 116) fluid(i), lij(I,:I)

	enddo 	
!
!	SIGMA-MATRIX
	DO I=1,NC
	  DO J=I,NC
	  
	    SIGAV(j,i) = (.5D0*(SIGMA(I)+SIGMA(J))*(1 - lij(j,i)))**3
	    sigav(i,j) = sigav(j,i)
	    SIGMAIJ(J,I) = SIGAV(j,i)*RM(j)*RM(I)
	    SIGMAIJ(I,J) = SIGMAIJ(J,I)
	    
	  ENDDO
	ENDDO
! 
! 	Read association contribution, if any
	read (NIn, *) NST
	if (NST > 0) then
	  
	  sigav(:NC,:NC) = sigav(:NC,:NC)/PI6
	  call SAFTPar (NIn, NOut, model, NC, NST, fluid(:NC), sigav(:NC,:NC))
	  
	endif
	return
14	format(3X,A10,X,F9.6,3X,F9.6,10X,F9.6,12X,F9.6,13X,F9.6,13X,F9.6,12X,F9.6)	
	endsubroutine readPCSAFT
!
!	Then Kij values will be called indicating the lower index first, e.g. Kij(1,3)
!--------------------------------------------------------------------------------------------------
!
!
!
	SUBROUTINE HelmPCSAFT (NDE, NTD, NC, X, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2, ArT, ArTT)
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	
	integer, PARAMETER :: NCM = 30
	real(8), parameter :: RGAS = 82.05D0 !0.08314472d0)
	
	dimension X(NC),Arn(NC),ArVn(NC),ArTn(NC),Arn2(NC,NC)
	DIMENSION AR_N (NC),AR_VN(NC),AR_TN(NC),AR_NN(NC,NC)
	DIMENSION XD1(NC),XD2(NC),XAR1(NCM,NCM),XAR2(NCM,NCM)
	DIMENSION DIA(NCM)
	DIMENSION DD(0:3,NCM), DDT(0:3,NCM),DTT(0:3,NCM)
	SAVE XAR1,XAR2,DIA,DD,DDT,DTT
	DIMENSION S(0:3),ST(0:3),STT(0:3)
	COMMON NHelm
	
	
!	GET T-DEPENDENT PROPERTIES
	NHelm = 1
	CALL DIAMET (NC, T, DIA, DD, DDT, DTT, NSUB)
	NHelm = 0
	CALL SAFT_COEF (NC, T, XAR1, XAR2)

!	GET BASE ZETA-COEFFICIENTS
	DO K = 0,3
      
	  S(K) = 0.D0
	  DO L=1,NC
	  
	    S(K) = S(K) + DD(K,L)*X(L)
	    
	  ENDDO
	  
	ENDDO

!	GET SUMS
	CALL XSUMS(NC,X,XAR1,XAR2,SS1,SS2,XD1,XD2)
	
	TINV2 = 1.D0/T**2
	
! 	A-derivatives
	RT=RGAS*T
	VR=V/RGAS
	CALL DERIVS (NSUB, NC, NDE, NTD, T, VR, X, DIA, DD, DDT, DTT, S, ST, STT, &
	             XAR1, XAR2, XD1, XD2, SS1, SS2, AR,                          &
	             AR_V, AR_T, AR_N, AR_VV, AR_VT, AR_TT, AR_VN, AR_TN, AR_NN)
	Ar=RT*AR
	ArV=RT*AR_V/RGAS
	Arn=RT*AR_N
	ArV2=RT*AR_VV/RGAS/RGAS
	ArVn=RT*AR_VN/RGAS
	IF (NDE == 2) THEN
	
	  Arn2=RT*AR_NN
	  
	END IF
!	Temperature derivatives
	IF (NTD == 1) THEN
	
	  ArT=Ar/T-RT*TINV2*AR_T
	  ArTV=ArV/T-RT*TINV2*AR_VT/RGAS
!	  ArTT= __?_ AR_TT:  1/T-DERIV. OF AR_T
	  ArTn=Arn/T-RT*TINV2*AR_TN
	  
	END IF

	endsubroutine
!--------------------------------------------------------------------------------------------------
!
!
!
	SUBROUTINE DERIVS (NSUB, NC, NDER, NTEMP, T, V, X, DIA, DD, DDT, DTT, &
	                   S, ST, STT, AR1, AR2, XD1, XD2, SS1, SS2,          &
	                   AR, AR_V, AR_T, AR_N, AR_VV, AR_VT, AR_TT,         &
	                   AR_VN, AR_TN, AR_NN)

!	The subroutine DERIVS calculaters the reduced residual
!	helmholz energy for the PC-SAFT EOS, and its derivatives wrt.
!	volume, temperature and composition

!	Input:
!	NC:          NO. of components in mixture
!	T:           Temperature (K)
!	V:           Volume/R --------------> PAY ATTENTION!!!
!	X:           Molar composition
!	DIA:
!	DD:
!	DDT:
!	S:
!	ST:
!	STT:
!	AR1:
!	AR2:
!	XD1:
!	XD2:
!	SS1:
!	SS2:
!
!
!	Output:
!
!	AR:          A^RES/RT
!	AR_V:        V-DERIV. OF AR
!	AR_T:        1/T-DERIV. OF AR
!	AR_N:        N-DERIV. OF AR
!	AR_VV:       V-DERIV. OF AR_V
!	AR_VT:       V-DERIV. OF AR_T
!	AR_TT:       1/T-DERIV. OF AR_T
!	AR_VN        V-DERIV OF AR_N
!	AR_TN        1/T-DERIV OF AR_N
!	AR_NN        N-DERIV OF AR_N
!
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (NCM = 30)
	
	DIMENSION AR_N(NC),AR_VN(NC),AR_TN(NC),AR_NN(NC,NC)
	DIMENSION X(NC)
	DIMENSION DIA(NCM)
	DIMENSION DD(0:3,NCM), DDT(0:3,NCM),DTT(0:3,NCM)
	DIMENSION S(0:3),ST(0:3),STT(0:3)
	DIMENSION AR1(NCM,NCM),AR2(NCM,NCM)
	DIMENSION XD1(NC),XD2(NC)
	DIMENSION DQDN(NC),DQNV(NC),DQNT(NC),DQ2DN2(NC,NC)
	dimension dFasn(NC), dFasTn(NC), dFasVn(NC), dFasnm(NC,NC)
	common /assoc/ NST
!	Repulsive contributions
	IF (NSUB  ==  0) THEN
      
	  CALL CARNAH (NC, NDER, NTEMP, X, DD, DDT, DTT, S, ST, STT, V,   &
	               Q, QV, QVV, QT, QVT, QTT, DQDN, DQ2DN2, DQNV, DQNT)
	  
	ELSE
	
	  CALL HARDCN (NSUB, NC, NDER, NTEMP, X, DIA, DD, DDT, DTT, S, ST, STT, V, &
	               Q, QV, QVV, QT, QVT, QTT, DQDN, DQ2DN2, DQNV, DQNT)
	  
	ENDIF

!	Attractive contribution
	CALL VATTR (NC, NDER, NTEMP, X, T, DD, DDT, S, ST, STT, SS1, SS2, XD1, XD2, &
	           AR1, AR2, V, AR, AR_V, AR_VV, AR_T, AR_VT, AR_TT, AR_N, AR_VN, AR_NN, AR_TN)

!	Add up
	AR = AR + Q
	AR_V = AR_V + QV
	AR_VV = AR_VV + QVV
	IF (NTEMP > 0) THEN
	
	  AR_T = AR_T + QT
	  AR_VT = AR_VT + QVT
	  
	  IF (NTEMP > 1) AR_TT = AR_TT + QTT
	  
	  DO I=1,NC
	  
	    AR_TN(I) = AR_TN(I) + DQNT(I)
	  
	  ENDDO
	  
	ENDIF
	
	DO I = 1,NC
	
	  AR_N(I) = AR_N(I) + DQDN(I)
	  AR_VN(I) = AR_VN(I) + DQNV(I)
	   
	ENDDO
	IF (NDER > 1) THEN
	  DO I = 1,NC
	    DO J=1,I
	    
	      AR_NN(J,I) = AR_NN(J,I) + DQ2DN2(J,I)
	
	    ENDDO
	  ENDDO
	ENDIF
	if (NST > 0) then
	  
	  call SAFT_Helmholtz (NSUB, NC, NST, NDER, 2, NTemp, T, V, x, &
	                       Fas, dFasV, dFasVV, dFasTV, dFasVn, dFasn, dFasnm, dFasTn, dFasT, dFasTT)
	  AR = AR + Fas
	  AR_V = AR_V + dFasV
	  AR_VV = AR_VV + dFasVV
	  if (NTemp > 0) then
	    
! 	    This subroutine works with derivatives wrt (1/T) instead of T.
! 	    AR_T = d(AR/RT)/d(1/T)
!	    dFasT = d(Aassoc/RT)/dT
	    T2 = T*T
	    AR_T = AR_T - T2*dFasT
	    AR_VT = AR_VT - T2*dFasTV
	    AR_TN(:NC) = AR_TN(:NC) - T2*dFasTn(:NC)
	    if (NTemp > 1) AR_TT = AR_TT + 2*T2*T*dFasT - T2*T2*dFasTT
	    
	  endif
	  if (NDer > 0) then	    
	    do i = 1, NC
		
	      AR_VN(i) = AR_VN(i) + dFASVn(i)
	      AR_N(i) = AR_N(i) + dFasn(i)
	      if (NDer > 1) then
	        do j = 1, i
		    
	          AR_NN(j,i) = AR_NN(j,i) + dFasnm(j,i)
	          AR_NN(i,j) = AR_NN(j,i)
		    
	        enddo
	      endif
	    enddo
	  endif
	  
	endif	
	endsubroutine DERIVS
!--------------------------------------------------------------------------------------------------
!
!
!	
	SUBROUTINE DIAMET (NC, T, DIA, DD, DDT, DTT, NSUB)
	
!	This subroutine calculates the diameters and their powers,
!	the factors M(I)*DIA(I)**K, K=0,1,2,3
!	and the first and second derivatives wrt. 1/T

!	INPUT:
!
!	NC:          NO. OF COMPONENTS IN MIXTURE
!	T:           TEMPERATURE (K)
!
!	OUTPUT:
!
!	DIA:         (VECTOR OF) DIAMETERS
!	DD:          DD(K,I) IS M(I)*D(I)**K,  K=0,1,2,3
!	DDT:         1ST 1/T-DERIVATIVE OF DD
!	             0-ELEMENT IS D DIA D (1/T)
!	DTT:         1ST 1/T-DERIVATIVE OF DDT
!	SUBMOD:      EOS TYPE, 0 = FULL CS, 1 = G FROM CS, 2 = ORG

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	
	PARAMETER (CF = 0.12D0, NCM = 30)
	
	DIMENSION DD(0:3,NCM),DDB(0:3,NCM),DDT(0:3,NCM),DTT(0:3,NCM),DIA(NCM)
	DOUBLE PRECISION M,KIJ
	INTEGER SUBMOD
	COMMON /forB/       DDB
	COMMON /PCSAFTDATA/ M(NCM),KIJ(NCM,NCM),ETA(NCM),SIGMA(NCM), SIGMAIJ(NCM,NCM),SUBMOD
	COMMON NHelm
	SAVE TOLD
	DATA TOLD /0.D0/
	
	NSUB = SUBMOD
	
	IF (T == TOLD) RETURN
	
	IF (NHelm == 1) TOLD = T
	TINF = 1.D0/T

!	Calculate T-dependent diameter and its powers (diam),
!	multiplied by M, and the derivatives wrt 1/T

	DO I=1,NC
	
	  EFAC = -3.D0*ETA(I)
	  CTRM = CF*EXP(EFAC*TINF)
!	  DIAMETER AND 1/T-DERIVATIVES
	  DI   =  SIGMA(I)*(1.D0-CTRM)
	  DIT  = -SIGMA(I)*CTRM*EFAC
	  DIT2 = DIT*EFAC
	  DITSQ =DIT*DIT
	  DIA(I) = DI
	  DDT(0,I) = DIT
	  DD(0,I) = M(I)
!  	  POWERS OF D
	  DO K=1,3
	  
	    DD(K,I) = DD(K-1,I)*DI
	    
	  ENDDO
!	  1ST 1/T-DERIVATIVES
	  DDT(1,I) = DIT*M(I)
	  DDT(2,I) = DDT(1,I)*2.0D0*DI
	  DDT(3,I) = DDT(2,I)*1.5D0*DI
! 	  2ND 1/T-DERIVATIVES
	  DTT(0,I) = DIT2
	  DTT(1,I) = DIT2*M(I)
	  DTT(2,I) = 2.D0*M(I)*(DI*DIT2 + DITSQ)
	  DTT(3,I) = 3.D0*DD(1,I)*(DI*DIT2 + 2.D0*DITSQ)
	  
	ENDDO
	DDB=DD
	endsubroutine DIAMET
!--------------------------------------------------------------------------------------------------
!
!
!	
	SUBROUTINE SAFT_COEF (NC, T, AR1, AR2)
!
!	THe subroutine saft_sums calculates the coefficients of the quadratic
!	sums of SAFT
!
!	Input:
!
!	NC:     No. of conponents
!	T:      Temperature (K)
!
!	Output:
!	AR1:   1st array
!	AR2:   2nd array
!
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	
	PARAMETER (NCM = 30)
	
	INTEGER SUBMOD
	DIMENSION AR1(NCM,*),AR2(NCM,*)
	DIMENSION ETAT(NC)
	DOUBLE PRECISION M,KIJ
	COMMON /PCSAFTDATA/ M(NCM),KIJ(NCM,NCM),ETA(NCM),SIGMA(NCM),SIGMAIJ(NCM,NCM),SUBMOD
	SAVE TOLD
	DATA TOLD/0.D0/
	
	IF (T == TOLD) RETURN
	
	TOLD = T
!	Coefficients for quadratic sums
	TINV = 1.D0/T
	DO I=1,NC
	
	  ETAT(I) =SQRT(ETA(I)*TINV)
	  
	ENDDO

!	Coefficient matrices for double sums

!	NOTE: SIGMAIJ is multiplied by M(I) M(J)
!	      and that twice the values of the coefficients are stored
	DO I=1,NC
	  DO J=I,NC
	  
	    COEF = ETAT(I)*ETAT(J)*(1.D0-KIJ(I,J))
	    AR1(J,I) = 2.D0*COEF*SIGMAIJ(J,I)
	    AR2(J,I) = COEF*AR1(J,I)
	    AR1(I,J) = AR1(J,I)
	    AR2(I,J) = AR2(J,I)
	    
	  ENDDO
	ENDDO
	
	endsubroutine SAFT_COEF
!--------------------------------------------------------------------------------------------------
!
!
!
	SUBROUTINE INTCON (M, ZETA, FA, FAZ, FAM, FAZZ, FAMM, FAZM, FB, FBZ, FBM, FBZZ, FBMM, FBZM)
!
!	This routine calculates the integral contributions
!	in PC-SAFT, "I1" and "M C_1_INV I2"
!	
!	Input:
!	
!	M:      segment m-parameter
!	ZETA:   zeta-parameter, S_3/V
!	
!	Output:
!	
!	FA:     1st integral, I1
!	FAZ:    zeta-derivative of FA
!	FAM:    M_INV derivative of FA
!	FAZZ:   zeta-derivative of FAZ
!	FAMM:   M_INV deribvative of FAM
!	FAZM:   M_INV Deribvative of FAZ
!	
!	FB:     2nd integral, I2, multiplied by M/C1
!	FBZ:    zeta-derivative of FB
!	FBM:    M_INV derivative of FB
!	FBZZ:   zeta-derivative of FBZ
!	FBMM:   M_INV deribvative of FBM
!	FBZM:   M_INV deribvative of FBZ
!	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION M
!	
!	Get integrals; FA
	CALL QPOLS (M, ZETA, FA, FAM, FAMM, FAZ, FAZZ, FAZM, FBX, FBMX, FBMMX, FBZX, FBZZX, FBZMX)

!	Get C-contribution
	CALL QDIFC (ZETA, M, F, FZ, FM, FZZ, FZM)

!	Restore values
	VALC = 1.D0/F
	FB = FBX*VALC
	FBZ = VALC*(FBZX - FB*FZ)
	FBM = VALC*(FBMX - FB*FM)
	FBMM = VALC*(FBMMX-2.D0*FM*FBM)
	FBZZ = VALC*(FBZZX - FB*FZZ - 2.D0*FZ*FBZ)
	FBZM = VALC*(FBZMX - FB*FZM - FZ*FBM - FM*FBZ)
	
	endsubroutine INTCON
!--------------------------------------------------------------------------------------------------
!
!
!
	SUBROUTINE INTCON_V (M_INV, VA, VB, ZETA, FA, FAZ, FAZZ, FB, FBZ, FBZZ)
!	
!	This routine calculates the integral contributions
!	in PC-SAFT, I1 and M C_1_INV I2
!	
!	Input:
!	
!	M_INV:  Inverse segment M-parameter
!	ZETA:   ZETA-parameter, S_3/V
!	
!	Output:
!	
!	FA:     1st integral, I1
!	FAZ:    ZETA-derivative of FA
!	FAM:    M_INV derivative of FA
!	FAZZ:   ZETA-derivative of FAZ
!	FAMM:   M_INV deribvative of FAM
!	FAZM:   M_INV deribvative of FAZ
!	
!	FB:     2nd integral, I2, multiplied by M/C1
!	FBZ:    ZETA-derivative of FB
!	FBM:    M_INV derivative of FB
!	FBZZ:   ZETA-derivative of FBZ
!	FBMM:   M_INV deribvative OF FBM
!	FBZM:   M_INV deribvative OF FBZ
!	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION M_INV
	DIMENSION VA(0:6),VB(0:6)
!	
!	Get coefficients for fixed M
!	
!	Calculate integral factors
!	
      CALL QPOLS_V (ZETA, VA, VB, FA, FAZ, FAZZ, FBX, FBZX, FBZZX)
!	
!	Get C-contribution
!	
      CALL QDIFC_V(ZETA, M_INV, F, FZ, FZZ)
!	
!	Restore values
!	
	VALC = 1.D0/F
	FB = FBX*VALC
	FBZ = VALC*(FBZX - FB*FZ)
	FBZZ = VALC*(FBZZX - FB*FZZ - 2.D0*FZ*FBZ)
	
	endsubroutine INTCON_V
!--------------------------------------------------------------------------------------------------
!
!
!
!	This routine calculates the derivatives of the power series
!	terms in PC-SAFT with respect to ZETA and 1/M

	SUBROUTINE QPOLS (M, ZETA, FA, FAM, FAMM, FAZ, FAZZ, FAZM,  FB, FBM, FBMM, FBZ, FBZZ, FBZM)

!	Input parameters:
!
!	M:      segment number
!	ZETA:   Volume fraction parameter
!
!	Output parameters:

!	FA:     Value of 'A' integral
!	FAM:    Derivative of FA wrt. 1/M
!	FAMM:   Derivative of FAM wrt. 1/M
!	FAZ:    Derivative of FA wrt. ZETA
!	FAZZ:   2nd derivative wrt. ZETA
!	FAZM:   Derivative oF FAM wrt. ZETA
!	FB:     Value of 'B' integral
!	FBM:    Derivative of FB wrt. 1/M
!	FBMM:   Derivative of FBM wrt. 1/M
!	FBZ:    Derivative of FB wrt. ZETA
!	FBZZ:   2nd derivative wrt. ZETA
!	FBZM:   Derivative of FBM wrt. ZETAC

	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	
	PARAMETER (SCAFA1 = -12.D0, SCAFA2 = -6.D0)
	
	DOUBLE PRECISION M,M_INV
	DIMENSION F(0:6),DF(0:6),DF2(0:6)
	
	COMMON /DISPTM/ A0(0:6),A1(0:6),A2(0:6),B0(0:6),B1(0:6),B2(0:6)
	
	F(0) = 1.D0
	DF(0) = 0.D0
	DF2(0) = 0.D0

!	Polynomials and 1st and second derivatives
	M_INV = 1.D0/M
	F(1) = ZETA
	DF(1) = 1.D0
	DF2(1) = 0.D0
	DO I=2,6
	
	  X = DBLE(I)
	  F(I) = ZETA*F(I-1)
	  DF(I) = X*F(I-1)
	  DF2(I) = X*DF(I-1)
	  
	ENDDO
	FAC1 = 1.D0 - M_INV
	FAC2 = 1.D0 - 2.D0*M_INV
	FTOT = FAC1*FAC2
	FAC1_M = -1.D0
	FAC2_M = -2.D0
	FTOT_M = FAC1_M*FAC2 + FAC1*FAC2_M
	FTOT_MM = 4.D0

!	Get sums, first integral
!
!	Do all sums for A first, then B
	S0A0 = A0(0)  + A0(1)*F(1)+A0(2)*F(2)+A0(3)*F(3)+A0(4)*F(4)   &
	              + A0(5)*F(5)+A0(6)*F(6)
	S0A1 = A1(0)  + A1(1)*F(1)+A1(2)*F(2)+A1(3)*F(3)+A1(4)*F(4)   &
	              + A1(5)*F(5)+A1(6)*F(6)
	S0A2 = A2(0)  + A2(1)*F(1)+A2(2)*F(2)+A2(3)*F(3)+A2(4)*F(4)   &
	              + A2(5)*F(5)+A2(6)*F(6)
	S1A0 = A0(1)  + A0(2)*DF(2)+A0(3)*DF(3)+A0(4)*DF(4)           &
	              + A0(5)*DF(5)+A0(6)*DF(6)
	S1A1 = A1(1)  + A1(2)*DF(2)+A1(3)*DF(3)+A1(4)*DF(4)           &
	              + A1(5)*DF(5)+A1(6)*DF(6)
	S1A2 = A2(1)  + A2(2)*DF(2)+A2(3)*DF(3)+A2(4)*DF(4)           &
	              + A2(5)*DF(5)+A2(6)*DF(6)
	S2A0 = 2.D0*A0(2) + DF2(3)*A0(3) + DF2(4)*A0(4)               &
	                  + DF2(5)*A0(5) + DF2(6)*A0(6)
	S2A1 = 2.D0*A1(2) + DF2(3)*A1(3) + DF2(4)*A1(4)               &
	                  + DF2(5)*A1(5) + DF2(6)*A1(6)
	S2A2 = 2.D0*A2(2) + DF2(3)*A2(3) + DF2(4)*A2(4)               &
	                  + DF2(5)*A2(5) + DF2(6)*A2(6)
	
	FA   = S0A0 + FAC1*S0A1 + FTOT*S0A2
	FAZ  = S1A0 + FAC1*S1A1 + FTOT*S1A2
	FAZZ = S2A0 + FAC1*S2A1 + FTOT*S2A2
	
	FAM  =  FAC1_M*S0A1 + FTOT_M*S0A2
	FAMM =              + FTOT_MM*S0A2
	FAZM =  FAC1_M*S1A1 + FTOT_M*S1A2
	
!	Add 'scale factors'
	FA = FA * SCAFA1
	FAZ = FAZ * SCAFA1
	FAZZ = FAZZ * SCAFA1
	FAM = FAM * SCAFA1
	FAMM = FAMM * SCAFA1
	FAZM = FAZM * SCAFA1
	
	
!	Get sums, second integral
	S0B0 = B0(0)  + B0(1)*F(1)+B0(2)*F(2)+B0(3)*F(3)+B0(4)*F(4)  &
	              + B0(5)*F(5)+B0(6)*F(6)
	S0B1 = B1(0)  + B1(1)*F(1)+B1(2)*F(2)+B1(3)*F(3)+B1(4)*F(4)  &
	              + B1(5)*F(5)+B1(6)*F(6)
	S0B2 = B2(0)  + B2(1)*F(1)+B2(2)*F(2)+B2(3)*F(3)+B2(4)*F(4)  &
	              + B2(5)*F(5)+B2(6)*F(6)
	S1B0 = B0(1)  + B0(2)*DF(2)+B0(3)*DF(3)+B0(4)*DF(4)          &
	              + B0(5)*DF(5)+B0(6)*DF(6)
	S1B1 = B1(1)  + B1(2)*DF(2)+B1(3)*DF(3)+B1(4)*DF(4)          &
	              + B1(5)*DF(5)+B1(6)*DF(6)
	S1B2 = B2(1)  + B2(2)*DF(2)+B2(3)*DF(3)+B2(4)*DF(4)          &
	              + B2(5)*DF(5)+B2(6)*DF(6)
	S2B0 = 2.D0*B0(2) + DF2(3)*B0(3) + DF2(4)*B0(4)              &
	                  + DF2(5)*B0(5) + DF2(6)*B0(6)
	S2B1 = 2.D0*B1(2) + DF2(3)*B1(3) + DF2(4)*B1(4)              &
	                  + DF2(5)*B1(5) + DF2(6)*B1(6)
	S2B2 = 2.D0*B2(2) + DF2(3)*B2(3) + DF2(4)*B2(4)              &
	                  + DF2(5)*B2(5) + DF2(6)*B2(6)
	
	FB   = S0B0 + FAC1*S0B1 + FTOT*S0B2
	FBZ  = S1B0 + FAC1*S1B1 + FTOT*S1B2
	FBZZ = S2B0 + FAC1*S2B1 + FTOT*S2B2
	
	FBM  =  FAC1_M*S0B1 + FTOT_M*S0B2
	FBMM =                FTOT_MM*S0B2
	FBZM =  FAC1_M*S1B1 + FTOT_M*S1B2
	
	
	
	FB = FB * SCAFA2
	FBZ = FBZ * SCAFA2
	FBZZ = FBZZ * SCAFA2
	FBM = FBM * SCAFA2
	FBMM = FBMM * SCAFA2
	FBZM = FBZM * SCAFA2
	
	endsubroutine QPOLS
!--------------------------------------------------------------------------------------------------
!
!
!	
      SUBROUTINE QPOLS_V (ZETA, AV, BV, FA, FAZ, FAZZ, FB, FBZ, FBZZ)

!	INPUT PARAMETERS:

!	ZETA:   Volume fraction parameter
!	AV:     1st coefficient array
!	BV:     2nd coeffcient array

!	Output parameters:

!	FA:     Value of 'A' integral
!	FAZ:    Derivative of FA wrt. ZETA
!	FAZZ:   2nd derivative wrt. ZETA
!	FB:     Value of 'B' integral
!	FBZ:    Derivative of FB wrt. ZETA
!	FBZZ:   2nd derivative wrt. ZETA

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION AV(0:6),BV(0:6)
	DIMENSION F(0:6),DF(0:6),DF2(0:6)
	F(0) = 1.D0
	DF(0) = 0.D0
	DF2(0) = 0.D0
	
!	POLYNOMIALS AND 1ST AND SECOND DERIVATIVES
	
	DO I=1, 6
	
	  X = DBLE(I)
	  F(I) = ZETA*F(I-1)
	  DF(I) = X*F(I-1)
	  DF2(I) = X*DF(I-1)
	  
	ENDDO
	
!	GET SUMS, FIRST INTEGRAL
	S0A0 = AV(0)  + AV(1)*F(1)+AV(2)*F(2)+AV(3)*F(3)+AV(4)*F(4)    &
	              + AV(5)*F(5)+AV(6)*F(6)
	S1A0 = AV(1)  + AV(2)*DF(2)+AV(3)*DF(3)+AV(4)*DF(4)            &
	              + AV(5)*DF(5)+AV(6)*DF(6)
	S2A0 = 2.D0*AV(2) + DF2(3)*AV(3) + DF2(4)*AV(4)                &
	                  + DF2(5)*AV(5) + DF2(6)*AV(6)
	
	FA   = S0A0
	FAZ  = S1A0
	FAZZ = S2A0
	
!	GET SUMS, SECOND INTEGRAL
	
	
	S0B0 = BV(0)  + BV(1)*F(1)+BV(2)*F(2)+BV(3)*F(3)+BV(4)*F(4)    &
	              + BV(5)*F(5)+BV(6)*F(6)
	S1B0 = BV(1)  + BV(2)*DF(2)+BV(3)*DF(3)+BV(4)*DF(4)            &
	              + BV(5)*DF(5)+BV(6)*DF(6)
	S2B0 = 2.D0*BV(2) + DF2(3)*BV(3) + DF2(4)*BV(4)                &
	                  + DF2(5)*BV(5) + DF2(6)*BV(6)
	
	FB   = S0B0
	FBZ  = S1B0
	FBZZ = S2B0
	endsubroutine QPOLS_V
	
!--------------------------------------------------------------------------------------------------
!
!
!	
!	THIS CALCULATES THE ZETA AND 1/M-DERIVATIVVES OF C1/M

	SUBROUTINE QDIFC (ZETA, M, F, FZ, FM, FZZ, FZM)

!	Input:

!	M:      Segment number
!	ZETA:   Volume fraction parameter
!
!	Output:
!
!	F:      C1-factor
!	FM:     M_INV-derivative of F
!	FZ:     ZETA-derivative of F
!	FZZ:    ZETA-derivative of FZ
!	FZM:    ZETA-derivative of FM
!
!
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION M_INV,M
	M_INV = 1.D0/M
	
!	Powers to be used subsequently
	X =1.D0/(2.D0-ZETA)
	Y =1.D0/(1.D0-ZETA) 
	X2 = X*X
	Y2 = Y*Y
	X3 = X2*X
	Y3 = Y2*Y
	X4 = X3*X
	Y4 = Y3*Y
	Y5 = Y3*Y2
	Y6 = Y3*Y3
	FM = -1.D0 - 4.D0*X2 + 3.D0*Y2
	F0 = 2.D0 - 5.D0*Y2 -4.D0*Y3 + 6.D0*Y4 + 4.D0*X2
	F = F0 + FM*M_INV
	FZM = -8.D0*X3 + 6.D0*Y3
	F0Z = -10.D0*Y3 - 12.D0*Y4 +24.D0*Y5 + 8.D0*X3
	FZ = F0Z + M_INV*FZM
	F0ZZ = -30.D0*Y4 - 48.D0*Y5 + 120.D0*Y6 + 24.D0*X4
	FZZ = F0ZZ + M_INV*(-24.D0*X4 + 18.D0*Y4)
	endsubroutine QDIFC
!--------------------------------------------------------------------------------------------------
!
!
!
!	This calculates the ZETA of C1/M

      SUBROUTINE QDIFC_V (ZETA, M_INV, F, FZ, FZZ)

!	Input:

!	M:      Segment number
!	ZETA:   Volume fraction parameter
!
!	OUTPUT:
!
!	F:      C1-factor
!	FM:     M_INV-derivative of F
!	FZ:     ZETA-derivative of F
!	FZZ:    ZETA-derivative of FZ
!	FZM:    ZETA-derivative of FM
!
!
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DOUBLE PRECISION M_INV,M
!
!	Powers to be used subsequently
	X =1.D0/(2.D0-ZETA)
	Y =1.D0/(1.D0-ZETA)
	X2 = X*X
	Y2 = Y*Y
	X3 = X2*X
	Y3 = Y2*Y
	X4 = X3*X
	Y4 = Y3*Y
	Y5 = Y3*Y2
	Y6 = Y3*Y3
	FM = -1.D0 - 4.D0*X2 + 3.D0*Y2
	F0 = 2.D0 - 5.D0*Y2 -4.D0*Y3 + 6.D0*Y4 + 4.D0*X2
	F = F0 + FM*M_INV
	FZM = -8.D0*X3 + 6.D0*Y3
	F0Z = -10.D0*Y3 - 12.D0*Y4 +24.D0*Y5 + 8.D0*X3
	FZ = F0Z + M_INV*FZM
	F0ZZ = -30.D0*Y4 - 48.D0*Y5 + 120.D0*Y6 + 24.D0*X4
	FZZ = F0ZZ + M_INV*(-24.D0*X4 + 18.D0*Y4)
	endsubroutine QDIFC_V
!--------------------------------------------------------------------------------------------------
!
!
!	Polynomial coefficients of the integrals I1 and I2
!
	BLOCK DATA DPARAM
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	COMMON /DISPTM/ A0(0:6),A1(0:6),A2(0:6),B0(0:6),B1(0:6),B2(0:6)
	DATA A0/0.9105631445D0, 0.6361281449D0, 2.6861347891D0, -26.547362491D0,97.759208784D0,-159.59154087D0,91.297774084D0/
	DATA A1/-0.3084016918D0, 0.1860531159D0,-2.5030047259D0, 21.419793629D0,-65.255885330D0, 83.318680481D0,-33.74692293D0/
	DATA A2/-0.0906148351D0, 0.4527842806D0, 0.5962700728D0, -1.7241829131D0,-4.1302112531D0, 13.776631870D0,-8.6728470368D0/
	DATA B0/0.7240946941D0, 2.2382791861D0, -4.0025849485D0, -21.003576815D0, 26.855641363D0, 206.55133841D0,-355.60235612D0/
	DATA B1/-0.5755498075D0, 0.6995095521D0, 3.892567339D0, -17.215471648D0, 192.67226447D0, -161.82646165D0, -165.20769346D0/
	DATA B2/0.0976883116D0, -0.2557574982D0, -9.155856153D0, 20.642075974D0, -38.804430052D0, 93.626774077D0, -29.666905585D0/
	ENDblockdata DPARAM
!--------------------------------------------------------------------------------------------------
!
!
	SUBROUTINE XSUMS (NC, X, AR1, AR2, SS1, SS2, XD1, XD2)
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (NCM = 30)
	DIMENSION AR1(NCM,NCM),AR2(NCM,NCM),XD1(NC),XD2(NC),X(NC)
	
	SS1 = 0.D0
	SS2 = 0.D0
	DO I = 1, NC
	
	  XD1I = 0.D0
	  XD2I = 0.D0
	  DO J = 1, NC
	  
	    XD1I = XD1I + X(J)*AR1(J,I)
	    XD2I = XD2I + X(J)*AR2(J,I)
	    
	  ENDDO
	  XD1(I) = XD1I
	  XD2(I) = XD2I
	  SS1 = SS1 + X(I)*XD1I
	  SS2 = SS2 + X(I)*XD2I
	  
	ENDDO
	SS1 = 0.5D0*SS1
	SS2 = 0.5D0*SS2
	endsubroutine XSUMS
!--------------------------------------------------------------------------------------------------
!
!	
	SUBROUTINE VATTR (NC, NDER, NTEMP, X, T, DD, DDT, S, ST, STT, &
	                  SS1, SS2, XD1, XD2, AR1, AR2, V,            &
	                  QA, QAV, QAVV, QAT, QAVT, QATT, QAN, QAVN, QAN2, QATN)
	
!	This subroutine calculates the 'attractive' contribution to the
!	reduced, residual Helmholz energy, and its derivatives wrt.
!	T, V and N

!	Input:

!	NC:          No. of components in mixture
!	NDER:        Level of composition derivatives (1 OR 2)
!	NTEMP:       Level of temperature derivatives (0,1 OR 2)
!	X:           Composition vector
!	T:           Temperature
!	DD:          M*DIA^K array , K=0,1,2,3
!	DDT:         1/T-derivatives of DD
!	S:           Sums OF DD*X
!	ST:          1/T-derivatives of S
!	STT:         1/T-derivatives of ST
!	SS1:         1st quadratic energy sum
!	SS2:         2nd quadratic energy sum
!	XD1:         SS1 N-derivave vector
!	XD2:         SS2 N-derivave vector
!	AR1:         XD1 N-derivative matrix
!	AR2:         XD2 N-derivative matrix
!	V:           Volume (/R)

!	Output:

!	QA           A^res/RT contribution
!	QAV:         V-deriv. of QA
!	QAVV:        V-deriv. of QAV
!	QAT:         1/T-deriv. of QA          NTEMP >= 1
!	QATV:        1/T-deriv. of QAV         NTEMP >= 1
!	QATT:        1/T-deriv. of QAT         NTEMP >= 2
!	QAN:         N-deriv. of QA
!	QAN2:        N-deriv. of QAN           NDER >= 2
!	QAVN:        N-deriv. of QAV
!	QATN:        N-deriv. of QAT           NTEMP >=1
	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	PARAMETER (NCM = 30)
	DIMENSION X(NC),DD(0:3,NCM),DDT(0:3,NCM),S(0:3),ST(0:3)
	DIMENSION STT(0:3),QAN(NC),QAVN(NC)
	DIMENSION QAN2(NC,NC),QATN(NC)
	DIMENSION AR1(NCM,NCM),AR2(NCM,NCM)
	DIMENSION QBN(NC),QBVN(NC),QBN2(NC,NC),QBTN(NC)
	DOUBLE PRECISION M,M_INV,MRD
	DIMENSION MRD(NC)
	DIMENSION XD1(NC),XD2(NC)

!	Get the quadratics	
	SUMX= 0.D0
	DO I=1,NC
	
	  SUMX = SUMX + X(I)
	  
	ENDDO
	SUMS3 = S(3)
	SUMMX = S(0)
	SUMXR = 1.D0/SUMX
	M = SUMMX*SUMXR
	SUMMXR = 1.D0/SUMMX
	M_INV = SUMX*SUMMXR
	DO I=1,NC
	
	  MRD(I) = (1.D0 - M_INV*DD(0,I))*SUMMXR
	  
	ENDDO
	VR = 1.D0/V
	ZETA = SUMS3*VR
	ZETAV = -ZETA*VR
	ZETAVV = -2.D0*ZETAV*VR
	VR2 = VR*VR
	
!	Basic derivatives
	CALL INTCON (M, ZETA, QA, FAZ, FAM, FAZZ, FAMM, FAZM, QB, FBZ, FBM, FBZZ, FBMM, FBZM)

!	Chain rule:
!	1st der.      F_X  = F_M M_X + F_A A_X
!	2nd der.      F_XY = F_M M_XY + F_MM M_X M_Y + F_MA M_X A_Y
!	                  +  F_A A_XY + F_AM A_X M_Y + F_AA A_X A_Y

	FAZV = FAZ*VR
	FBZV = FBZ*VR
	QAV = FAZ*ZETAV
	QBV = FBZ*ZETAV
	
	ZETAT = VR*ST(3)
	ZETATT =VR*STT(3)
	ZETATV = -VR*ZETAT
	QAT = FAZ*ZETAT
	QBT = FBZ*ZETAT

!	1/T-derivatives of DQ/DV
	QAVT = FAZ*ZETATV + FAZZ*ZETAV*ZETAT
	QBVT = FBZ*ZETATV + FBZZ*ZETAV*ZETAT
	QATT = FAZ*ZETATT + FAZZ*ZETAT**2
	QBTT = FBZ*ZETATT + FBZZ*ZETAT**2
	
	QAVV= FAZZ*ZETAV**2 + FAZ*ZETAVV
	QBVV= FBZZ*ZETAV**2 + FBZ*ZETAVV
	
!	Coefficients for V-cross derivatives
	CNVA1 = FAZM*ZETAV
	CNVA2 = VR*(ZETAV*FAZZ - VR*FAZ)
	CNVB1 = FBZM*ZETAV
	CNVB2 = VR*(ZETAV*FBZZ - VR*FBZ)
	FAZVT = VR*(FAZZ*ZETAT)
	FAMT = FAZM*ZETAT
	FBZVT = VR*(FBZZ*ZETAT)
	FBMT = FBZM*ZETAT
	
	DO I=1,NC
	
	  QAN(I)  = FAZV*DD(3,I) + FAM*MRD(I)
	  QAVN(I) = CNVA1*MRD(I) + CNVA2*DD(3,I)
	  QBN(I)  = FBZV*DD(3,I) + FBM*MRD(I)
	  QBVN(I) = CNVB1*MRD(I) + CNVB2*DD(3,I)
	  
	ENDDO
	
	IF (NTEMP > 0) THEN
	  DO I = 1,NC
	  
	    QATN(I) = FAZVT*DD(3,I) + FAMT*MRD(I) + FAZV*DDT(3,I)
	    QBTN(I) = FBZVT*DD(3,I) + FBMT*MRD(I) + FBZV*DDT(3,I)
	    
	  ENDDO
	ENDIF
	
	IF (NDER > 1) THEN
	
!	  The double loop
	  CON1A= FAMM + 2.D0*M*FAM
	  CON2A= -FAM*SUMXR
	  CON3A= VR*FAZM
	  CON4A= FAZZ*VR2
	  CON1B= FBMM + 2.D0*M*FBM
	  CON2B= -FBM*SUMXR
	  CON3B= VR*FBZM
	  CON4B= FBZZ*VR2
    
!	  Compact results
	  DO I=1,NC
	  
	    FACG_JA = CON1A*MRD(I) + CON3A*DD(3,I) + CON2A
	    FACD_JA = CON3A*MRD(I) + CON4A*DD(3,I)
	    ADDIA =   CON2A*MRD(I)
	    FACG_JB = CON1B*MRD(I) + CON3B*DD(3,I) + CON2B
	    FACD_JB = CON3B*MRD(I) + CON4B*DD(3,I)
	    ADDIB =   CON2B*MRD(I)
	    DO J=1,I
	    
	      QAN2(J,I) = FACG_JA*MRD(J) + FACD_JA*DD(3,J) + ADDIA
	      QBN2(J,I) = FACG_JB*MRD(J) + FACD_JB*DD(3,J) + ADDIB
	      
	    ENDDO
	    
	  ENDDO
!
!	  Update with sums
	  DO I=1,NC
	    DO J=1,I
	    
	      QAN2(J,I) = SS1*QAN2(J,I) + QA*AR1(J,I) + QAN(J)*XD1(I) + QAN(I)*XD1(J)
	
	    ENDDO
	  ENDDO
	  
	ENDIF
	DO I=1,NC
	
	  QAVN(I) = QAVN(I)*SS1 + QAV*XD1(I)
	  
	ENDDO
!
!	For T-derivatives, SS1 and XD1 are proportional to 1/T,
!	                   SS2 and XD2 TO 1/T^2
	T2 = T + T
	DO I=1,NC
	  QAN(I) = QAN(I)*SS1 + QA*XD1(I)
	ENDDO
	IF (NTEMP > 0) THEN
	  DO I=1,NC
	  
	    QATN(I) = SS1*QATN(I) + QAT*XD1(I) + T*QAN(I)
	    
	  ENDDO
	ENDIF
	QAVV = QAVV * SS1
	QAV  = QAV  * SS1
	QAVT = QAVT*SS1 + QAV*T
	QATT = SS1*(QATT + T2*QAT)
	QA   = QA   * SS1
	QAT = QAT*SS1 + QA*T
	
!	Update with sums for second term and add scale factor, 1/V
	IF (NDER > 1) THEN
	  DO I=1,NC
	    DO J=1,I
	    
	      QAN2(J,I) = VR*(QAN2(J,I) + SS2*QBN2(J,I) + QB*AR2(J,I) + QBN(J)*XD2(I) + QBN(I)*XD2(J))
	    
	    ENDDO
	  ENDDO
	ENDIF
	DO I=1,NC
	
	  QBNI = QBN(I)*SS2 + QB*XD2(I)
	  QAN(I) =  VR*(QAN(I) + QBNI)
	  IF (NTEMP > 0) THEN
	  
	    QBTNI = QBTN(I)*SS2 + T2*QBNI + QBT*XD2(I)
	    QATN(I) = VR*(QATN(I) + QBTNI)
	    
	  ENDIF
	  
	ENDDO
	DO I=1,NC
	
	  QAVN(I) = VR*(QAVN(I) + QBVN(I)*SS2 + QBV*XD2(I)-QAN(I))
	  
	ENDDO
	
!	Incorporate 1/T-derivatives
	QAX    = QA  + QB   * SS2
	QAXT   = QAT + (QBT +T2*QB)*SS2
	SS2T = T2*SS2
	SS2TT = T*SS2T
	QAVXT = QAVT + (QBVT  + QBV*T2)*SS2
	QAXTT  = QATT + QBTT*SS2 + 2.D0*SS2T*QBT + QB*SS2TT
	QATT = VR*QAXTT
	QAVVX = QAVV + QBVV * SS2
	QAVX  = QAV  + QBV  * SS2
	QA = VR*QAX
	QAV = VR*(QAVX - QA)
	QAVV = VR*(QAVVX - 2.D0*QAV)
	QAT = VR*QAXT
	QAVT = VR*(QAVXT - QAT)
	endsubroutine VATTR
!--------------------------------------------------------------------------------------------------
!
!
	SUBROUTINE CARNAH (NC, NDER, NTEMP, X, DD, DDT, DTT, S, ST, ST2, V, &
	                   Q, QV, QVV, QT, QTV, QTT, DQDN, DQ2DN2, DQNV, DQDNT)
	
!	The subroutine calculates the repulsive contribution to the
!	reduced, residual Helmholz energy, and its derivatives wrt. V,
!	T AND N of the "Simplified PC-SAFT".
!	"sPC-SAFT" uses Carnahan-Starling free-volume term, instead 
!	of Boublík-Mansoori-Leeland.

!	Input:
!	NC:          No. of components in mixture
!	NDER:        Level of composition derivatives (1 or 2)
!	NTEMP:       Level of temperature derivatives (0,1 or 2)
!	X:           Composition vector
!	DD:          M*DIA^K array , K=0,1,2,3
!	DDT:         1/T-derivatives of DD
!	DTT:         1/T-derivatives of DDT
!	S:           Sums of DD*X
!	ST:          1/T-Derivatives of S
!	ST2:         1/T-Derivatives of ST
!	V:           Volume (/R)
!
!	Output:
!
!	Q:           A^res/RT contribution
!	QV:          V-deriv. of Q
!	QVV:         V-deriv. of QV
!	QT:          1/T-deriv. of Q          NTEMP >= 1
!	QTV:         1/T-deriv. of QV         NTEMP >= 1
!	QTT:         1/T-deriv. of QT         NTEMP >= 2
!	DQDN         N-deriv. of Q
!	DQ2DN2       N-deriv. of QN           NDER >= 2
!	DQNV         N-deriv. of QV
!	DQDNT        N-deriv. of QT           NTEMP >=1

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	
	PARAMETER (NCM = 30)
	DIMENSION X(NC),DD(0:3,NCM),DDT(0:3,NCM),DTT(0:3,NCM)
	DIMENSION DQDN(NC), DQNV(NC), DQ2DN2(NC,NC),DQDNT(NC)
	DIMENSION S(0:3),ST(0:3),ST2(0:3)
	DIMENSION CUFG(NC)
	DOUBLE PRECISION M
	
	DO K = 0, 3
	
	  ST(K) = 0.D0
	  ST2(K) = 0.D0
	  DO I=1,NC
	  
	    ST(K) = ST(K) + X(I)*DDT(K,I)
	    ST2(K) = ST2(K) + X(I)*DTT(K,I)
	    
	  ENDDO
	  
	ENDDO
	M = S(0)
	VR = 1.D0/V
	ZETA = S(3)*VR
	FINV1 = 1.D0/(1.D0-ZETA)
	ZT2 = 2.D0 - ZETA
	FINV2 = 1.D0/(2.D0-ZETA)
	
!	Powers of X
	P1 = FINV1
	P2 = FINV1*FINV1
	P3 = FINV1*P2
	
!	Functions    
	G = -3.D0 + 2.D0*P1 + P2 	!Carnahan-Starling A^hs/RT (repulsive without chain)
	P4 = FINV1 * P3
	F = dLOG (0.5D0*ZT2*P3)		!A^chain / [RT(1-m)] = ln(g^hs)
	H = G - F
	GZ = 2.D0*(P2 + P3)
	GZZ = 4.D0*P3 + 6.D0*P4
	FZ = - FINV2 + 3.D0*P1
	FZZ = -FINV2**2 + 3.D0*P2
	HZ = GZ - FZ
	HZZ = GZZ - FZZ
	ZETAV = -ZETA*VR
	ZETAVV =-(ZETAV+ZETAV)*VR
	Q = M*H + F				!A^r/RT (sphere + chain formantion)
	ZD1 = M*HZ + FZ
	QV= ZD1*ZETAV
	QZZ = M*HZZ+FZZ
	QVV= ZD1*ZETAVV + QZZ*ZETAV**2
	ZETAT = ST(3)*VR
	ZETATT = ST2(3)*VR
	ZETAVT = -VR*ZETAT
	QT = ZD1*ZETAT
	QTT= ZD1*ZETATT + QZZ*ZETAT**2
	QTV= ZD1*ZETAVT + QZZ*ZETAT*ZETAV
	
!	Composition derivatives
	QZ = M*HZ + FZ
	QZ_VR = VR*QZ
	DO I = 1,NC
	
	  DQDN(I) = DD(0,I)*H + F + QZ_VR*DD(3,I)
	  
	ENDDO
	HV = HZ * ZETAV
	FV = FZ * ZETAV
	HVV2 = (M*HZZ+FZZ)*ZETAV
	
!	Volume cross derivatives    
	FD3I = (HVV2 - QZ_VR)*VR
	DO I = 1, NC
	
	  DQNV(I) = HV*DD(0,I) + FV + FD3I*DD(3,I)
	  
	ENDDO
	
!	Temperature cross derivatives	
	IF (NTEMP > 0) THEN
	
	  HT = HZ * ZETAT
	  FT = FZ * ZETAT
	  HTT2R= (M*HZZ+FZZ)*ZETAT*VR
	  DO I=1,NC
	  
	    DQDNT(I) =HT*DD(0,I) + FT + QZ_VR*DDT(3,I) + HTT2R*DD(3,I)
	    
	  ENDDO
	  
	ENDIF
	IF (NDER >= 2) THEN
	
	  CONA = VR*(M*HZZ+FZZ)*VR
	  DO I=1,NC
	  
	    CUFG(I) = VR*(DD(0,I)*HZ+FZ)
	    ADDI = CUFG(I)+CONA*DD(3,I)
	    DO J=1,I
	    
	      DQ2DN2(J,I) = ADDI*DD(3,J)+DD(3,I)*CUFG(J)
		
	    ENDDO
	    
	  ENDDO
	  
	ENDIF
	endsubroutine CARNAH
!--------------------------------------------------------------------------------------------------
!
!
	SUBROUTINE GZL_AVC (NC, NDER, NTEMP, DD, DDT, X, S, ST, ST2, V, &
	                    Q, QV, QVV, QT, QTV, QTT, DQDN, DQ2DN2, DQNV, DQNT)

!	This subroutine calculates the attractive chain-term and its
!	1st and second derivatives wrt T, V and N, for the sPC-SAFT. 

!	Input:

!	NC:          No. of components
!	NDER:        Level of comp.derivatives
!	NTEMP:       Level of temp.derivatives
!	DD:          Power of diameter array
!	DDT:         Derivative of DD wrt 1./T
!	X:           Molar feed
!	S:           Diameter power sum
!	ST:          1/T-derivative of S
!	ST2:         1/T-derivative of ST
!	V:           Volume
!
!	OUtput
!
!	Q:           Residual A- contribution
!	QV:          Derivative of Q wrt V
!	QVV:         Derivative of QV wrt V
!	QT:          Derivative of Q wrt 1/T
!	QTV:         Derivative of QV wrt 1/T
!	QTT:         Derivative of QT wrt 1/T
!
!	DQDN:        N-derivatives of Q
!	D2QDN2:      N-derivatives of DQDN
!	DQNV:        N-derivatives of QV
!	DQNT:        N-derivatives of QT
!
!	We use here the modified expression
!
!	SUM_I X(I)(M(I)-1)ln g, 
!
!	with g = (1-.5*ZETA)/(1-ZETA)^3, ZETA = S(3)/V
!	
!	being "g" the hard-sphere radial distribution function at contact
!	from Carnahan-Starling.
!
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	
	PARAMETER (NCM = 30)
	
	DIMENSION DD(0:3,NCM),DDT(0:3,NCM)
	DIMENSION X(NC),S(0:3),ST(0:3),ST2(0:3)
	DIMENSION DQDN(NC), DQNV(NC)
	DIMENSION DQ2DN2(NC,NC),DQNT(NC)
	DIMENSION XFACT(NC)
	common /assoc/ NST	
	common /RDF_sPCSAFT/ GZ, GZZ, ZM2, zeta, zetaT, zetaTT, zetaTV, zetaV, zetaVV
	
	
!	Basic factors
	VR = 1.D0/V
	VR2 = VR*VR
	ZETA = S(3)*VR
	XMSUM = 0.D0
	DO I=1,NC
	
	  XFACT(I) = DD(0,I) - 1.D0 !m(i) - 1
	  XMSUM = XMSUM + X(I)*XFACT(I)
	  
	ENDDO
	ZETINV = 1.D0/(1.D0-ZETA) !Y
	ZM2 = 2.D0 - ZETA
	ZETF2= 1.D0 / ZM2

!	Radial distribution function
	RDF = 0.5D0*ZM2 *ZETINV**3
	
	G = dlog(RDF)
! 	G = dLOG (0.5D0*ZM2 *ZETINV**3)
	Q = XMSUM*G
!
!	G-derivatives
	GZ   = -ZETF2 + 3.D0*ZETINV
	GZZ  =3.D0*ZETINV**2 -ZETF2*ZETF2
	Q_M = G
	Q_Z = XMSUM*GZ
	Q_ZZ = XMSUM*GZZ
	Q_ZM = GZ
	ZETAV = -VR*ZETA
	ZETAVV = -2.D0*VR*ZETAV
	ZETAT = ST(3)*VR
	ZETATV = -VR*ZETAT
	ZETATT =ST2(3)*VR
	QV = Q_Z*ZETAV
	QT = Q_Z*ZETAT
	QVV = Q_Z*ZETAVV + Q_ZZ*ZETAV**2
	QTV = Q_Z*ZETATV + Q_ZZ*ZETAV*ZETAT
	QTT = Q_Z*ZETATT + Q_ZZ*ZETAT*ZETAT
	QZF1 = VR*Q_Z
	QZMZ = Q_ZM*ZETAV
	QZV2 = Q_ZZ*ZETAV*VR - Q_Z*VR2
	
	DO  I= 1,NC
	
	  DQDN(I) = Q_M*XFACT(I) + QZF1*DD(3,I)      !DD(3,i) = pi*di^3/(k6)
	  DQNV(I) = QZMZ*XFACT(I) + QZV2*DD(3,I)
	  
	ENDDO
	IF (NTEMP > 0) THEN
	
	  QZT1 = Q_ZM*ZETAT
	  QZT2 = Q_ZZ*ZETAT*VR
	  QZT3 = Q_Z*VR    
	  DO I=1,NC
	  
	    DQNT(I) = QZT1*XFACT(I)+QZT2*DD(3,I)+QZT3*DDT(3,I)
	    
	  ENDDO
	  
	ENDIF
	IF (NDER > 1) THEN
	
	  QFF1=Q_ZM*VR
	  QFF2 = Q_ZZ*VR2
	  DO I=1,NC
	  
	    QFI1X= QFF1*XFACT(I)
	    QFI1D= QFF1*DD(3,I)
	    QF2I = QFF2*DD(3,I)+QFI1X
	    DO J=1,I
	    
	      DQ2DN2(J,I) =  QFI1D*XFACT(J)+QF2I*DD(3,J)
	    
	    ENDDO
	    
	  ENDDO
	   
	ENDIF

! 	if (NST > 0) then
! 
! ! 	  Variables used for RDF in association strenght
! 	  do i = 1, NC
! 	    
! 	    common_d(i) = dia(i)
! 	    common_dd(0:3,i) = DD(0:3,i)
! 	    common_ddT(0:3,i) = DDT(0:3,i)
! 	    common_S(0:3) = S(0:3)
! 	    common_ST(0:3) = ST(0:3)
! 	    common_ST2(0:3) = ST2(0:3)
! 	    
! 	  enddo	  
! 	
! 	endif
	endsubroutine GZL_AVC
!--------------------------------------------------------------------------------------------------
!
!
	SUBROUTINE HARDCN (NTYP, NC, NDER, NTEMP, X, DIA, DD, DDT, DTT, S, ST, ST2, V,  &
	                   Q, QV, QVV, QT, QTV, QTT, DQDN, DQ2DN2, DQNV, DQDNT)
!
!	The subroutine calculates the repulsive contribution to the
!	reduced, residual Helmholz energy, and its derivatives wrt. V,
!	T and N
!
!
!	Input:
!
!	NTYP:        Submodel, 1 = CS in G, 2 = org. PC-SAFT
!	NC:          No. of components in mixture
!	NDER:        Level of composition derivatives (1 OR 2)
!	NTEMP:       Level of temperature derivatives (0,1 OR 2)
!	X:           Composition vector
!	DIA:         Diameters
!	DD:          M*DIA^K array , K=0,1,2,3
!	DDT:         1/T-derivatives of DD
!	DTT:         1/T-derivatives of DDT
!	S:           Sums OF DD*X
!	ST:          1/T-derivatives of S
!	ST2:         1/T-derivatives of ST
!	V:           Volume (/R)
!
!	Output:
!
!	Q:           A^res/RT contribution
!	QV:          V-deriv. of Q
!	QVV:         V-deriv. of QV
!	QT:          1/T-deriv. of Q          NTEMP >= 1
!	QTV:         1/T-deriv. of QV         NTEMP >= 1
!	QTT:         1/T-deriv. of QT         NTEMP >= 2
!	DQDN         N-Deriv. of Q
!	DQ2DN2       N-Deriv. of QN           NDER >= 2
!	DQNV         N-Deriv. of QV
!	DQDNT        N-Deriv. of QT           NTEMP >=1

	IMPLICIT DOUBLE PRECISION (A-H, O-Z)

	PARAMETER (NCM = 30)
	
	DIMENSION X(NC),DD(0:3,NCM),DDT(0:3,NCM),DIA(NC),DTT(0:3,NCM)
	DIMENSION DQDN(NC), DQNV(NC), DQ2DN2(NC,NC),DQDNT(NC)
	DIMENSION DQDN_R(NC), DQNV_R(NC), DQ2DN2_R(NC,NC)
	DIMENSION DQNT_R(NC)
	DIMENSION S(0:3),ST(0:3),ST2(0:3)
!
!	Calculate contriubtion to repulsive A
!
!	Initialize sums
	DO K = 0, 3
	
!	  S(K) = 0.D0
	  ST(K) = 0.D0
	  ST2(K) = 0.D0
	  DO I=1,NC
	  
!	    S(K) = S(K) + X(I)*DD(K,I)
	    ST(K) = ST(K) + X(I)*DDT(K,I)
	    ST2(K) = ST2(K) + X(I)*DTT(K,I)
	    
	  ENDDO
	  
	ENDDO
	ST (0) = 0.D0
	ST2(0) = 0.D0

!	Hard sphere
	CALL HARDSP (NC, NDER, NTEMP, DD, DDT, S, ST, ST2, V,                &
	             Q, QV, QVV, QT, QTV, QTT, DQDN, DQ2DN2, DQNV, DQDNT)

!	Distriubtion function .... discount edition !
	IF (NTYP == 1) THEN
	
	  CALL GZL_AVC (NC, NDER, NTEMP, DD, DDT, X, S, ST, ST2, V,          &
	                Q_R, QV_R, QVV_R, QT_R, QTV_R, QTT_R, DQDN_R, DQ2DN2_R, DQNV_R, DQNT_R)
	  
	ELSE
	
	  CALL GZLCAL (NC, NDER, NTEMP, DD, DDT, DTT, X, DIA, S, ST, ST2, V, &
	               Q_R, QV_R, QVV_R, QT_R, QTV_R, QTT_R, DQDN_R, DQ2DN2_R, DQNV_R, DQNT_R)
	
	ENDIF
	
!	Subtract
	Q = Q - Q_R
	QV = QV - QV_R
	QVV = QVV - QVV_R
	QT = QT - QT_R
	QTV = QTV - QTV_R
	QTT = QTT - QTT_R
	DO I = 1, NC
	
	  DQDN (I) = DQDN (I) - DQDN_R(I)
	  DQNV (I) = DQNV (I) - DQNV_R(I)
	  
	ENDDO
	
	IF (NTEMP > 0) THEN
	
	  DO I=1, NC
	  
	    DQDNT(I) = DQDNT(I) - DQNT_R(I)
	    
	  ENDDO
	  
	ENDIF
	IF (NDER >= 2) THEN
	
	  DO I=1, NC
	  
	    DO J=1, I
	    
	      DQ2DN2(J,I) = DQ2DN2(J,I) - DQ2DN2_R(J,I)
	  
	    ENDDO
	    
	  ENDDO
	   
	ENDIF
	endsubroutine HARDCN
!--------------------------------------------------------------------------------------------------
!
!
	SUBROUTINE HARDSP (NC, NDER, NTEMP, DD, DDT, S, ST, ST2, V,    &
	                   Q, QV, QVV, QT, QTV, QTT, DQDN, DQ2DN2, DQNV, DQDNT)
!
!	Purpose: to calculate the hard-sphere contribution to the
!	Residual Helmholz energy, and its derivatives wrt V and N
!
!	Input:
!	NC:          No. of components
!	NDER:        Level of comp.derivatives
!	NTEMP:       Level of temp.derivatives
!	DD:          Diameter coefficient array
!	DDT:         1/T-derivative of DD
!	S:           Vector SUM (X(I)*D(K)(I)
!	ST:          1st 1/T-derivative of S
!	ST2:         1st 1/T derivative of ST
!	V:           Volume

!	Output:
!	Q:           Residual A/RT
!	QV:          Volume derivative of Q
!	QVV:         Volume derivative of QV
!	QT:          1/T-derivative of Q
!	QTV          V-derivative of QT
!	QT2          2ND 1/T-derivative of Q
!	DQDN:        N-derivative of Q
!	DQ2DN2:      N-derivative of DQDN
!	DQNV:        N-derivative of QV
!	DQDNT:       1/T-derivative of DQDN

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	
	PARAMETER (NCM = 30)
	DIMENSION DD(0:3, NCM), DDT(0:3, NCM)
	DIMENSION DQDN(NC),  DQNV(NC),  DQ2DN2(NC, NC)
	DIMENSION DQDNT(NC)
	DIMENSION S(0:3),  ST(0:3), QS(0:3), QVS(0:3), QSS(0:3, 0:3)
	DIMENSION QST(0:3), ST2(0:3)

!	Calculate derivatives with respect to the S-elements and V
	CALL QDERS (V, S, Q, QS, QVS, QSS, DQV, DQVV)

!	T-derivatives from chain rule
      QT = 0.D0
      QTT = 0.D0
      QTV = 0.D0
      DO K=1, 3
      
	  QT = QT + QS(K)*ST(K)
	  QTT = QTT + QS(K)*ST2(K)
	  QTV = QTV + QVS(K)*ST(K)
	  
	ENDDO
	QV  = DQV
	QVV = DQVV
!
!	Generate T-derivatives of the QS, using chainrule
	IF (NTEMP > 0) THEN
	
	  DO K=0,3
	  
	    QST(K) = 0.D0
	    DO L=1,3
	    
	      QST(K) = QST(K) + QSS(K,L)*ST(L)
		
	    ENDDO
!	    ...and add in contribution to the second derivative
	    QTT = QTT + QST(K)*ST(K)
	    
	  ENDDO
	  DO I=1,NC
	  
	    DQDNT(I) = QST(0)*DD(0,I)+QST(1)*DD(1,I)+QST(2)*DD(2,I) +QST(3)*DD(3,I) &
	              +QS(1)*DDT(1,I)+QS(2)*DDT(2,I)+ QS(3)*DDT(3,I)
	    
	  ENDDO
	
	ENDIF
!
!	Chain rule based comp-derivatives
	DO I=1,NC
	
	  DQDN(I) = QS(0)*DD(0,I)+QS(1)*DD(1,I)+QS(2)*DD(2,I)+ QS(3)*DD(3,I)
!	  DQDNT(I) = QST(0)*DD(0,I)+QST(1)*DD(1,I)+QST(2)*DD(2,I)+QST(3)*DD(3,I) &
!	            +QS(1)*DDT(1,I)+QS(2)*DDT(2,I)+ QS(3)*DDT(3,I)
	  DQNV(I) = QVS(0)*DD(0,I)+QVS(1)*DD(1,I) +QVS(2)*DD(2,I)+QVS(3)*DD(3,I)
	  
	ENDDO
	
!	2nd derivs; use chain rule and utilize that some elements of QSS
!	are identically zero
	IF (NDER > 1) THEN
	
	  DO I=1,NC
	  
	    C0= DD(3,I)*QSS(3,0)
	    C1= DD(2,I)*QSS(2,1)+DD(3,I)*QSS(3,1)
	    C2= DD(1,I)*QSS(1,2)+DD(2,I)*QSS(2,2)+DD(3,I)*QSS(3,2)
	    C3= DD(0,I)*QSS(0,3)+DD(1,I)*QSS(1,3)+DD(2,I)*QSS(2,3)+DD(3,I)*QSS(3,3)
	    DO J=1,I
	    
	      DQ2DN2(J,I) = C0*DD(0,J)+C1*DD(1,J)+C2*DD(2,J)+C3*DD(3,J)
		
	    ENDDO
	    
	  ENDDO
	  
	ENDIF
	endsubroutine  HARDSP
!--------------------------------------------------------------------------------------------------
!
!
	SUBROUTINE QDERS (V, S, Q, QS, QVS, QSS, DQV, DQVV)

!	QDERS calculates the hard-sphere contriubtion to the helmholz
!	residual energy
!
!	Input:
!
!	V:      Volume
!	S:      Vector of sums, SUM_I X_I M_I D(I)^K, K=0,1,2,3
!
!	Output:
!
!	Q:      Helmholz energy (reduced)
!	QS:     S-derivative of Q
!	QVS:    V-derivative of  QS
!	QSS:    Double S-derivative
!	DQV:    V-derivative of Q
!	DQVV:   V-derivative of DQV

	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION S(0:3), QS(0:3),QVS(0:3), QSS(0:3,0:3)

!	Initialize
	DO I=0,3
	  DO J=0,3
	  
	    QSS(J,I) = 0.D0
	    
	  ENDDO
	ENDDO

!	'ZETA-factors'
	VR = 1.D0/V
	S1R = S(1)*VR
	S2R = S(2)*VR
	S3R = S(3)*VR
	
	ZETA = S3R
	Y = 1.D0/(1.D0-ZETA)
	S1RY = S1R*Y
	S2RY = S2R*Y
	VRY =  VR*Y
	Y2= Y*Y
	ZLOG = dLOG(1.D0 - ZETA)
	S22 = S(2)*S(2)
	S23 = S22*S(2)
!
!	Get factor function
	CALL DCALC (V, S(3),F,FV,F3,FVV,F33,F3V)
	S1RY3 = 3.D0*S1RY
	S2RY3 = 3.D0*S2RY

!	Base Q; 3 contributions
	Q = S1RY3*S2R*V + S23*F - S(0)*ZLOG
!
!	Q-derivatives
	DQV = -S1RY3*S2RY + S23*FV - VR*S(0)*(Y-1.D0)
	QS(0) = -ZLOG
	QS(1) = S2RY3
	QS(2) = 3.D0*(S1RY + S22*F)
	PROD3 = S1RY3*S2RY
	QS(3) = PROD3 + S23*F3 + S(0)*VRY

!	Higher derivatives
	QSS(0,3) = VRY
	QVS(0) = VR*(1.D0-Y)
	QSS(1,1) = 0.D0
	QSS(1,2) = 3.D0*VRY
	QSS(1,3) = S2RY3*VRY
	QVS(1) = -QSS(1,3)
	QSS(2,2) = 6.D0*S(2)*F
	S22M3 = 3.D0*S22
	QSS(2,3) = S22M3*F3 + S1RY3*VRY
	QVS(2) =   S22M3*FV - S1RY3*VRY
	PROD33 = 2.D0*PROD3*VRY
	S0VR2 = S(0)*VRY**2
	QSS(3,3) = PROD33 + S23*F33 + S0VR2
	QVS(3) =  -PROD33 + S23*F3V - S0VR2
	DQVV   =   PROD33 + S23*FVV + S(0)*VR**2 *(-1.D0 + Y2)
	DO K=1,3
	  DO I=0,K-1
	  
	    QSS(K,I) = QSS(I,K)
	    
	  ENDDO
	ENDDO
	endsubroutine QDERS
!--------------------------------------------------------------------------------------------------
!
! 
	SUBROUTINE DCALC (V, S3, F, FV, F3, FVV, F33, F3V)

!	Input:
!
!	V:      Volume/R
!	S3:     Sum function S3
!
!	Outout:
!
!	F:      Multiplier function
!	FV:     1st partial derivative of F wrt V
!	F3:     1st partial derivative of F wrt V
!
!	FVV:    2nd partial derivative wrt V
!	F33:    2nd partial derivative wrt S3
!	F3V:    2nd partial cross derivative
!
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	
	VR = 1.D0/V
	ZETA = S3*VR
	X = 1.D0/ZETA
	Y = 1.D0/(1.D0-ZETA)
	VR2 = VR*VR
	X2 = X*X
	Y2 = Y*Y
	VR3 = VR*VR2
	X3 = X2*X
	Y3 = Y2*Y
	VR4 = VR*VR3
	X4 = X3*X
	Y4 = Y3*Y
	XY2 = X*Y2
	ZLOG = dLOG (1.D0-ZETA)
	F = VR2*(XY2 + X2*ZLOG)
	XY = X*Y
	XY3=X*Y3
	XY4=X*Y4
!
!	V-DERIVATIVES
	FV= VR3*(XY2 -2.D0*XY3 + XY)
	X2Y2 = X2*Y2
	FVV = VR4*(-4.D0*XY3+6.D0*XY4+X2-X2Y2)
	X2Y3 = X2*Y3
	X2Y=X2*Y
!
!	S- and cross-derivatives
	F3= VR3*(-X2Y2+2.D0*(XY3-X3*ZLOG)-X2Y)
	F33 = VR4*(2.D0*X3*Y2 + 6.D0*(XY4+X4*ZLOG) +4.D0*(X3*Y-X2Y3) - X2Y2)
	F3V = VR4*(-X2Y2+2.D0*(X2Y3+XY3-X2Y) -6.D0*XY4 +X2Y2)
	
	endsubroutine DCALC
!--------------------------------------------------------------------------------------------------
!
! 
	SUBROUTINE GZLCAL (NC, NDER, NTEMP, DD, DDT, DTT, X, DIA, S, ST, ST2, V,   &
	                   Q, QV, QVV, QT, QTV, QTT, DQDN, DQ2DN2, DQNV, DQNT)
!
!	This subroutine calculates the attractive chain-term and its
!	1st and second derivatives wrt T, V, and N
!
!	Input:
!
!	NC:          No. of components
!	NDER:        Level of comp.derivatives
!	NTEMP:       Level of temp.derivatives
!	DD:          Power of diameter array
!	DDT:         Derivative of DD wrt 1./T
!	DTT:         Derivative of DDT wrt 1./T
!	X:           Molar feed
!	DIA:         Vector of hard-sphere diameters
!	S:           Diameter power sum
!	ST:          1/T-derivative of S
!	ST2:         1/T-derivative of ST
!	V:           Volume
!
!	Output
!
!	Q:           Residual A-contribution
!	QV:          Derivative OF Q wrt V
!	QVV:         Derivative OF QV wrt V
!	QT:          Derivative OF Q wrt 1/T
!	QTV:         Derivative OF QV wrt 1/T
!	QTT:         Derivative OF QT wrt 1/T

!	DQDN:        N-derivatives of Q
!	D2QDN2:      N-derivatives of DQDN
!	DQNV:        N-derivatives of QV
!	DQNT:        N-derivatives of QT
!
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	
	PARAMETER (NCM = 30)
	
	DIMENSION DD(0:3,NCM),DDT(0:3,NCM),DTT(0:3,NCM)
	DIMENSION X(NC),DIA(NCM),S(0:3),ST(0:3),ST2(0:3)
	DIMENSION QS(3),QSS(3,3)
	DIMENSION QST(3)
	DIMENSION QD1_2(NC),QD1_3(NC),FLOG(NC),DQDN(NC),DQNV(NC)
	DIMENSION DQ2DN2(NC,NC),DQNT(NC)
	common /assoc/ NST
	common /PCSAFT_FV/   common_d(NCM), common_DD(0:3,NCM), common_DDT(0:3,NCM), &
	                     common_DTT(0:3,NCM), common_S(0:3), common_ST(0:3),       &
	                     common_ST2(0:3), VR, Y
	
!	Here, S3 = (pi/6)*sum( n(i)*m(i)*d(i)^3 ) = (pi/6)*lambda3
!	that is, S3 = ntot*(zeta3/rho) = zeta3*V
	
	S2 = S(2)           !lambda2*(pi/6)^0.66 = zeta2*V
	S3 = S(3)           !lambda3*pi/6 = zeta3*V
	VMS = V - S3        !
	VMSR = 1.D0/VMS     ! 1/(V - V*zeta3) 
	Y = VMSR*V
	VMSR3 = VMSR**3     ! 1/(V - V*zeta3)^3
	VMSR3F= VMSR3*V     ! V/(V - V*zeta3)^3
	VLOG = dLOG(VMSR3F) ! ln[ V/(V - lambda3)^3 ]
	VR = 1D0/V
	Q = 0.D0
	DO I=1,3
	
	  QS(I) = 0.D0
	  QST(I) = 0.D0
	  DO K=1,3
	  
	    QSS(K,I) = 0.D0
	    
	  ENDDO
	  
	ENDDO
	QT = 0.D0
	QTV = 0.D0
	QTT = 0.D0

!	This is somewhat messy as we need the component loop
      DO I = 1, NC
      
	  XMFAC = DD(0,I)-1.D0 !(m(i) - 1)
	  XFAC = X(I)*XMFAC
	  DIAH = .5D0*DIA(I)
	  D2 = DIAH*S2
	  FN1I = 1.D0/(VMS + D2)  !1/(V - V*zeta3 + d(i)/2*lambda2)
	  FN2I = 1.D0/(VMS +2.D0*D2)
	  FN1I2 = FN1I*FN1I
	  FN2I2 = FN2I*FN2I
	  GLOG = VLOG - dLOG (FN1I*FN2I)
	  FLOG(I) = GLOG*XMFAC    !ln[g(ii)]]
	  Q = Q + X(I)*FLOG(I)
	  ATCI =   DIA(I)*ST(2) + DDT(0,I)*S(2)

!	Contributions to T-derivatives
	  IF (NTEMP > 0) THEN
	  
	    AT1 = -ST(3) + .5D0*ATCI
	    AT2 = -ST(3) +      ATCI
	    AT3 = -ST(3)
	    FR1 = FN1I*AT1
	    FR2 = FN2I*AT2
	    FR3 = VMSR*AT3
	    FACT3 = 3.D0*VMSR*AT3
	    DQNT(I)  = XMFAC*( FR1 +FR2 -3.D0*FR3)
	    QT = QT + DQNT(I)*X(I)

!	    2nd 1/T-derivatives here
	    ATCIT = DIA(I)*ST2(2) + 2.D0*DDT(0,2)*ST(2) + DTT(0,I)*S(2)
	    AT1T = -ST2(3) + .5D0*ATCIT
	    AT2T = -ST2(3) +      ATCIT
	    AT3T = -ST2(3)
	    FR1T = -FR1*FR1 + FN1I*AT1T
	    FR2T = -FR2*FR2 + FN2I*AT2T
	    FR3T = -FR3*FR3 + VMSR*AT3T
	    QTT = QTT + XMFAC*X(I)*(FR1T+FR2T-3.D0*FR3T)
	    CONT2_T = XMFAC * (DDT(0,I)*(.5D0*FN1I + FN2I) + DIA(I)*(-.5D0*FN1I2*AT1 - FN2I2*AT2))
	    QTVA =  XMFAC*(-FN1I2*AT1-FN2I2*AT2+VMSR*FACT3)
	    QTV = QTV + X(I)*QTVA
	    CONT3_T = -QTVA
	    QST(2) = QST(2) + X(I)*CONT2_T
	    QST(3) = QST(3) + X(I)*CONT3_T
	    
	  ENDIF
!
!	  End T-derivatives
	  TV = FN1I + FN2I -3.D0*VMSR
!
!	  1st derivatives
	  CONT1 = XMFAC*(TV+VR)
	  CONT3 =-XMFAC*TV
	  CONT2 = XMFAC*(DIAH*FN1I+DIA(I)*FN2I)
	  QS(1) = QS(1) + X(I)*CONT1
	  QS(2) = QS(2) + X(I)*CONT2
	  QS(3) = QS(3) + X(I)*CONT3
	  DQNV(I) = CONT1
	  QD1_2(I) = CONT2
	  QD1_3(I) = CONT3

!	  2ND COMP. DERIVATIVES
	  DI2 = DIA(I)**2
	  TERM = FN1I2 + FN2I2  - 3.D0*VMSR**2
	  QSS(1,1) = QSS(1,1) -XFAC*(TERM + VR*VR)
	  QSS(2,2) = QSS(2,2) -XFAC*DI2*(0.25D0*FN1I2 + FN2I2)
	  QSS(1,2) = QSS(1,2) -XFAC*DIA(I)*(0.5D0*FN1I2 + FN2I2)
	  QSS(3,3) = QSS(3,3) -XFAC*TERM
	  
	ENDDO
	QV = QS(1)
	QSS(2,3) = -QSS(1,2)
	QSS(1,3) = -QSS(3,3)
	QVV = QSS(1,1)
	DO J=1,3
	  DO K=J+1,3
	  
	    QSS(K,J) = QSS(J,K)
	    
	  ENDDO
	ENDDO
!
!	1st composition derivatives
	DO I=1,NC
	
	  DQDN(I) = FLOG(I) + QS(2)*DD(2,I) + QS(3)*DD(3,I)
	  DQNV(I) = DQNV(I) + QSS(1,2)*DD(2,I) + QSS(1,3)*DD(3,I)
	  
	ENDDO
	
!	T - if activ
	IF (NTEMP > 0) THEN
	  DO I=1,NC
	  
	    DQNT(I) = DQNT(I) + QST(2)*DD(2,I) + QST(3)*DD(3,I) + QS(2)*DDT(2,I) + QS(3)*DDT(3,I)
	    
	  ENDDO
	ENDIF

!	2nd composition derivatives
      IF (NDER > 1) THEN
	  DO J=1,NC
	    DO K=1,J
	    
	      DQ2DN2(K,J)= QD1_2(J)*DD(2,K) + QD1_3(J)*DD(3,K) + QD1_2(K)*DD(2,J) + QD1_3(K)*DD(3,J)
	      CF1 = DD(2,J)*QSS(2,2) + DD(3,J)*QSS(2,3)
	      CF2 = DD(2,J)*QSS(2,3) + DD(3,J)*QSS(3,3)
	      DQ2DN2(K,J) = DQ2DN2(K,J) + CF1*DD(2,K) + CF2*DD(3,K)
		
	    ENDDO
	  ENDDO
	ENDIF
	if (NST > 0) then

! 	  Variables used for RDF in association strenght
	  common_S(0:3) = S(0:3)
	  common_ST(0:3) = ST(0:3)
	  common_ST2(0:3) = ST2(0:3)
	  do i = 1, NC
	    
	    common_d(i) = dia(i)
	    common_dd(0:3,i) = DD(0:3,i)
	    common_ddT(0:3,i) = DDT(0:3,i)
	    common_DTT(0:3,i) = DTT(0:3,i)
	    
	  enddo	  
	
	endif	
	endsubroutine GZLCAL

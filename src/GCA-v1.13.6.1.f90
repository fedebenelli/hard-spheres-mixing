subroutine GCA_version (versGCA, thermo_name)

	character(20)  ::thermo_name, versGCA

	versGCA = "1.13.5.0 (2.0RC13)"
	thermo_name = "GCA-EOS"
	return
endsubroutine GCA_version
!---------------------------------------------------------------------------------------
!
!
!
!
!
!	Group Contribution with Association Equation of State (GCA-EoS)
!	Basada en la GC-EoS de Steen Skjold-Jørgensen con el agregado de una contribución
!	asociativa derivada de la ecuación de Chapman (1989,1990), pero formulada a con-
!	tribución grupal.
!
!
SUBROUTINE GCEOS (NC, NG, NST, iDer, iTemp, T, P, XN, phi, dLPhi, dLPhiT, DLPHIP, CZ, IGz, MTyp, IC)
!-----------------------------------------------------------------------
!	SUBRUTINA GCEOS calcula el logaritmo de los coeficientes de fugaci-
!	dad y sus derivados con respecto a la temperatura, presión y compo-
!	sición una mezcla de componentes NC.
!	XN especifica el número de moles de cada componente.
!	IOPT especifica si el cálculo debe ser hecho para la fase gaseosa
!	(MTYP =- 1) o la fase lí­quida (MTYP = 1).
!
!	iDer = 1: derivadas 1ª respecto del número de moles
!	     = 2:     "     2ª
!
!	iTemp = 1: derivadas respecto de T
!
!	iGz: inicialización del cálculo. Detalles en subr. ZMAX
!
!
!	PHI(I)    = ln del coeficiente de fugacidad del comp i
!     DLPHIT(I) = derivada del ln(phi i) respecto de T
!     DLPHIP(I) = derivada del ln(phi i) respecto de P
!     DLPHI(I,K)= derivada del ln(phi i) respecto de nk
!-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z)
!
	parameter (NCM = 30, NGM = 30, NGAM = 14, NSM = 24)
	integer sigma

!
	dimension x(NCM),xms(NCM),Ps(NCM,NGM),xN(NCM),phi(NCM),            &
     &          help3(NCM,NGM),help7(NCM,NGM),dLam1(NCM),dLam2(NCM),     &
     &          dLam3(NCM),Z(2),dYdn(NCM),dYvdn(NCM)
!
	DIMENSION DLPHI(NCM,NCM),HELP8(NCM,NGM),HELP9(NCM,NGM),            &
     &          HELP10(NCM,NGM),DFDNN(NCM,NCM),DPDN(NCM),DFDN(NCM)
!
	DIMENSION DH3DT(NCM,NGM),DH7DT(NCM,NGM),DH2DT(NGM),DH4DT(NGM),     &
     &          DH5DT(NGM),DH6DT(NGM),DLAMT1(NCM),DLAMT2(NCM),           &
     &          DLAMT3(NCM),DYDTN(NCM)
!
	real*8 LAMBDA2,LAMBDA,invXs_m(NSM,NSM)
	dimension aux(NSM),dXs_dn(NSM,NCM),dXs_ds(NSM,NSM),                &
     &          d2Fassoc_dsdT(NSM)


	DIMENSION DLPHIT(NCM),DLPHIP(NCM),DFDNT(NCM)
!
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),rNyT(NGM),QT,    &
     &           TETA(NGM),HELP6(2,NGM),HELP11(2,NGM),HELP12(2,NGM),     &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3
!
	COMMON/COORD/ZZ
!
!	Variables especí­ficas del término dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)
!
!	Versión asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Versión asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!	Variables ya calculadas en la subrutina ZMAX
	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM), omega(NCM)
	
	common/ScndDer/dPV, dPDT, dPDn
	

	PI=3.1415926536D0
	RT=T*R
	dPdV = 0
	dPdT = 0
	dPdn = 0
	dXs_ds(:NST,:NST) = 0.D0
!
!	Calcula el número total de moles, fracciones molares, moles
!	de superficie total, y las fracciones superficie de area
	xNtot = sum(xN(:NC))
!
!	Define la fase y su factor de compresibilidad
	IPHAS=1+(MTYP+1)/2
	if (IGZ /= 0 .AND. MTYP /= 0)Z(IPHAS) = CZ
	CALL ZMAX (Z, MTYP, T, P, XN, IGZ, NC, NG ,NST, IC)
	IPHAS = 1 + (MTYP*IC + 1)/2
	if (MTYP == 0) IPHAS=1+(IC+1)/2
	CZ=Z(IPHAS)
!
!	Cálculo de coeficientes de fugacidad.
	rNyT(:NG) = rNyT(:NG)*xNtot

	QT = QT*XNTOT               	!Área total de los grupos atractivos
	V = Z(IPHAS)*RT/P*XNTOT     	!Volumen de la fase
	DPV = -DPDV(IPHAS)*XNTOT/V/V	!(dP/dV)_(T,n): viene del COMMON "ZZZ"
!
	DO 22 I=1,NC

	  XSUM=0.D0
	  DO 2 K=1,NG

	    PS(I,K)=0.D0             	!PS(i,k) = ny(k,i)*q(k)
!	    if (NYN == 0) GOTO 2
	    if (NY(I,K) /= 0) then

!	      PSIK=DFLOAT(NYN)*Q(K)		!Necesita convertir Ny en número decimal, sino FORTRAN es capaz de convertirle Ps(i,k) en un entero...
	      PS(I,K)=DFLOAT(NY(I,K))*Q(K)	!(PS)ik	!Necesita convertir Ny en número decimal, sino FORTRAN es capaz de convertirle Ps(i,k) en un entero...
	      XSUM=XSUM+PS(I,K)

	    endif

2	   enddo
	   XMS(I)=XSUM                   	!MS(i) = SUM(ny(k,i)*q(k), k=1,NG)

22	enddo
!
!	Variable auxiliar para tau_ij:
      DNOM=QT/V/RT
!
	dPdTat = 0.D0
	if (NST  >  0) then
!
!	  Esta subrutina trabaja con volumen total en vez de molar como
!	  la sub. ZMAX
	  sm(:NST) = xnTot*sm(:NST)
	  dXs_dV(iPhas,:NST) = dXs_dV(iPhas,:NST)/xnTot
	  sm_Xs(iPhas,:NST) = sm_Xs(iPhas,:NST)*xnTot

	endif
!
!     Variables auxiliares para el calculo del ln(coef.fug)i y sus derivadas
!     con respecto a T,P y nro. de moles del comp. j
!
!     Contribución atractiva
	DO 6 I=1,NC

	  XMSI=XMS(I)
	  DO 61 K=1,NG

	    HEL3K=0.D0
	    HEL7K=0.D0
	    HEL9=0.D0
	    HEL10=0.D0
	    DO 23 J=1,NG

	      if (ny(i,j) > 0) then

	        TAU=E(IPHAS,J,K)
	        AJK=A(J,K)*DNOM
	        PSIJ=PS(I,J)*TAU
	        HEL7K=HEL7K+PSIJ         !H7ik
	        PTIJ=PSIJ*DLOG(TAU)
	        PSIJ=PSIJ*AJK
	        HEL3K=HEL3K+PSIJ         !H3ik
	        HEL10=HEL10+PTIJ
	        HEL9=HEL9+PTIJ*AJK

	      endif
!
23	    enddo
	    H4=HELP4(IPHAS,K)
	    HEL3K=HEL3K/H4
	    HEL7K=HEL7K/H4
	    HEL9=HEL9/H4
	    HEL10=HEL10/H4
	    HELP3(I,K)=HEL3K        !H3ik/H4k
	    HELP7(I,K)=HEL7K        !H7ik/H4k
	    HELP9(I,K)=HEL9
	    HELP10(I,K)=HEL10
	    HELP8(I,K)=HEL7K-XMSI*(1.D0-HELP6(IPHAS,K))

61	  enddo

6	enddo
!
	if (iTemp /= 0) then
!
!  Derivadas del término atractivo con respecto de T
	  DH3DT(:NC,:NG)=0.D0
	  DH7DT(:NC,:NG)=0.D0
	  DHEL=0.D0
	  DO 300 J=1,NG

	    H30=0.D0
	    H31=0.D0
	    H32=0.D0
	    H33=0.D0
	    H34=0.D0
	    H35=0.D0
	    H36=0.D0
	    AA=A(J,J)
	    DAT=DADT(J,J)
	    DO 301 K=1,NG

	      TAU=E(IPHAS,K,J)
	      TET=TETA(K)
	      DAKJ=DADT(K,J)
	      ARG=DLOG(TAU)
	      AKJ=A(K,J)
	      DDEL=DAKJ-DAT
	      DARG=0.D0
	      DCRT=DABS(AKJ-AA)
	      if (DCRT  >  1.D-2) DARG=ARG*(DDEL/(AKJ-AA)-1.D0/T)
	      ALF=ALFA(K,J)
	      TETAU=TAU*TET
	      HELPI=TAU*DARG
	      HELP=TETAU*DARG
	      H30=H30+HELP
	      HELPI1=HELPI*AKJ*DNOM
	      HELP1=HELP*AKJ*DNOM
	      H31=H31+HELP1
	      H32=H32+HELP*ARG
	      H33=H33+HELP1*ARG
	      HEL2=TETAU*DAKJ*DNOM
	      HELPI2=TAU*DAKJ*DNOM
	      H34=H34+HEL2
	      H35=H35+HEL2*ARG
	      H36=H36+TETAU*DDEL*ALF*DNOM
	      DO 3010 I=1,NC

	        if (ny(i,k) /= 0) then

	          PSIK=PS(I,K)
	          DH3DT(I,J)=DH3DT(I,J)+HELPI1*PSIK+HELPI2*PSIK  !d(H3ij/H4j)/dT
	          DH7DT(I,J)=DH7DT(I,J)+HELPI*PSIK               !d(H7ij/H4j)/dT

	        endif

3010	      enddo

301	    enddo
	    H2=HELP2(IPHAS,J)
	    H4=HELP4(IPHAS,J)
	    H5=HELP5(IPHAS,J)
	    DH2T=(H31+H34)/H4
	    DH4T=H30/H4
	    DH5T=(H33+H35+H31)/H4
	    H6=HELP6(IPHAS,J)
	    DH6T=-H6/T+(H32+H36)/H4
	    DH2DT(J)=DH2T             !d(H2j/H4j)/dT
	    DH4DT(J)=DH4T             !d(H4j)/dT
	    DH5DT(J)=DH5T             !d(H5j/H4j)/dT
	    DH6DT(J)=DH6T             !d(H6j/H4j)/dT
!
!---------(dP/dT)att = -RT*(d2(F)att/dVdT) - R*(d(F)att/dV)
!
	    HEL=DH5T+DH2T-DH2T*H6-H2*DH6T-(H5+H2-2.D0*H2*H6 )*DH4T
	    DHEL = DHEL + HEL*rNYT(J)*Q(J)    !- 2*V/(ZZ*RT)*(dP/dT)att

300	  enddo
!
	  DPDT=XNTOT*R/V-ZZ/2.D0*RT/V*DHEL  !n*R/V + (dP/dT)att (falta contribucion asociativa y repulsiva)
	  dPatdT = -ZZ/2.D0*RT/V*DHEL
!
!	  Acá termina el if (iTemp == 0)goto 306
!--------------------------------------------------------------------------

	endif
	dXs_dn(:NST,:NC) = 0D0
	dPdTas = 0D0
306	if (NST > 0) then

	  if (iDer == 2 .OR. iTemp == 1) then
!
!  Calculation of the number of associating site moles derivatives of the non-bonded fraction. This is required for the
!  calculation of d2[ln(phi_k)]/[dni dnj], and for d[ln(phi_k)]/dT and dP
!
	    if (NST == 1) then
	    
	      dXs_ds(1,1) = -dXs_dV(iPhas,1)*V/sm(1)
	
	    elseif (NST == 2 .AND. Delta(1,1) <= 0 .AND. Delta(2,2) <= 0) then
	    
	      dXs_ds(major,major) = -Xs(iPhas,major)**2 * Delta(1,2)/V/2D0*( (2D0 - b_aux(iPhas))/root(iPhas) - 1D0)
	      dXs_ds(major,3-major) = -Xs(iPhas,major)**2 *Delta(1,2)/V/2D0*(1D0 + b_aux(iPhas)/root(iPhas))
	      if (sm(major) == sm(3-major)) then
	  
	        dXs_ds(3-major,3-major) = dXs_ds(major,major)
	        dXs_ds(3-major,major) = dXs_ds(major,3-major)
! 	        dXs_ds(major,3-major) = dXs_ds(major,major)
	  
	      else
	  
	        dXs_ds(3-major,major) = -Xs(iPhas,major)**2 *Delta(1,2)/V*(Xs(iPhas,major) + sm(major)*dXs_ds(major,major))
	        dXs_ds(3-major,3-major) = -Xs(iPhas,major)**2 *Delta(1,2)/V*sm(2)*dXs_ds(major,3-major)
		  
	      endif
	
	    else
	
	      do i = 1, NST
	
	        if (sm(i) >= 1D-16) dXs_ds(i,:NST) = -sm(i)*Xs(iPhas,:NST)*Delta(i,:NST)/V/xnTot
! 	        dXs_ds(i,:NST) = -sm(i)*Xs(iPhas,:NST)*Delta(i,:NST)/V/xnTot
	        dXs_ds(i,i) = dXs_ds(i,i) + (1D0/Xs(iPhas,i) - 1D0 - dot_product( sm_Xs(iPhas,:NST), Delta(i,:NST) )/V)/xnTot
	      
	      enddo
!
!	      Conviene hallar uno u otro primero, dependiendo qué sea más
!	      grande: NC o NST
! 	      dXs_ds(:NST,:NST) = dXs_ds(:NST,:NST)/xnTot
! 	      if (NST < NC) then

	      do i = 1,NST

	        dXs_ds(:NST,i) = -dXs_ds(:NST,i)
	        call LUBksb (H(iPhas,:NST,:NST), NST, NST, indx(iPhas,:NST), dXs_ds(:NST,i))

	      enddo
!
!	      dXs/dn = dXs/ds*ds/dn, ds/dn = sigma

! 	    else
! !
! !	      matriz dg(i)/dn(l):
! 	      dXs_dn(:NST,:NC) = matmul(-dXs_ds(:NST,:NST),dfloat(sigma(:NST,:NC))) !/XNTOT
! 	      do i = 1,NC
! 
! 	        call LUBksb (H(iPhas,:NST,:NST), NST, NST, indx(iPhas,:NST), dXs_dn(:NST,i))
! 
! 	      enddo
! 
	    endif
	    
	    dXs_dn(:NST,:NC) = matmul(dXs_ds(:NST,:NST),dfloat(sigma(:NST,:NC)))

	    if (iTemp == 1) then
!
!  Calculation of the temperature derivative of the Helmholtz free energy:
!  (dP/dT)Vn = -(d2A/dTdV)n
!  (d2A/dTdV)n = RT (d2(A/RT)/dTdV)n + R (d(A/RT)/dV)Tn
!  dFVasoc = d(Aasoc/(nRT))/dv = d(Aasoc/RT)/dV
	      dPdTas = - R*dFVas(iPhas)
	      do i = 1, NST  
! 
	        dPdTas = dPdTas + RT*sm_Xs(iPhas,i)*SUM( sm_Xs(iPhas,:NST)*dDeldT(i,:NST)*(dXs_dV(iPhas,:NST)/Xs(iPhas,:NST) - 5D-1/V) )/V
	     
	      enddo
	      do k = 1, NST

	        aux(:NST) = sm(:NST)*dXs_ds(:NST,k)
	        aux(:NST) = matmul( dDeldT(:NST,:NST), aux(:NST) )
	        d2Fassoc_dsdT(k) = dot_product( sm_Xs(iPhas,:NST), aux(:NST) )
	        d2Fassoc_dsdT(k) = -(d2Fassoc_dsdT(k) + Xs(iPhas,k)*dot_product( sm_Xs(iPhas,:NST), dDeldT(:NST,k) ))/V
	        
	        
! 	        d2Fassoc_dsdT(k) = 0.D0
! 	        do i = 1, NST
! 	          if (sm(i) >= 1E-16) then
! ! 	          if (sm(i) >= 1E-16) then
! 	    
! 	            d2Fassoc_dsdT(k) = d2Fassoc_dsdT(k) - sm_Xs(iPhas,i)*Xs(iPhas,k)*dDeldT(k,i)
! 	            do j = 1, NST
! 	              if (sm(j) >= 1E-16) then
! 	    
! 	                d2Fassoc_dsdT(k) = d2Fassoc_dsdT(k) - sm_Xs(iPhas,i)*sm(j)*dXs_ds(j,k)*dDeldT(i,j)
! 	                
! 	              endif
! 	            enddo
! 	            
! 	          endif
! 	        enddo	  
! 	        d2Fassoc_dsdT(k) = d2Fassoc_dsdT(k)/V

	      enddo	      
!
	    endif

	  endif

	endif
!
!  Contribución repulsiva
	XLAM1=XLAM1*XNTOT    !lambda1
	XLAM2=XLAM2*XNTOT    !lamdda2
	XLAM3=XLAM3*XNTOT    !lamdda3
	TLAM1=0.D0
	TLAM2=0.D0
	TLAM3=0.D0
	DO 9 I=1,NC

	  DIA=D(I)
!	  if (ITEMP == 0) GOTO 311
	  if (iTemp /= 0) then
	    
!  d(dc)/dT y var. aux. para d(lambda_k)/dT
	    DIAT=DT(I)
	    DELDT=DIAT
	    XDT=XN(I)*DIAT
	    DLAMT1(I)=DELDT
	    TLAM1=TLAM1+XDT
	    DELDT=2.D0*DELDT*DIA
	    DLAMT2(I)=DELDT
	    XDT=XDT*DIA
	    TLAM2=TLAM2+XDT
	    DLAMT3(I)=1.5D0*DELDT*DIA
	    TLAM3=TLAM3+XDT*DIA
!
	  endif
311	  continue
	  DLAM1(I)=DIA/XLAM1         !d(lambda1)/dni / lambda1
	  DIAV=DIA*DIA
	  DLAM2(I)=DIAV/XLAM2        !d(lambda2)/dni / lambda2
	  DLAM3(I)=DIAV*DIA

9	enddo
	PI6 = PI/6.D0/V
	Y = 1.D0/(1.D0-PI6*XLAM3)
	Y2 = Y*Y
	DYDV = -Y2*PI6*XLAM3/V
!	if (ITEMP == 0) GOTO 312
	if (iTemp /= 0) then
!
	  TLAM2=TLAM2*2.D0
	  TLAM3=TLAM3*3.0D0
	  DYDT=Y2*PI6*TLAM3        !dY/dT
	  TLAM3=TLAM3/XLAM3        !d(lambda1)/dt / lambda1
	  TLAM2=TLAM2/XLAM2        !d(lambda2)/dt / lambda2
	  TLAM1=TLAM1/XLAM1        !d(lamdda3)/dt / lambda3
	  DYDVT=DYDV*TLAM3+2.D0/Y*DYDV*DYDT !d2Y/dTdV
!
	endif
312	continue
	DO 98 I = 1, NC

	  DYDN(I)=Y2*PI6*DLAM3(I)     	!dY/dni
	  if (iTemp /= 0) then
!
	    DYDTN(I) = 2.D0/Y*DYDN(I)*DYDT + DLAMT3(I)/DLAM3(I)*DYDN(I) !d2Y/dTdni
!
	  endif
298	  DYVDN(I) = (2.D0/Y*DYDV - 1.D0/V)*DYDN(I)	!d2Y/dVdni
	  DLAM3(I) = DLAM3(I)/XLAM3       		!d(lambda3)/dni / lambda3

98	enddo
	R0=XLAM2/XLAM3
	R1=XLAM1*R0
	R3=R0*R0*XLAM2
	R4=Y2-Y-DLOG(Y)
	R6=2.D0*Y-1.D0-1.D0/Y
	R33=2.D0+1.D0/Y2
	R34=Y-1.D0
	if (iTemp /= 0) then

	  R30=TLAM1+TLAM2-TLAM3
	  R31=3.D0*TLAM2-2.D0*TLAM3
!
!	  (dP/dT)rep = -R*d(F)rep/dV - R*T*d2(F)rep/dVdT
	  DPDTR=PREP(IPHAS)/T+DYDVT/DYDV*PREP(IPHAS)-RT*DYDV*(3.D0*R1*R30+R3*R31*R6+R3*R33*DYDT-XNTOT/Y2*DYDT)   !(dP/dT)rep
	  DPDT = DPDT + DPDTR + dPdTas  !n*R/V+(dP/dT)att+(dP/dT)rep+(dP/dT)asoc

	endif
314	continue
!-------------------------------------------------------------------------------
!
!
!  Cálculo del ln(coef.fug)i y sus derivadas con respecto a T,P y nro. de moles del comp. j
	DO 8 I=1,NC

	  dFdnas = 0D0
! 	  if (maxval(sigma(:NST,i)) > 0) dFdnas=1.D0	!lo hallo por productoria, no sumatoria
	  dPdnas=0.D0
	  dFdnTas=0.D0
	  DFDTA=0.D0
	  DFDNI=0.D0
	  DPDNI=0.D0
	  DYDNI=DYDN(I)
	  XMSI=XMS(I)
	  NC1=I
	  DO 54 K = NC1,NC

	    DFDNN(I,K) = 0.D0

54	  enddo
	  if (iTemp /= 0) then
!
	    TLAM1I=DLAMT1(I)/XLAM1
	    TLAM2I=DLAMT2(I)/XLAM2
	    TLAM3I=DLAMT3(I)/XLAM3
	    DYDTNI=DYDTN(I)
!
	  endif
299	  CONTINUE
!
!  Contribucion atractiva
	  DO 50 J = 1,NG

	    TET=TETA(J)
	    PSIJ=PS(I,J)
	    H3IJ=HELP3(I,J)
	    H7IJ=HELP7(I,J)
	    H2=HELP2(IPHAS,J)
	    H4=HELP4(IPHAS,J)
	    H5=HELP5(IPHAS,J)
	    H6=HELP6(IPHAS,J)
	    H8IJ=HELP8(I,J)
	    if (iDer > 1 .OR. iTemp /= 0) then
!
	      H12=HELP12(IPHAS,J)
	      H10IJ=HELP10(I,J)
	      H9IJ=HELP9(I,J)
	      H11=HELP11(IPHAS,J)
	      H21=RNYT(J)*Q(J)
	      H20=H5+H2-H2*H6

	    endif
51	    HDER=-H2*(H7IJ-XMSI*(1.D0-H6))+H3IJ+XMSI*H5
	    DELI=-(HDER*TET+PSIJ*H2)
	    DFDNI=DFDNI+DELI                 !2./zz*d(F)att/dni
	    if (iTemp > 0) then
!
	      DH2T=DH2DT(J)
	      DH4T=DH4DT(J)
	      DH5T=DH5DT(J)
	      DH6T=DH6DT(J)
	      HELPT=PSIJ*(DH2T-H2*DH4T)+TET*(DH3DT(I,J)/H4+XMSI*DH5T-(H3IJ+XMSI*H5)*DH4T)
	      HELPT=HELPT-TET*H8IJ*(DH2T-2.D0*H2*DH4T)-TET*H2*(DH7DT(I,J)/H4-DH4T*XMSI+DH6T*XMSI)
	      DFDTA=DFDTA-DELI/T-HELPT         !2./zz*d2(F)att/dnidT
!
	    endif
316	    CONTINUE
	    if (iDer > 1 .OR. iTemp > 0) then
!
	      DP = -PSIJ*H20
	      DP=DP-H21/QT*(XMSI*(3.D0*H5+H2+H11)+H9IJ+H3IJ)
	      DP=DP+H21/QT*(H6*(H3IJ+2.D0*XMSI*H5)+H2*(H10IJ+XMSI*H12)+H7IJ*(H5+H2)+3.D0*XMSI*H2*H6)
	      DP=DP-H21/QT*2.D0*(H7IJ*H2*H6+XMSI*H2*H6*H6)
	      DPDNI=DPDNI+ZZ/2.D0*DP*RT/V      !-R*T*d2(F)att/dVdni
!
	      if (iDer > 1) then

	        H3H5=(H3IJ+XMSI*H5)/QT
	        H2H8=H2*H8IJ/QT
	        H7H6=TET*H2*(H7IJ+H6*XMSI)
	        H8=TET*H8IJ
	        H2H4=TET*H2
	        DO 53 K=NC1,NC

	          XMSK=XMS(K)
	          PSKJ=PS(K,J)
	          H9KJ=HELP9(K,J)
	          H3KJ=HELP3(K,J)
	          H7KJ=HELP7(K,J)
	          H10KJ=HELP10(K,J)
	          D2=(H3KJ+XMSK*H5)/QT
	          D3=XMSK/QT*(H9IJ+H3IJ)
	          D4=(H7KJ-XMSK*(1.D0-H6))/QT
	          D5=(H9KJ+XMSK*(H11+H5))/QT
	          D6=(H10KJ+XMSK*H12)/QT
	          D7=XMSK*H10IJ/QT
	          DELIK=-(PSIJ*D2+(PSKJ-TET*XMSK)*H3H5)
	          DELIK=DELIK+(PSKJ-TET*XMSK)*H2H8-D4*H7H6+D2*H8-DELI*D4
	          DELIK=DELIK+H2H4*(D7+XMSI*D6)-TET*(D3+XMSI*D5)
	          DFDNN(I,K)=DFDNN(I,K)+DELIK      !2./zz*d2(F)att/dnidnj

53	        enddo

	      endif

	    endif

50	  enddo
	  if (iDer == 2) then

	    do k = NC1, NC

!  Antes lo multiplicaba al sumar contribuciones. Desde la aparición de asociación, multiplica 
!  ahora y luego suma... cosas heredadas...
	      dFdnn(i,k) = dFdnn(i,k)*zz/2.D0      !d2(F)att/dnidnj.

	    end do

	  end if

!  Contribución asociativa
	  if (maxval(sigma(:NST,i)) > 0) then

!  Coeficiente de fugacidad:
	    do j = 1, NST
!
!  dF/dn(i)  = SUM(sigma(j,i)*ln(Xs(j), j=1..NST) = ln [PROD( Xs(j)^sigma(j,i), j=1..NST)]
	      dFdnas = dFdnas + dlog(Xs(iPhas,j))*dfloat(sigma(j,i))
! 	      dFdnas = dFdnas*Xs(iPhas,j)**dfloat(sigma(j,i))

	    enddo
! 	    dFdnas = dlog(dFdnas)
	    if (iTemp == 1 .OR. iDer == 2) then
!
!  dPasoc/dn(i) = -RT d2F/(dn(i)dV)
	      dPdnas = -RT*sum(dfloat(sigma(:NST,i))*dXs_dV(iPhas,:NST)/Xs(iPhas,:NST))


	    endif
!
!  d2(F)ass/dnidT
	    if (iTemp == 1) then

! 	      dFdnTas = sum(dfloat(sigma(:NST,i))*dXs_dT(:NST)/Xs(iPhas,:NST))
	      dFdnTas = dot_product( dfloat(sigma(:NST,i)), d2Fassoc_dsdT(:NST) )
! 	      dFdnTas = 0.D0
! 	      do k = 1, NST
! 	
! 	        d2Fassoc_dsdT_k = 0.D0
! 	        if (sigma(k,i) > 0) then
! 	          do j = 1, NST
! 	            if (sm(j) >= 1E-16 .AND. Delta(k,j) >= 1E-16) then
! 	     
! 	              d2Fassoc_dsdT_k = d2Fassoc_dsdT_k + sm_Xs(iPhas,j)*Xs(iPhas,k)*dDeldT(k,j)
! 	              do l = 1, NST
! 	                if (sm(l) >= 1E-16 .AND. Delta(j,l) >= 1E-16) then
! 	      
! 	                  d2Fassoc_dsdT_k = d2Fassoc_dsdT_k + sm_Xs(iPhas,j)*sm(l)*dXs_ds(l,k)*dDeldT(j,l)
! 	                  
! 	                endif
! 	              enddo
! 	              
! 	            endif
! 	          enddo	  
! 	          dFdnTas = dFdnTas - d2Fassoc_dsdT_k*dfloat(sigma(k,i))/V
! 	        endif
! 	      enddo
	      
	    endif
	    if (iDer == 2) then
	      do k = NC1, NC

	        aux_d2fn = sum(dfloat(sigma(:NST,i))*dXs_dn(:NST,k)/Xs(iPhas,:NST))
!
!  d2(F)att/dnidnj+d2(F)asoc/dnidnj
	        dFdnn(i,k) = dFdnn(i,k) + aux_d2fn
!
	      enddo
	    endif

	  endif
!
!  Contribucion repulsiva
	  DL1=DLAM1(I)
	  DL2=DLAM2(I)
	  DL3=DLAM3(I)
	  R2=DL1+DL2-DL3
	  R20=R2+DYDNI/R34
	  R20=R20*R1*3.D0
	  R5=3.D0*DL2-2.D0*DL3
	  DFDNR= R34*R20+R3*R4*R5+R3*R6*DYDNI+DLOG(Y)+XNTOT/Y*DYDNI                                      !d(F)rep/dni
!	  if (ITEMP == 0) GOTO 317
	  if (iTemp > 0) then

	    DL1TN=TLAM1I-TLAM1*DL1
	    DL2TN=TLAM2I-TLAM2*DL2
	    DL3TN=TLAM3I-TLAM3*DL3
	    DFDTR=DYDT*R20+R34*R30*R20+3.D0*R34*R1*(DL1TN+DL2TN-DL3TN+(DYDTNI-DYDT*DYDNI/R34)/R34)
	    DHEL=R31*R4*R5+DYDT*R6*R5+R4*(3.D0*DL2TN-2.D0*DL3TN)+R6*R31*DYDNI + R33*DYDT*DYDNI+R6*DYDTNI
	    DFDTR=DFDTR+R3*DHEL+DYDT/Y-XNTOT/Y*(DYDNI*DYDT/Y-DYDTNI) !d2(F)rep/dnidT
!
!  d2(A^R/RT)/dn(i)dT:
	    DFDNT(I)=DFDTA*ZZ/2.D0 + DFDTR + dFdnTas

	  endif
317	  CONTINUE
!	  IF (IDER == 1 .AND. ITEMP == 0) GOTO 8
	  if (iDer > 1 .OR. iTemp > 0) then
!
!  -R*T*d2(F)rep/dnidV
	    DPREP = PREP(IPHAS)*DYVDN(I)/DYDV - DYDV*RT*(3.D0*R1*R2 + R3*R6*R5 + R3*R33*DYDNI + 1.D0/Y - XNTOT/Y2*DYDNI)
!  dP/dni
	    DPDN(I) = DPDNI + DPREP + RT/V + dPdnas
!	    if (IDER == 1) GOTO 8
	    if (iDer > 1) then
!
	      DO 56 K = NC1, NC

	        DYDNK=DYDN(K)
	        DLAM2K=DLAM2(K)
	        DLAM1K=DLAM1(K)
	        DLAM3K=DLAM3(K)
	        R21=DLAM1K+DLAM2K-DLAM3K
	        R25=3.D0*DLAM2K-2.D0*DLAM3K
	        DREP=R20*(DYDNK+R34*R21)
	        DL1IK=DLAM1K*DL1
	        DL2IK=DLAM2K*DL2
	        DL3IK=DLAM3K*DL3
	        DYIK=DYDNI*DYDNK
	        DREP=DREP + 3.D0*R34*R1*( -DL1IK - DL2IK + DL3IK) - 3.D0*R1*(DYIK/R34 - DYIK*2.D0/Y)
	        DREP=DREP+R4*R5*R25*R3
	        DREP=DREP+R3*R5*R6*DYDNK+R3*R4*(2.D0*DLAM3K*DL3 - 3.D0*DLAM2K*DL2)
	        DREP=DREP+R3*R25*R6*DYDNI+R3*(R33+2.D0*R6/Y)*DYIK
	        DREP=DREP+(DYDNI+DYDNK)/Y+XNTOT*DYIK/Y2          !d2(F)rep/dnidnj
!
!  d2(A^r/RT)/dn(i)dn(k):
	        DFDNN(I,K)=DFDNN(I,K)+DREP   !(dFdnn_at + dFdnn_as) + dFdnn_rep

56 	      enddo
!
	    endif

	  endif
!  d(A^r/RT)/dn(i):
	  DFDN(I)=DFDNI*ZZ/2.D0+DFDNR +dFdnas

8	enddo
	DO 60 I=1,NC

	  PHI(I)=-DLOG(Z(IPHAS))+DFDN(I)          !ln(coef.de fugac)i
!	  if (ITEMP == 0) GOTO 62
	  if (iTemp > 0) then

	    DLPHIT(I)=DFDNT(I)+DPDN(I)*DPDT/DPV/RT+1.D0/T  !d(ln(coef.fug.)i)/dT
	    DLPHIP(I)=-DPDN(I)/DPV/RT - 1.D0/P               !d(ln(coef.fug.)i)/dP

	  endif
!62	  if (IDER == 1) GOTO 60
62	  continue
	  if (iDer > 1) then
	    DO 6100 K=I,NC

	      DLPHI(I,K)=DFDNN(I,K)+DPDN(I)*DPDN(K)/DPV/RT
	      DLPHI(I,K)=1.D0+XNTOT*DLPHI(I,K)         !n*d(ln(coef.fug.)i)/dnj
	      DLPHI(K,I)=DLPHI(I,K)                    !n*d(ln(coef.fug.)j)/dni

6100	    enddo
	  endif
	!por maxwell son iguales
!
!	      ln(coef fug i) = d(ln coef fug mezcla)/dni
!	      ln(coef fug j) = d(ln coef fug mezcla)/dnj
!	      Entonces la derivada de ln(coef fug)i respecto de j o viceversa [derivada de
!	      ln(coef fug)j respecto de i] termina dando la derivada cruzada de la misma funcion
!           LN(COEF FUG MEZCLA) y por Maxwell deben ser iguales.(Ec dif. exacta)
!
60	enddo
99	RETURN

ENDSUBROUTINE
!-----------------------------------------------------------------------
!
!	NTLPY Calcula la entalpí­a residual a P y T de una mezcla de NC com-
!	ponentes.
!	Está definida como
!
!	      r     conf          conf,GI
!	     H   = H   (T,P,n) - H       (T,P,n)
!
!	donde Hconf es la entalpí­a configuracional.
!	xN es el vector de número de moles de la mezcla (como siempre).
!	mTyp especifica la fases: -1 gas
!	                           1 lí­quido
!	                           0 fase con menos energí­a libre de Gibbs
!
!	Las unidades, como en toda la GCA, están en [atm cm^3], ya que
!	R = 82.05 atm cm^3/(mol K).
!
SUBROUTINE NTLPY (NC, NG, NST, T, P, XN, HRES, CZ, IGZ, MTYP, IC)
!
	IMPLICIT REAL*8(A-H,O-Z)
	parameter (NCM = 30, NGM = 30, NSM = 24)
	integer sigma
!
	DIMENSION XN(NCM),Z(2)
!
	real*8 LAMBDA2,LAMBDA,Id(NST,NST)
	dimension assoc_aux(NSM)
!
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),rNyT(NGM),QT,    &
     &           TETA(NGM),HELP6(2,NGM),HELP11(2,NGM),HELP12(2,NGM),     &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3
!
	COMMON/COORD/ZZ
!
!	Variables especí­ficas del término dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
!	Versión asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Versión asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!	Variables ya calculadas en la subrutina ZMAX
	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
!
!
	PI=3.1415926536D0
	RT=T*R
!
!	Nro. total de moles
	xnTot = sum(xN(:NC))
!
!	Cálculo del factor de compresibilidad
	IPHAS=1+(MTYP+1)/2
	if (IGZ /= 0 .AND. MTYP /= 0) Z(IPHAS)=CZ
	CALL ZMAX(Z,MTYP,T,P,XN,IGZ,NC,NG,NST,IC)
!
!	Variable auxiliar para identificar tipo de fase:
	IPHAS=1+(MTYP*IC+1)/2
	if (MTYP == 0) then

	  IPHAS=1+(IC+1)/2
	  IPHBIS=IPHAS

	END IF
!	Factor compresibilidad de la fase elegida
	CZ=Z(IPHAS)
!
!	Esta subrutina trabaja con volumen total en vez de molar como
!	la sub. ZMAX
!	Nro. de moles del grupo att. k
	rNyT(:NG) = xNtot*rNyT(:NG)
!	Área grupal total
	QT=QT*XNTOT
	if (NST > 0) then

	  sm(:NST) = xnTot*sm(:NST)
	  sm_Xs(iPhas,:NST) = sm_Xs(iPhas,:NST)*xnTot	  

	endif
!	Volumen de la fase
	V=Z(IPHAS)*RT/P*XNTOT
!-----------------------------------------------------------------------
!	Cálculo de la ent.residual:
!
!	                r       2  dF
!	               H  = -R T  ----  + P V - n R T
!	                           dT
!
!	Contribución atractiva
	DNOM=QT/V/RT
	DADTA=0.D0
	ARES=0.D0
	DO 300 J = 1, NG

	  H30=0.D0
	  H31=0.D0
	  H34=0.D0
	  AA=A(J,J)
	  DAT=DADT(J,J)
	  DO 301 K=1,NG

	    TAU=E(IPHAS,K,J)
	    TET=TETA(K)
	    DAKJ=DADT(K,J)
	    ARG=DLOG(TAU)
	    AKJ=A(K,J)
	    DARG=0.D0
	    DDEL=DAKJ-DAT
	    DCRT=DABS(AKJ-AA)
	    if (DCRT > 1.D-2) DARG=ARG*(DDEL/(AKJ-AA)-1.D0/T)
	    TETAU=TAU*TET
	    HELP=TETAU*DARG
	    H30=H30+HELP
	    HELP1=HELP*AKJ*DNOM
	    H31=H31+HELP1
	    HEL2=TETAU*DAKJ*DNOM
	    H34=H34+HEL2

301	  enddo
	  DH2T=H31+H34
	  DH4T=H30
	  ARES=ARES+RNYT(J)*Q(J)*HELP2(IPHAS,J)               !Fatt=(Aatt)/R/T
	  DADTA=DADTA+RNYT(J)*Q(J)*(DH2T-DH4T*HELP2(IPHAS,J))/HELP4(IPHAS,J)

300	enddo
!	d(Fatt)/dT
	dAdTa=-ZZ/2.D0*(dAdTa-Ares/T)
!	write (*, *) dAdTa
!
!	Contribución asociativa
	dFasdT = 0.D0
	if (NST > 0) then
!
	  assoc_aux(:NST) = matmul( dDeldT(:NST,:NST), sm_Xs(iPhas,:NST) )
	  dFasdT = -dot_product( assoc_aux(:NST),sm_Xs(iPhas,:NST) )/2D0/V

	endif
!	write (*, *) dFasdT
!
!	Contribución repulsiva
	XLAM1=XLAM1*XNTOT        !lambda1
	XLAM2=XLAM2*XNTOT        !lambda2
	XLAM3=XLAM3*XNTOT        !lambda3
	TLAM1=0.D0
	TLAM2=0.D0
	TLAM3=0.D0
	DO 9 I=1,NC

	  DIA=D(I)
	  DIAT=DT(I)
	  XDT=XN(I)*DIAT
	  TLAM1=TLAM1+XDT
	  XDT=XDT*DIA
	  TLAM2=TLAM2+XDT
	  TLAM3=TLAM3+XDT*DIA

9	ENDDO
	PI6=PI/6.D0/V
	Y=1.D0/(1.D0-PI6*XLAM3)
	TLAM2=TLAM2*2.D0
	TLAM3=TLAM3*3.0D0
	DYDT=Y*Y*PI6*TLAM3
	R0=XLAM2/XLAM3
	R1=XLAM1*R0
	R3=R0*R0*XLAM2
	R4=Y*Y-Y-DLOG(Y)
	R6=2.D0*Y-1.D0-1.D0/Y
	R34=Y-1.D0
	R30=TLAM1/XLAM1+TLAM2/XLAM2-TLAM3/XLAM3
	R31=3.D0*TLAM2/XLAM2-2.D0*TLAM3/XLAM3
	!write (*, *) 3.D0*R1*R30*R34+3.D0*R1*DYDT+R3*R6*DYDT+R3*R31*R4+XNTOT/Y*DYDT
!
	DARDT=DADTA + dFasdT + 3.D0*R1*R30*R34+3.D0*R1*DYDT+R3*R6*DYDT+R3*R31*R4+XNTOT/Y*DYDT               !dF/dT
!
	HRES=(-DARDT*T+CZ*XNTOT-XNTOT)*RT    !Hres
!
	RETURN
	ENDsubroutine
!-----------------------------------------------------------------------
!
!
!
subroutine NTLPY_GI (NC, T, Hgi)
!
!	Esta subrutina calcula la entalpí­a del gas ideal a T, para dar entalpí­-
!	as reales complementando con las entalpí­as residuales de NTLPY.
!
!	La subrutina almacena las entalpí­as de GI molares dentro del vector Hgi
!
!	La Hgi se halla por correlación. Los parámetros de la misma se leen en
!	en la subrutina PARMOL.
!
!	Francisco, 14/05/11. Durante la versión 1.3.3.
!
!	           14/07/11. Limpieza, v-1.3.4
!
	implicit real*8 (a-h, o-z)
	integer, parameter :: NCM = 30, NGM = 30, NGAM = 24
	character(10)      :: Cname, GName, GAName
	dimension Hgi(NC)
	common/NAME/CName(NCM), GName(NGM), GAName(NGAM)
	common/GCPROM/PMM(ncm), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM),     &
&	              HHG(NCM)
!
!-----De momento planteo ingresar los parámetros como están en la correlación de Passut-Danner: BTU/(molLB R)
	TRan = T*1.8D0	!K -> ºR (
	Hgi = 0.D0
	Hgi(:NC) = HHA(:NC) + HHB(:NC)*TRan + HHC(:NC)*TRan**2 + HHD(:NC)*TRan**3 + HHE(:NC)*TRan**4 + HHF(:NC)*TRan**5
!
!	Traspaso de unidades de la GCA: atm*cm3/mol-g
	Hgi(:NC) = Hgi(:NC)*PMM(:NC)*1.0412238D4/454.D0	!(BTU/lb)*(10412 atm cm3/BTU)*(PM lb/mol-lb)*(454 mol-g/mol-lb)
!
	return

endsubroutine
!-----------------------------------------------------------------------
!
!
!
subroutine NTRPY_GI (NC, T, Sgi)
!
!	Esta subrutina calcula la entropí­a del gas ideal a T y 1 atm.
!
!	La subrutina almacena las entalpí­as de GI molares dentro del vector Sgi
!
!	La Sgi se halla por correlación. Los parámetros de la misma se leen en
!	en la subrutina PARMOL.
!
!
!	Francisco, 10/03/14. Durante la versión 1.9.26
!
	implicit real(8) (a-h,o-z)
	parameter(NCM = 30)
	character(10) CName
	real(8) Sgi(NC)
	common/NAME/CName(NCM)
	common/GCPROM/PMM(NCM), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM),     &
&	              HHG(NCM)
!
!	Correlación de Passut-Danner: BTU/(mol-lb °R)
	TRan = T*1.8D0	!K -> ºR
	Sgi = 0.D0
	do i = 1, NC
	  
	  Sgi(i) = HHB(i)*dlog(TRan) + 2.D0*HHC(i)*TRan + 1.5D0*HHD(i)*TRan**2 + 4.D0/3.D0*HHE(i)*TRan**3 &
&	           + 1.25D0*HHF(i)*TRan**4 + HHG(i)
!
!	  Traspaso de unidades de la GCA: atm*cm3/(mol-g K)
	  Sgi(i) = Sgi(i)*PMM(i)*1.0412238D4/454.D0*1.8D0
	  !(Sgi BTU/lb °R)*(10412 atm cm3/BTU)*(PM lb/mol-lb)/(454 mol-g/mol-lb)*(1.8 °R/K)
	 
	enddo
!
	return

endsubroutine

!-----------------------------------------------------------------------
!
!	La subrutina PARAGC evalúa los parámetros dependientes de la tem-
!	peratura
!
!
!	d = diámetros de esfera "blanda"
!     g(i,i) y k(i,j) = parámetros energéticos y de interacción
!	Delta(i,j) = fuerza de asociación
!
SUBROUTINE PARAGC (T, NC, NG, NST, NTEMP)
!
	IMPLICIT REAL*8(A-H,O-Z)
	parameter (NCM = 30, NGM = 30, NGAM = 14, NSM=24)
	real*8 kappa
	integer sigma
!
!	Variables especí­ficas del término dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
!
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM)
!
	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!
! 	C = 1D0/(1 - .12D0*dexp(-11D0/9))
	DO 1 I=1,NC

	  ARG=-2.D0*TC(I)/T/3.D0

! 	  ARG=-11*TC(I)/T/9	!mod
	  
	  DARG=DEXP(ARG)
!	  if (NTEMP == 0)GOTO 1
	  if (nTemp > 0) then

	    DT(I)=DC(I)*DARG*.12D0*ARG/T*1.065655D0

! 	    DT(I)=DC(I)*DARG*.12D0*ARG/T*C	!mod

	  endif
	  D(I)=DC(I)*(1.D0-.12D0*DARG)*1.065655D0
	  
! 	  D(I)=DC(I)*(1.D0-.12D0*DARG)*C

1	enddo
	DO 2 I=1,NG

!	  if (T > TSPL(I))GOTO 4
	  if (T < Tspl(i)) then

	    TR=T/TSTR(I)
	    if (nTemp > 0) then
!	    if (NTEMP == 0) GOTO 5

	      DADT(I,I)=GSTR(I)*(G1(I)/TSTR(I)+G2(I)/T)

	    endif
5	    A(I,I)=GSTR(I)*(1.D0+G1(I)*(TR-1.D0)+G2(I)*DLOG(TR))
	    GOTO 2

	  else

	    A(I,I)=GSTR(I)/4.D0/(T/TSPL(I))**EPX(I)
!4	  A(I,I)=GSTR(I)/4.D0/(T/TSPL(I))**EPX(I)
!	  if (NTEMP == 0)GOTO 2
	    if (nTemp > 0) then

	      DADT(I,I)=-EPX(I)*A(I,I)/T

	    endif

	  endif

2	enddo
!2	continue

!	Si hay 1 sólo grupo, no necesita k(i,j)
	
	if (NG >= 2) then
	  DO 3 I = 1, NG - 1

	    NG1=I+1
	    DO 30 J=NG1,NG

	      TR=T/(TSTR(I)+TSTR(J))*2.D0
	      XKIJ1=XKIJ(I,J)*(1.D0+AKIJ(I,J)*DLOG(TR))
	      HELP=DSQRT(A(I,I)*A(J,J))
	      A(I,J)=HELP*XKIJ1
!	      if (NTEMP == 0) GOTO 3
	      if (nTemp > 0) then

	        DKIJ=XKIJ(I,J)*AKIJ(I,J)/T
	        DADT(I,J)=DKIJ*HELP+XKIJ1/2.D0/HELP*(A(I,I)*DADT(J,J)+A(J,J)*DADT(I,I))
	        DADT(J,I)=DADT(I,J)

	      endif
!3	  A(J,I)=A(I,J)
	      A(J,I)=A(I,J)

30	    enddo

3	  enddo
	endif
!
	
	Delta(:NST,:NST) = 0.D0
	if (NST > 0) then
	  do i = 1,NST
	    do j = i,NST

	      Delta(i,j) = kappa(i,j)*(dexp(eps_R(i,j)/T) - 1.D0)
	      Delta(j,i) = Delta(i,j)
	      if (nTemp == 1) then

	        dDeldT(i,j) = -(Delta(i,j) + kappa(i,j))*eps_R(i,j)/T**2
	        dDeldT(j,i) = dDeldT(i,j)

	      endif

	    enddo
	  enddo
	endif
!
      RETURN
ENDsubroutine
!-----------------------------------------------------------------------
!
!	La subrutinas PARGR lee los parámetros grupales dispersivos y aso-
!	ciativos
!
!
!
SUBROUTINE PARGR (NC, NG, NGA, NST, inputFile, outputFile)
!
	IMPLICIT REAL*8(A-H,O-Z)
	parameter (NCM = 30, NGM = 30, NGAM = 14, NSM = 24)

	real(8), dimension(:), allocatable :: X, X1, X2, Xa, Xa1, Xa2
	integer, dimension(:), allocatable :: spaces, MAssoc
	integer, allocatable :: aux_ny(:,:) 
	
	integer NNY(NCM,NGM+1),IGROUP(NGM),INT(2*(NGM+1)), &
     &          IDG(NCM,NGM+1),ID(NGM), strlen
     
!
!	Energí­a y volumen de asociación con 4 subí­ndices: sólo están en
!	esta subrutina.
	dimension enAs(NSM,NGAM,NSM,NGAM),volAs(NSM,NGAM,NSM,NGAM),        &
     &          nyAss(NCM,NGAM)
	integer sigma, outputFile
	real(8) kappa
	
	character(10) :: CName, GName, GAName, string
! 	character(20) :: versGCA, thermoName
	
! 	common /versSUB/ versGCA, thermoName

	COMMON/COORD/ZZ
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)

	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)

	common/NAME/CName(NCM), GName(NGM), GAName(NGAM)
	
! 	versGCA = '1.13.3.2 (2.0RC13)'
! 	thermoName = "GCA-EOS"
!
!	Comprobación de errores:
	if (NC > NCM) then

	  write (*, *) 'Error: número de compuestos máximo sobrepasado'
	  stop

	elseif (NG > NGM) then

	  write (*, *) 'Error: número máximo de grupos atractivos sobrepasado'
	  stop

	elseif (NGA > NGAM) then

	  write (*, *) 'Error: número máximo de grupos asociativos sobrepasado'
	  stop

	endif
!
!  Número de coordinación:
	ZZ=10.D0
!
!  Constante de los gases ideales
	R=82.05D0	!atm cm3/(mol K)

! 	Count number of groups in database:
	NGROUP = 0
	NGroupA = 0
	OPEN (UNIT = 10, FILE = 'FOR010.DAT', STATUS = 'OLD')
	do
	  read (10, *, iostat = ios)
	  if (ios /= 0) exit
	  NGROUP = NGROUP + 1
	enddo
	close (unit = 10)
	allocate (X(NGROUP), X1(NGROUP), X2(NGROUP), spaces(NG))
	
	if (NGA > 0) then
	  
	  OPEN (UNIT = 14, FILE = 'FOR014.DAT', STATUS = 'OLD')
	  do
	    read (14, *, iostat = ios)
	    if (ios /= 0) exit
	    NGROUPA = NGROUPA + 1
	  enddo
	  close (unit = 14)
	  allocate (Xa(NGROUPA), Xa1(NGROUPA), Xa2(NGROUPA), Massoc(NGROUPA))
	  
	endif
!
      read (inputFile, *) NOWNII, NOWNIJ
!
	allocate ( aux_ny(NC, max(NG,NGA)) )
	call group_reading (inputFile, NGROUP, NC, NGM, NG, aux_ny, iGroup)
	ny(:NC,:NG) = aux_ny(:NC,:NG)
!-----------------------------------------------------------------------
!
!	Lectura parametros grupales atractivos puros del banco de datos
	I=0
	J=1
	OPEN (UNIT = 10, FILE = 'FOR010.DAT', STATUS = 'OLD')
!	do while(iGroup(j) > i)

!	  I=I+1
!	  do while(j <= NG)

11	    READ (10, *) (X(K),K=1,5) , k, string
	    i = i + 1
	  if (IGROUP(J) /= I) GOTO 11

	    ID(J)=I
	    Q(J)=X(1)
	    TSTR(J)=X(2)
	    GSTR(J)=X(3)
	    G1(J)=X(4)
	    G2(J)=X(5)
	    GName(j) = string
	    spaces(j) = strlen(GName(j)) !non-blanc spaces in the name
	    J=J+1
	    if (J <= NG) GOTO 11

!	  enddo

!	enddo
	CLOSE (UNIT = 10)
!-----------------------------------------------------------------------
!
!	Lectura parametros grupales atractivos puros ingresados por el usuario
	if (NOWNII > 0) then	!EQ.0) GOTO 12

	  DO 13 I = 1,NOWNII

	    read (inputFile, *) IOWN, (X(J), J = 1,5), string
	    DO 1300 J=1,NG
	      if (IGROUP(J) == iOwn) then	!NE.IOWN) GOTO 13

	        Q(J)=X(1)   		! area q del grupo
	        TSTR(J)=X(2)		! temperatura t* del grupo
	        GSTR(J)=X(3)		! parametro g* del grupo
	        G1(J)=X(4)		! parametro g' del grupo
	        G2(J)=X(5)		! parametro g'' del grupo
	        GName(j) = string	! name
	        spaces(j) = strlen(GName(j))	        

	      endif
1300	    enddo

13	  enddo

	endif
!-----------------------------------------------------------------------
!
!	Halla temperatura de spline con Newton
12	write (outputFile, 105) (J, J = 1,NC)
	write (outputFile, '(X,86("-"),<NC>("-----"))')
	do I = 1,NG

	  if (G1(I) < -1.D-4) then	!GT.-1.D-4) GOTO 18

	    TR=1.D0-1.D0/G1(I)
	    IT=0
	    abs_AX = 1.
	    do while(abs_AX > 1.D-15)

!  Función a ser igualada a 0 = AX(Tspl)
16	      AX=1.D0+G1(I)*(TR-1.D0)+G2(I)*DLOG(TR)-1.D0/4.D0
!  Derivada de dicha función: d(AX)/d(Tspl)
	      DAX=G1(I)+G2(I)/TR
!  Cociente de Newton:
	      DTR=-AX/DAX
	      TR=TR+DTR
	      IT=IT+1
!	    if (TR < 0.D0 .OR. IT > 10) GOTO 18
	      abs_AX=DABS(AX)
!	      if (AX > 1.D-15) GOTO 16
	      if (TR >= 0D0 .AND. it <= 10) then

	        TSPL(I)=TR*TSTR(I)
	        EPX(I)=-4.D0*(G1(I)*TR+G2(I))
!	        GOTO 15

	      else

! Si Ts está dando < 0 o no converge, Ts = 1000.
	        TSPL(I)=1.D3
	        epx(i) = 0
	        exit

	      endif

	    enddo

	  else
!
!  Si g' > 0, no puede haber una Ts > 0 (o difí­cilmente la haya)
18	   TSPL(I)=1.D3
	   epx(i) = 0

	  endif
15	  write (outputFile, 106)ID(I), GName(i), TSTR(I), Q(I), GSTR(I), G1(I), G2(I), Tspl(i), epx(i), (NY(J,I), J = 1, NC)

	enddo
14	continue
!-----------------------------------------------------------------------
105	FORMAT(///, X, 'Pure group atractive parameters:', /, 91X, 'Group configuration per component', /, 3X, &
&	       'No', 3X, "Name", 9X, 'T*(K)', 5X, 'q', 3X, 'g*(atm cm6/mol2)', 5X, 'g`', 8X, 'g"', 6X, "Ts(K)", 4X,        &
&	       "exp", X, <NC>I5)
106	FORMAT(2X, I3, 3X, A10, X, F7.2, F8.4, 4X, F10.1, 2X, 2F10.5, 3X, F7.2, 2X, F5.2, <NC>I5)

!	Lectura parametros de interaccion binaria grupales atractivos del banco de datos
	if (NG > 1) then

	  I=0
	  J=0
!       
!            -------        -------       -------
!	  do while(j <= NG)
        
29	      J=J+1
!       
	  if (J > NG) GOTO 28
        
!	    do while(i <= iGroup(j))	! > i)
        
20	      I=I+1
	      II=NGROUP+1-I
        
	      OPEN (UNIT=11, FILE='FOR011.DAT', STATUS='OLD', access = 'sequential')
	      open (unit = 12, file = 'FOR012.DAT', STATUS='OLD', access = 'sequential')
	      OPEN (UNIT=13, FILE='FOR013.DAT', STATUS='OLD', access = 'sequential')
	      READ(11, *) (X(K),K=1,II)
	      READ(12, *) (X1(K),K=1,II)
	      READ(13, *) (X2(K),K=1,NGROUP)
111	      FORMAT(8F11.5)	    
        
	    if (IGROUP(J) /= I) GOTO 20
!	    enddo
	    XKIJ(J,J)=X(1)
	    AKIJ(J,J)=0.D0
	    L=0
	    K=0
        
27	      L=L+1
        
	    if (L > NG) GOTO 29
        
25	        K=K+1
	      if (IGROUP(L) /= K) GOTO 25
        
	      CABS=DABS(88.888D0-X2(K))
	      if (CABS < 1.D-3) X2(K)=0.D0
	      ALFA(J,L)=X2(K)
!       
!            -------        -------       -------
!       
	    if (K <= I) GOTO 27
        
	    KK=K-I+1
	    CABS=DABS(88.888D0-X(KK))
	    IF (CABS > 1.D-3) GOTO 400
	    XKIJ(J,L)=1.D0
	    AKIJ(J,L)=0.D0
	    write (outputFile, 130) IGROUP(J), GName(j)(:spaces(j)), IGROUP(L), GName(l)(:spaces(l))
130	    FORMAT(/, ' Binary interaction parameters of groups ',I3,' (', A, ")  and ", I3,' (', A, ') are not available in the database.',/ &
&	            , ' Standard (predictive) values will be used unless they are specified by the user.')	  
! 130	    FORMAT(/, ' Los parámetros de interacción binaria entre los grupos ',I3,' (', A, ")  y ", I3,' (', A, ') no están disponibles.',/ &
! &	            , ' Se usarán valores estándar a menos que estos sean ingresador por el usuario')	  
	    GOTO 401
400	    XKIJ(J,L)=X(KK)
	    AKIJ(J,L)=X1(KK)
401	    AKIJ(L,J)=AKIJ(J,L)
	    XKIJ(L,J)=XKIJ(J,L)
	    GOTO 27
        
28	  CONTINUE
	  CLOSE (UNIT=11)
	  CLOSE (UNIT=12)
	  CLOSE (UNIT=13)
!-----  ------------------------------------------------------------------
!       
!	  Lectura parametros de interaccion binaria grupales atractivos
!       ingresados por el usuario
        if (NOWNIJ > 0) then	!EQ.0) GOTO 33
	    DO 21 I=1,NOWNIJ
        
	      read (inputFile, *) IDI,IDJ,XK,AK,ALF1,ALF2
	      DO 22 J=1,NG
        
	        if (IGROUP(J) == IDI) INI=J
	        if (IGROUP(J) == IDJ) INJ=J
        
22	      enddo
	      XKIJ(INI,INJ)=XK			! (k*)ij
	      XKIJ(INJ,INI)=XK			! (k*)ji
	      AKIJ(INI,INJ)=AK			! (k')ij
	      AKIJ(INJ,INI)=AK			! (k')ji
	      ALFA(INI,INJ)=ALF1		! (alfa)ij
	      ALFA(INJ,INI)=ALF2		! (alfa)ji
        
21	    enddo
	  endif
!-----  ------------------------------------------------------------------
!       
!	  Impresión de parámetros binarios
33	  write (outputFile, 115) GName(:NG)
115	  FORMAT (/,1X,'Binary interaction parameter at T* matrix, [k*]', //, 18X, <NG>A10)
116	  FORMAT (3X, A10, 2X, '|',<NG>(F8.4, 2X),'|')
	  do i = 1, NG
        
	    write (outputFile, 116) GName(i), XKIJ(I,:I)
        
	  enddo
	  write (outputFile, 117) GName(:NG)
117	  FORMAT(/,1X,'Binary interaction parameter temperature dependence matrix, [k`]', //, 18X, <NG>A10)	
	  do i = 1,NG
        
	    write (outputFile, 116)  GName(i), AKIJ(I,:I)
        
	  enddo
	  write (outputFile, 118) GName(:NG)
118	  FORMAT(/,1X,'Binary dumping factors matrix, [alpha]', //, 18X, <NG>A10)	
	  write (outputFile, '(3X, A10, 2X, "| -------  ", <NG-1>(F8.4, 2X), "|")') GName(1), ALFA(1,2:NG)
	  DO 23 I = 2, NG - 1
        
	    write (outputFile, 119)  GName(i), ALFA(I,:I-1), ALFA(I,I+1:NG)
        
23	  enddo
119	  format (3X, A10, 2X, '|', <I-1>(F8.4, 2X), " -------  ", <NG-I>(F8.4, 2X), '|')
	  write (outputFile, '(3X, A10, 2X, "|", <NG-1>(F8.4, 2X), " -------  |")') GName(NG), ALFA(NG,1:NG-1)
	  
	endif
!-----------------------------------------------------------------------
!
	NST = 0
	if (NGA > 0) then
	  
	  deallocate (spaces)
	  allocate (spaces(NGA))
	  iGroup = 0
	  read (inputFile, *) NOWNII, NOWNIJ ! The same but now for the assoc contr.
!
	  call group_reading (inputFile, NGROUPA, NC, NGAM, NGA, aux_ny, iGroup)
	  nyAss(:NC,:NGA) = aux_ny(:NC,:NGA)
	  
	  open (unit = 14, file = "FOR014.DAT", access = "sequential", status = "old")
	  open (unit = 15, file = "FOR015.DAT", access = "sequential", status = "old")
	  X(1) = maxval(iGroup(:NGA))
	  maxIndex = X(1)
	  i2 = 0
! 	  Reading of the sites per group in the database. "i" is a group desired for the user, and "i2" represents
!	  a database group.
	  i = 1
	  do i2 = 1, NGROUPA !(databse group)
	      
! 	    if (i > NGA) exit
	    
	    read (14, *) Xa(i2), string !it's necessary to read _all_ de database, in order of being being able to 
	                                !run over each row of the interaction matrix.
	    if (i2 == iGroup(i)) then
	    
	      MAssoc(i) = Xa(i2)
	      GAName(i) = string
	      spaces(i) = strlen(GAName(i))
	      i = i + 1
	      
	    endif
	    
	  enddo
	  write (outputFile, 217) NGA
217	  FORMAT(//,X,"Number of associating groups =",I3)
!
	  write (outputFile, 218)
218	  FORMAT(/,X,'Associating group configuration ',/,34X,'Compound No')

	  write (outputFile, '("   No.  Name", 7X,  "No of sites ",<NC>(I4,X))')((I), I=1,NC)
!
	  write (outputFile, '(X,29("-"),<NC>(5("-")))')
	  DO I=1,NGA

	    write (outputFile, 219) iGroup(i), GAName(i), MAssoc(I), (nyass(j,i), J=1,NC)
219	    FORMAT(2X, I3, 3X, A10, 4X, I3, 6X, <NC>(I4,X))

	  END DO	  
	  
! 	  Reading of the binary association database
	  i = 1
	  j = 1
	  do i2 = 1, NGROUPA !(databse group)

	    do k = 1, Xa(i2)
	
	      j = i
	      do j2 = i2, NGROUPA !(databse group)
	    
	        if (j2 == i2) then
	
	          do l = k, Xa(j2)
	      
	            read (15, *) Xa1(:2)
	            if (i2 == iGroup(i) .AND. j2 == iGroup(j)) then
	
	               if (dabs(Xa1(1) - 88.888D0) < 1D-3) then                        
	                 write (outputFile, 130) IGROUP(i), GAName(i)(:spaces(i)), IGROUP(j), GAName(j)(:spaces(j))
	                 Xa1(1) = 0
	                 Xa1(2) = 0
	               endif
	               enAs(k,i,l,j) = Xa1(1)
	               enAs(l,j,k,i) = Xa1(1)
	               volAs(k,i,l,j) = Xa1(2)
	               volAs(l,j,k,i) = Xa1(2)
	                
	            endif
	            
	          enddo
	          
	        else
	
	          do l = 1, Xa(j2)
	      
	            read (15, *) Xa1(:2)
	            if (i2 == iGroup(i) .AND. j2 == iGroup(j)) then
	
	               if (dabs(Xa1(1) - 88.888D0) < 1D-3) then                        
	                 write (outputFile, 130) IGROUP(i), GAName(i)(:spaces(i)), IGROUP(j), GAName(j)(:spaces(j))
	                 Xa1(1) = 0
	                 Xa1(2) = 0
	               endif
	               enAs(k,i,l,j) = Xa1(1)
	               enAs(l,j,k,i) = Xa1(1)
	               volAs(k,i,l,j) = Xa1(2)
	               volAs(l,j,k,i) = Xa1(2)   
	                
	            endif
	            
	          enddo	          
	          
	        endif
	        if (j2 == iGroup(j)) then
		          
	          j = j + 1 !time to assign in user matrix
	          
	        endif	            
		  
	      enddo
	      
	    enddo
	    if (i2 == iGroup(i)) then
		
	      i = i + 1 !time to assign in user matrix
	      j = i
	      
	    endif
	    if (i2 == maxIndex) exit !it's not necessary to continue reading
	    
	  enddo
	  close (unit = 14); close (unit = 15)                     
!
!  Read user specified association parameters	
	  do i2 = 1, NOWNII
	    
	    read (inputFile, *) iOwn, j, string
	    do i = 1, NGA
	    
	      if (iGroup(i) == iOwn) then
	  
	        MAssoc(i) = j
	        GAName(i) = string
	        
	      endif
	  
	    enddo
	    
	  enddo
	  do i2 = 1, NOWNIJ
	    
	    read (inputFile, *) k, iOwn, l, jOwn, Xa1(:2)
	    do i = 1, NGA
	      if (iGroup(i) == iOwn) then
		  
	        do j = i, NGA
	          if (iGroup(j) == jOwn) then
                
	            enAs(k,i,l,j) = Xa1(1)
	            enAs(l,j,k,i) = Xa1(1)
	            volAs(k,i,l,j) = Xa1(2)
	            volAs(l,j,k,i) = Xa1(2)  
	          
	          endif
	        enddo
	        
	      endif
	    enddo
	    
	  enddo

!  Conversion de variables a notacion por sitio "m" (y no, por sitio
!  "k" en el grupo "i")
	  m1 = 0
	  m2 = 0
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      m2 = m1
	      m1 = m1 + 1      	!Defino actual sitio "k" como sitio "m1"
	      do j = i,NGA      !esta definicion solo es valida gracias a que hasta ahora, solo ha habido
	        if (i == j) then  !un tipo de sitio en cada grupo.
	          do l = k,Massoc(j)

	            m2 = m2 + 1 !defino actual sitio "l" como sitio "m2"
	            eps_R(m1,m2) = enAs(k,i,l,j)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = volAs(k,i,l,j)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        else
	          do l = 1,Massoc(j)

	            m2 = m2 + 1  !defino actual sitio "l" como sitio "m2"
	            eps_R(m1,m2) = enAs(k,i,l,j)
	            eps_R(m2,m1) = eps_R(m1,m2)
	            kappa(m1,m2) = volAs(k,i,l,j)
	            kappa(m2,m1) = kappa(m1,m2)

	          enddo
	        endif
	      enddo

	    enddo
	  enddo
	  NST = sum(Massoc(:NGA))
!
        m1 = 0
!  Reasignación en notación 1..NST:
	  do i = 1,NC

	    m1 = 0
	    do j = 1,NGA
	      do k = 1, MAssoc(j)

	        m1 = m1 + 1
	        sigma(m1,i) = nyAss(i,j)

	      enddo
	    enddo

	  enddo
!
!  Impresion parametros asociativos
!
	  write (outputFile, '(/, X, "Associating energy matrix, [epsilon/R] (K)", /)')

	  write (outputFile, 220) ((GAName(i), k = 1, MAssoc(i)), i = 1, NGA)
	  write (outputFile, 221) ((k, k = 1, MAssoc(i)), i = 1, NGA)	  
220	  FORMAT (23X, <NST>(A10, 1X) )
221	  format (23X, <NST>(I1, 10X))

!  Escritura de parámetros energéticos de asociación en forma de matriz
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write (outputFile, 222) GAName(i), k, ((enAs(k,i,l,j), l = 1, MAssoc(j)), j = 1,i)

	    enddo
	  enddo
222	  format(3X, A10, X, I2,' |',<NST>(2X,F7.1,2X),'|')
!  Escritura del volumen de asociación
	  write (outputFile, '(/," Associating volume matrix, [kappa] (cm3/mol)", /)')
	  write (outputFile, 220) ((GAName(i), k = 1, MAssoc(i)), i = 1, NGA)
	  write (outputFile, 221) ((k, k = 1, MAssoc(i)), i = 1, NGA)	  
	  do i = 1,NGA
	    do k = 1,Massoc(i)

	      write (outputFile, 223) GAName(i), k, ((volAs(k,i,l,j), l = 1, Massoc(j)), j = 1, i)

	    enddo
	  enddo
223	  format(3X, A10, X, I2,' |',<NST>(X,G10.4),' |')

	endif
!-----------------------------------------------------------------------


100	FORMAT(20I3)
101	FORMAT(5F10.2)
102	FORMAT(I3,5F10.2)

! ! 119	FORMAT(2I3, 4F10.2)
120	FORMAT(1X,90("-"))

	RETURN
ENDsubroutine
!------------------------------------------------------------------------------
!  
!  This subroutine reads group configuration per compound and sorts them ascending.
!
!  It works for any "contribution".
!
subroutine group_reading (inputFile, database_size, NC, NGM, NG, ny, iGroup)

	implicit none
	integer                                :: i, ii, ij, iTop, j, j1, j2, jj
	
	integer, intent(in)                    :: database_size, inputFile, NC, NG, NGM
	integer, dimension(2*(NG+1))           :: idg_ny
	integer, dimension(NC,NGM+1)           :: idG, NNy
	
	integer, dimension(NGM), intent(out)   :: iGroup
	integer, dimension(NC,NG), intent(out) :: ny
	
	
	idg_ny = 0
	iGroup(:NGM) = database_size + 1
	DO 1 I = 1, NC

	  read (inputFile, *) (idg_ny(J), J = 1, 2*NG)
	  DO 2 J = 1, NG + 1

	    IDG(I,J) = idg_ny(2*J - 1)  !identified which is the group J in compound I 
	    NNY(I,J) = idg_ny(2*J)      !number of groups J in compound I. Note that variable NNY is different fron NY, but related.

2	  enddo

	  II = 1
	  IJ = 1
	  
!  Sorts ascending the attractive group ID numbers.

7	  if (idG(i,ij) /= 0) then !idG(i,ij) = 0 at the end of the read.

	      if (idG(i,ij) <= iGroup(ii)) then

	        if (idG(i,ij) /= iGroup(ii)) then

	          iTop = NGM - II
	          if (iTop /= 0) then

	            DO 6 J = 1, ITOP

	              J1 = NGM - J ! exchanges iGroup(j) with j+1
	              J2 = J1+1
	              IGROUP(J2)=IGROUP(J1)

6	            enddo

	          endif
8	          IGROUP(II)=IDG(I,IJ)

	        endif
5	        IJ=1+IJ

	      endif
4	      II = 1 + II

	    if (II <= NGM) goto 7

	  endif

1	enddo
!
!  Generate compound group contribution for the NG groups.
!  Result is stored in NY matrix.
!
	DO 10 I = 1, NC
	  DO 31 J = 1, NG

	    NY(I,J) = 0
	    JJ = 1
9	    if (idG(i,jj) /= 0) then

	      do while (jj <= NG)

	        if (idG(i,jj) == iGroup(j)) then

!	          Number of j groups in the I-th compound
	          NY(I,J) = NNY(I,JJ)
	          exit

	        else

32	          JJ = JJ + 1

	        endif

	      enddo

	    endif

31	  enddo
10 	enddo
!-----------------------------------------------------------------------
	return
endsubroutine group_reading
!-----------------------------------------------------------------------
!
!	La subrutinas PARMOL lee los parámetros moleculares de cada compo-
!	nente:
!	Tc, en K
!	Pc, en atm
!	omega (factor acéntrico
!	Tsat, en K
!	Psat, en atm
!	dc, en cm/mol^1/3
!	PM (peso molec), en g/mol o lb/mol-lb
!	c (Peneloux), en cm3/mol
!	A, B, C, D, E, F (coef de Passut&Danner para Hgi, en BTU/mol-lb/ºR
!
!
SUBROUTINE PARMOL (NC, NG, NST, inputFile, outputFile)
!
	IMPLICIT REAL*8 (A-H, O-Z)
	parameter (NCM = 30, NGM = 30, NGAM = 14, NSM = 24)
	character*10 :: CName, GName, GAName
!
	DIMENSION TSAT(NCM), PS(NCM), Z(2), NYOLD(NGM), EXPOL(NGM)
	DIMENSION ID(NGM),GSTOL(NGM),G1OL(NGM),G2OL(NGM),TSTROL(NGM),      &
     &          TSPLOL(NGM),EPXOL(NGM),QOL(NGM),XKIJOL(NGM,NGM),         &
     &          XK1OL(NGM,NGM),ALFOL(NGM,NGM)
!
	common/versSUB/versGCA
!
	COMMON/GROUP1/GSTR(NGM),G1(NGM),G2(NGM),TSTR(NGM),TSPL(NGM),       &
     &              XKIJ(NGM,NGM),AKIJ(NGM,NGM),EPX(NGM)
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)
!
	common/GrupAs1/kappa(NSM,NSM),eps_R(NSM,NSM)
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM), omega(NCM)
!
	common/NAME/CName(NCM), GName(NGM), GAName(NGAM)
	common/GCPROM/pmm(ncm),Pen(NCM),HHA(ncm),HHB(ncm),hhc(ncm),HHD(ncm),HHE(NCM),HHF(NCM),     &
&	              HHG(NCM)
	common/ENTALP/iEntalp

	dimension eps_R_bak(NST,NST)
	real*8 kappa_bak(NST,NST),kappa
	integer sigma, sigma_bak(NST,NC), outputFile
!
!  Iniciación de constantes de Antoine
	HA(:NC) = 0D0
	HB(:NC) = 0D0
!
!	write (outputFile, '(//,"Diámetros críticos",/)')
	read (inputFile, *) nUnit    !nUnit 0 = unidades en MPa, nUnit=1 en atm
	DO I=1,NC

	  write (*, *) i
	  read (inputFile, *) CNAME(I)
	  read (inputFile, *) TC(I), PC(I), OMEGA(I), TSAT(I), PS(I), DC(I)
	  read (inputFile, *) PMM(I), Pen(i), HHA(I), HHB(I), HHC(I), HHD(I), HHE(i), HHF(i), HHG(i)
	  if (nUnit == 0) then

	    PC(I) = PC(I)/0.101325D0
	    Ps(i) = Ps(i)/0.101325D0

	  endif
	  if (dc(I) == 0) then

	    DC(I)=7.337449944D0*TC(I)/PC(I)
	    DC(I)=DC(I)**.3333333333D0

	  end if
!	  write (outputFile, *) CName(i),' dc = ',dc(i),'cm/mol^(1/3)'

	enddo
1	continue
!
!  Guarda la información del componente "1", porque para ajustar el dc(i)
!  sobrescribe esta información, debido al que'l traspase de la misma se
!  hace mediantes "COMMON"s y no mediante entrada de subrutinas.
	DO 8 I=1,NG

	  NYOLD(I)=NY(1,I)

8	enddo
!
!
!  Almacenamiento de variables grupales asociativas
!  Entiendo que es un derroche talvez, almacenar TODO el paquete de
!  variables, pero también es muy sencillo... y seguro.
!			Francisco, 18/07/2011, v-1.8.5
	
	sigma_bak(:NST,:NC) = sigma(:NST,:NC)
	
	kappa_bak(:NST,:NST) = kappa(:NST,:NST)
	
	eps_R_bak(:NST,:NST) = eps_R(:NST,:NST)
!
!  Almacenamiento de variables moleculares (término repulsivo)
	TCOLD=TC(1)
	PCOLD=PC(1)
	DCOLD=DC(1)
	DO 3 I=1,NC

! 	  write (*, *) 'comp=', i
!  Ajusta el dc(i) para reproducir el dato de presión de vapor
	  NG1=0
	  if (I == 2) DCOLD=DC(1)
	  DO 400 J=1,NG

!  Almacenamiento de variables grupales atractivas
	    if (NY(I,J) > 0) then	!EQ.0) GOTO 400

	      NG1=NG1+1
	      NY(1,NG1)=NY(I,J)
	      ID(NG1)=J
	      GSTOL(NG1)=GSTR(NG1)
	      G1OL(NG1)=G1(NG1)
	      G2OL(NG1)=G2(NG1)
	      TSTROL(NG1)=TSTR(NG1)
	      TSPLOL(NG1)=TSPL(NG1)
	      EXPOL(NG1)=EPX(NG1)
	      QOL(NG1)=Q(NG1)
	      GSTR(NG1)=GSTR(J)
	      G1(NG1)=G1(J)
	      G2(NG1)=G2(J)
	      TSTR(NG1)=TSTR(J)
	      TSPL(NG1)=TSPL(J)
	      EPX(NG1)=EPX(J)
	      Q(NG1)=Q(J)

	    endif

400	  enddo
        DO 401 J=1,NG1

	    IJ=ID(J)
	    DO 4010 K=1,NG1

	      IK=ID(K)
	      XKIJOL(J,K)=XKIJ(J,K)
	      XKIJ(J,K)=XKIJ(IJ,IK)
	      XK1OL(J,K)=AKIJ(J,K)
	      AKIJ(J,K)=AKIJ(IJ,IK)
	      ALFOL(J,K)=ALFA(J,K)
	      ALFA(J,K)=ALFA(IJ,IK)

4010	    enddo

401	  enddo
!
!
!  Traspase a la posición "1" de parámetros asociativos
	  m = 0
	  do k = 1,NST

	    if (sigma_bak(k,i) /= 0) then

	      n = m
	      m = m + 1
	      sigma(m,1) = sigma_bak(k,i)

	      do l = k,NST

	        if (sigma_bak(l,i) /= 0) then

	          n = n + 1
	          kappa(m,n) = kappa_bak(k,l)
	          kappa(n,m) = kappa(m,n)
	          eps_R(m,n) = eps_R_bak(k,l)
	          eps_R(n,m) = eps_R(m,n)

	        endif

	      enddo

	    endif

	  enddo
	  NST1 = m !?????   Hasta acá... 01:00 hrs. v1.8.2
!
!
!  Carga las propiedades del compuesto "i" en la primera posición.
!  Parámetros moleculares y repulsivo:
	  TC(1)=TC(I)
	  PC(1)=PC(I)
	  DC(1)=DC(I)
	  if ((TSAT(I) > 1.D-1) .OR. (omega(i) > 1.d-8)) then ! GOTO 7

	    T=TSAT(I)
	    P1=PS(I)
	    PFIT=P1
	    if (Tsat(i) < 1D-7) then
!	    if (TSAT(I) > 1.D-1) GOTO 7


	      T=TC(I)*.7D0
	      P1=PC(I)*DEXP(-(1.D0+OMEGA(I))*2.3026D0)
	      PFIT=P1

	    endif
7	    CALL PARAGC(T,1,NG1,NST1,0)
	    IGP=1
	    IGZ=0
!---------Lazo de ajuste del dc----------------------------------------------------
	    fCrt = 1.0
! 	    write(*,*)i
	    iter_dc_fit = 0
	    do while(fCrt >= 1.E-8)

	      iter_dc_fit = iter_dc_fit + 1
5	      D1=D(1)
	      
	      CALL PSAT(T,P1,Z,1,NG1,NST1,1,IGP,IGZ,IER)
	      if (iter_dc_fit == 3) IGZ=5
	      F1=P1-PFIT
	      FCRT=DABS(F1)

!	    if (FCRT < 1.D-8) GOTO 4

	      P2=P1
	      D(1)=D1+1.D-4
	      CALL PSAT(T,P2,Z,1,NG1,NST1,1,IGP,IGZ,IER)
	      F2=P2-PFIT
	      DFDD=(F2-F1)/1.D-4
	      DPDD=(P2-P1)/1.D-4
! 	      write (*, *) d(1),dFdd
	      if (DFDD == 0) then
	        write (*, *) "comp = ", i
	        stop
	      endif
	      DEL=-F1/DFDD
	      DLC=DEL/D(1)
	      if (DLC > .05D0) DEL=.05D0*D(1)
	      if (DLC < -.05D0) DEL=-.05D0*D(1)
	      D(1)=D1+DEL
	      P2=P1+DEL*DPDD
	      Z(2)=Z(2)/P1*P2
	      P1=P2

!	    GOTO 5
	    enddo
!-----------------------------------------------------------------------
!
!	    Nuevo valor del diámetro crí­tico del compuesto "i":
4	    DC(I)=D(1)/(1.D0-.12D0*DEXP(-2.D0*TC(I)/3.D0/T))/1.065655D0
	    TIN=TSAT(I)
	    PIN=PS(I)

	  endif
	  if (TSAT(I) < 1.D-1) TIN=.7D0*TC(I)
	  if (TSAT(I) < 1.D-1) PIN=PC(I)*DEXP(-(1.D0+OMEGA(I))*2.3026D0)
	  HB(I)=DLOG(PC(I)/PIN)/(1.D0/TC(I)-1.D0/TIN)
	  HA(I)=DLOG(PC(I))-HB(I)/TC(I)
!
!	  Reasignación de parámetros grupales atractivos:
	  DO 402 J=1,NG1

	    GSTR(J)=GSTOL(J)
	    G1(J)=G1OL(J)
	    G2(J)=G2OL(J)
	    TSTR(J)=TSTROL(J)
	    TSPL(J)=TSPLOL(J)
	    EPX(J)=EXPOL(J)
	    Q(J)=QOL(J)
	    DO 4020 K=1,NG1

	      XKIJ(J,K)=XKIJOL(J,K)
	      AKIJ(J,K)=XK1OL(J,K)
	      ALFA(J,K)=ALFOL(J,K)

4020	    enddo

402	  enddo

3	enddo	!CONTINUE
!
!  Reasignación de configuracion asociativa y parámetros
	
	sigma(:NST,:NC) = sigma_bak(:NST,:NC)
	
	kappa(:NST,:NST) = kappa_bak(:NST,:NST)
	
	eps_R(:NST,:NST) = eps_R_bak(:NST,:NST)
!
!  Reasignación de parámetros moleculares del comp. 1:
	if (NC > 1) then
	  
	  TC(1)=TCOLD
	  PC(1)=PCOLD
	  DC(1)=DCOLD	!if NC = 1, this would overwrites dc previous fit
	  
	endif

!  Reasignación de configuración grupal dispersiva del comp. 1
!      DO 9 I=1,NG
      do i = 1,NG
!    9 NY(1,I)=NYOLD(I)
	  ny(1,i) = nyOld(i)

	enddo

      if (nUnit == 1) then

	  write (outputFile, 10)

	elseif (nUnit == 0) then

	  write (outputFile, 12)

	endif
	write (outputFile, 13)
	do i = 1,NC

	  write (outputFile, 11) CName(i),TC(I),PC(I),OMEGA(I),TSAT(I),PS(I),DC(I),PMM(i),Pen(i)

	enddo
!
!	Constantes de Passut-Danner: impresiónn sólo si al menos se le ingresa valores a un compuesto:
	iEntalp = NC
	do i = 1,NC
	  if ((HHA(i) == 0) .AND. (HHB(i) == 0) .AND. (HHC(i) == 0) .AND. (HHD(i) == 0) .AND. (HHE(i) == 0) .AND. (HHF(i) == 0)) then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	if (iEntalp /= 0) then

	  write (outputFile, '(/," Passut-Danner constant for ideal gas enthalpy calculation:",/,"   Compound   A(BTU/lb)  B(BTU/(lb ºR))  ", &
     &             "C(BTU/(lb ºR2)·1E3  D(BTU/(lb ºR3))·1E6  E(BTU/(lb ºR4))·1E10  F(BTU/(lb ºR5))·1E14     G(BTU/(lb ºR))",/,148("-"))')
	  do i = 1,NC

	    write (outputFile, 14)CName(i),HHA(i),HHB(i),1.D3*HHC(i),1.D6*HHD(i),1D10*HHE(i),1D14*HHF(i),HHG(i)

	  enddo

	endif

   91 FORMAT(8F10.4)
   10 FORMAT(//,X,'Pure component properties',//,3X,'Compound',6X,'Tc(K)',1X,'Pc(atm)',3X,'omega',2X,'Tsat(K)',2X,'Psat(atm)',2X, &
     &       'dc(cm/mol^1/3)',3X,"M(g/mol)",3X,"c(cm3/mol)")
   11 FORMAT(4X,A10,1X,2F7.1,3X,F6.4,F9.2,2X,F8.4,5X,F8.4,7X,F6.2,8X,F6.2)
   12 FORMAT(//,2X,'Pure component properties',/,3X,'Nro.',2X,'Tc(K)',2X,'Pc(MPa)',3X,'omega',2X,'Tsat(K)',2X,'Psat(MPa)',2X,  &
     &       'dc(cm/mol^1/3)',X,"M(g/mol)",X,"c(cm3/mol)")
13	format("  ",100("-"))
14	format(3X,A10,X,F9.6,3X,F9.6,10X,F9.6,12X,F9.6,13X,F9.6,13X,F9.6,12X,F9.6)
	RETURN
ENDSUBROUTINE PARMOL
!-----------------------------------------------------------------------
!
!
!
SUBROUTINE PSAT (T, PS, Z, NC, NG, NST, IC, IGP, IGZ, IER)
!
!  Esta subrutina provee las presiones de vapor de compuestos puros.
!  Los volumenes de las fases son obtenidos llamando a ZMAX.
!
!  Mensajes de error (Valores de iEr)
!      1    : ningún error
!     -1    : solución trivial (vL = vV)
!
	IMPLICIT REAL*8(A-H,O-Z)
	parameter (NCM = 30, NGM = 30, NSM = 24)
	integer sigma
	real*8 LAMBDA
	DIMENSION X(NCM),Z(2)
	COMMON/SAT/PHI3
	COMMON/COORD/ZZ
	COMMON/EXTREM/IEXT(2)
!
!  Variables especí­ficas del término dispersivo
	COMMON/GROUP2/Q(NGM),A(NGM,NGM),DADT(NGM,NGM),ALFA(NGM,NGM),R,     &
     &              NY(NCM,NGM)

!  Versión asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!

!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),rNyT(NGM),QT,    &
     &           TETA(NGM),HELP12(2,NGM),HELP15(2,NGM),HELP7(2,NGM),     &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),XLAM1,XLAM2,XLAM3

!	Versión asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM), omega(NCM)
!
!
	IGOLD=IGZ
	if (IGZ > 1) IGZ=-1
	EPSP = 1.D-10
	PI = 3.1415926536D0
	IER = 1
	x(:NC) = 0D0
	X(IC) = 1.D0
	if (IGP == 0) PS = DEXP(HA(IC) + HB(IC)/T)
!
!-----------------------------------------------------------------------
!
!  Cálculo de iterativo de la presion de vapor.
!  Utiliza el criterio de igualdad de isofugacidad
	Pcrt = 1.
    dxi = 0
!	do while(PCRT > EPSP)

3	  CALL ZMAX(Z,1,T,PS,X,IGZ,NC,NG,NST,IF)
	  if (IEXT(1) /= 0) then	!EQ.0) GOTO 20

	    PS=PS*.75D0
	    GOTO 3

	  endif
20	  if (IEXT(2) /= 0) then	!EQ.0) GOTO 21

	    PS=PS*1.25D0
	    GOTO 3

	  endif
21	  ZCR=DABS(Z(1)-Z(2))
!
!  Fase supercritica:
	  if (ZCR < 1.D-6) IER=-1
	  if (IER < 0) GOTO 99
!	  if (iEr >= 0) then
!
!  Contribución repulsiva de la energia de Helm. (fase vapor y lí­quida)
	    XSIV = PS/T/R/Z(1)*PI/6.D0*PHI3
	    XSIL = PS/T/R/Z(2)*PI/6.D0*PHI3
	    CSV=3.D0*XSIV/(1.D0-XSIV)+XSIV/(1.D0-XSIV)**2
	    CSL=3.D0*XSIL/(1.D0-XSIL)+XSIL/(1.D0-XSIL)**2
!
!  Contribución atractiva de la energia de Helm. (fase vapor y lí­quida)
	    ATTV=0.D0
	    ATTL=0.D0
	    DO 2 K=1,NG

	      XNY=DFLOAT(NY(IC,K))
	      if (XNY >= 1D-10) then	!LT.1.D-10) GOTO 2

	        ATTV=ATTV+XNY*Q(K)*HELP2(1,K)
	        ATTL=ATTL+XNY*Q(K)*HELP2(2,K)

	      endif

2	    enddo
	    ATTV=ATTV*ZZ/2.D0
	    ATTL=ATTL*ZZ/2.D0
!
!  Contribución asociativa de la energia de Helm. (fase vapor y lí­quida)
	    AasocL = 0D0
	    AasocV = 0D0
	    if (NST > 0) then

!	      Vapor:
	      AasocV =  dot_product(sm(:NST),dlog(Xs(1,:NST))) + (Sum_S - sum(sm_Xs(1,:NST)))*.5D0
!	      Lí­quido:
	      AasocL =  dot_product(sm(:NST),dlog(Xs(2,:NST))) + (Sum_S - sum(sm_Xs(2,:NST)))*.5D0

	    endif
!
!	    F = CSV - ATTV - CSL + ATTL - DLOG(Z(1)/Z(2)) - AasocL + AasocV
	    F = CSV - ATTV - CSL + ATTL + AasocV - AasocL + Zcr - DLOG(Z(1)/Z(2)) !isofugacity = ln(phiV) - ln(phiL)
!	    PNY=-F/(Z(1)-Z(2))*PS !This is a Newton update in P (Pny = Ps - F/(dF/dP)).
	    xi = dlog(Ps)
	    dxi_old = dxi
	    dxi = F/Zcr !Newto step in xi = ln(P)
	    test = dabs(dxi + dxi_old)/EPSP
	    if (test < 10) dxi = dxi/2
	    xi_new = xi - min(max(dxi,-1D0),1d0)
	    Pny = dexp(xi_new)
	    PCRT=dABS(PNY-PS)
!	    PLIM=(PS-PNY)/PS
!	    if (PLIM > .1D0) PNY=.9D0*PS
!	    if (PLIM < -.1D0) PNY=1.1D0*PS
	    Z(2)=Z(2)*PNY/PS
	    PS=PNY

	if (PCRT > EPSP) GOTO 3
!	  else
!
!  Si la fase es supercritica sale del lazo.
!	    exit

!	  endif
99	  CONTINUE

!	enddo
	IGZ=IGOLD
	RETURN

ENDsubroutine PSAT
!-----------------------------------------------------------------------
!
!
SUBROUTINE ZMAX (Z, iTyp, T, P, XN, iGues, NC, NG, NST, iC)
!
! 	Esta subrutina encuentra el(los) factor(es) de compresibilidad del
!	gas y/o la fase lí­quida a determinada composición, temperatura y
!	presión.
!	El(los) factor(es) de compresibilidad es(son) encontrado por solu-
!	ción iterativa de
!
!           P(EXP)=-RT(D/DV(Ares/RT)-rho)
!
!     La solución es devuelta en Z.
!
!     Si iGuez > ó = 1, el usuario debe proporcionar un buen valor ini-
!	cial del factor de compresibilidad deseado (S). De lo contrario se
!	genera automáticamente otra estimación inicial.
!
!	Las opciones siguientes especifican la solución deseada
!
!	Opciones
!     iTyp   iGues
!      -1  cualquiera   Fase gaseosa
!       1    0 OR 1     Fases lí­quida y gaseosa
!       1      >1       Fase lí­quida
!       0       0       Fases lí­quida y gaseosa verificando energí­a libre
!                       de Gibbs mí­nima
!
!     Mensajes de error (valores de iC)
!      1    : fase correcta
!     -1    : fase incorrecta
!     Si iType = 0 el valor devuelto es el que de menor energí­a libre de
!	Gibbs
!
!-----------------------------------------------------------------------
	IMPLICIT REAL*8 (A-H, O-Z)

	parameter (NCM = 30, NGM = 30, NSM = 24)
	parameter (maxit=200)
!
	DIMENSION Z(2), XN(NCM), X(NCM)
	DIMENSION AA(NGM,NGM), IER1(2)

	integer, dimension(:),   allocatable :: pivot_vector
	real(8), dimension(:,:), allocatable :: Delta_aux, H_aux
	real(8), dimension(:),   allocatable :: Xs_aux, sm_aux

	logical calc


!
!  Variables auxiliares de asociación
	dimension PSI_V(NSM)
	integer sigma
!
!  Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM), HELP4(2,NGM), HELP5(2,NGM), rNyT(NGM), QT,    &
     &           TETA(NGM), HELP12(2,NGM), HELP15(2,NGM), HELP7(2,NGM),     &
     &           E(2,NGM,NGM), PREP(2), DPDV(2), XLAM1, XLAM2, XLAM3
!
!  Variables especí­ficas del término dispersivo
	COMMON/GROUP2/Q(NGM), A(NGM,NGM), DADT(NGM,NGM), ALFA(NGM,NGM), R,     &
     &              NY(NCM,NGM)
!
!  Propiedades moleculares
	COMMON/MOL/DC(NCM), D(NCM), DT(NCM), HA(NCM), HB(NCM)
	COMMON/CRIT/TC(NCM), PC(NCM), omega(NCM)
!
	COMMON/COORD/ZZ
	COMMON/SAT/PHI3
	COMMON/EXTREM/IEXT(2)
!
!  Versión asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM), sm(NSM), dXs_dV(2,NSM)

!  Versión asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM), dDeldT(NSM,NSM)

	common/GrupAs3/sm_Xs(2,NSM), Sum_S, dFVas(2)
      common/GrupAs5/H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
!
	common/itervolumen/nroIter(2)
	common/GCPROM/PMM(NCM), Pen(NCM), HHA(NCM), HHB(NCM), HHC(NCM), HHD(NCM), HHE(NCM), HHF(NCM)
!
	IOPT=ITYP
	IMIN=0
	if (ITYP == 0) IMIN = 1
	if (ITYP == 0) IGUES = 0
	if (ITYP == 0) IOPT = 1
	IEXT(:) = 0
	IER1(:) = 0
	ICYC = 0
	EPS = 1.D-16
	PI = 3.1415926536D0
	RT = R*T
!-----------------------------------------------------------------------
!
!	Cantidades independientes del volumen
!
!	Las siguientes variables sirven a la subrutina VIRIAL:
	xNtot = sum(xN(:NC))
	x(:NC) = xN(:NC)/xNtot	!fraccion molar componente i
!
!	Contribución atractiva
!	Fraccion molar grupo k
	rNyT(:NG) = matmul(x(:NC),dfloat(ny(:NC,:NG)))
!	Área total grupos
	qT = dot_product(rNyT(:NG),q(:NG))
!	Fraccion area grupo k
	teta(:NG) = q(:NG)/qT*rNyT(:NG)
!
!	Contribución asociativa
	if (NST > 0) then
	  
	  sm_Xs(:2,:NST) = 0 !initializate values.

!	  Allocate internal (ZMAX only) variables which are passed here in COMMON vectors, but they are also input variables
!	  in the subroutine OptiNewton.
	  allocate (sm_aux(NST), Xs_aux(NST), Delta_aux(NST,NST), H_aux(NST,NST), pivot_vector(NST))
	  sm(:NST) = matmul(dfloat(sigma(:NST,:NC)),x(:NC))
	  sm_aux = sm(:NST)
	  
	  Delta_aux = Delta(:NST,:NST)

	  if (NST == 2) then
	    if (sm(2) > sm(1)) then
	  
	      major = 2	!this is only for the Analytical solution of 2 cross associating sites
	        
	    else
	
	      major = 1
	      
	    endif
	  endif
!
	endif
	Sum_S = sum(sm(:NST))
!
!     Contribucion repulsiva
	XLAM1 = 0.D0; XLAM2 = 0.D0; XLAM3 = 0.D0
	DO 3 I=1,NC

	  DIA = D(I)
	  DLA = X(I)*DIA
	  XLAM1 = XLAM1+DLA        !lambda 1
	  DLA = DLA*DIA
	  XLAM2 = XLAM2+DLA        !lambda 2
	  XLAM3 = XLAM3+DLA*DIA    !lambda 3

3	enddo
	PHI3 = XLAM3
!



!	write (*, *) B,Bvl,Bat,Bas
!  	Inicializacion de las fracciones no asociadas y deriv respecto de V
	Xs(:2,:NST) = 1D0
	dXs_dV(:2,:NST) = 0.D0
	nroIter = 0
!-------------------------------------------------------------------------------------
!
!	Inicialición de los volumenes de cada fase
	if (iGues >= 0 .AND. iGues <= 1) then

	  if (iOpt /= -1 .OR. iGues /= 1) then

!  	    Coeficiente virial: iniciación del volumen del gas, y algunas va-
!	    riables que necesita ZMAX.	  
	    call GCAvirial(NC, NG, NST, T, xn(:NC), x(:NC), B, Bvl, Bat, Bas, .TRUE.)	    
	    if (B > 0 .OR. (P < -RT/4D0/B)) then	!máximo en presión para B < 0

!	      Virial para la fase vapor: Z = 1 + B/v (más exacto que "1 + B*P/RT"... supuestamente)
	      Z(1) = 1D0 + 2D0*B*P/(R*T + dsqrt(RT**2 + 4D0*RT*B*P))
!
	    elseif (B*P/RT > -0.9) then
	
! 	      Try the common expression:
	      Z(1) = 1D0 + B*P/RT
		
	    else

!	      Gas ideal para la fase vapor	    
	      Z(1) = 1.D0	   
	    
	    endif

	  endif
	  expo = 0
6 	  if (iOpt /= 1 .OR. iGues /= 1) then

	    VTOT = 0.D0
	    DO 9 I = 1, NC

	      if (X(I) >= 1E-2) then

	        TR = T/TC(I)
	        TCRIT = .99D0*TC(I)
	        if (T > TCRIT) TR = 1.D0
	        
	        if (T < TCRIT) then
	          EXPO = (1.D0-TR)**.2857D0
	          FAC=.29D0**EXPO
	        else
	          FAC = 1.D0
	        endif
!	        Racket for the liq. molar volume of each compound
	        V = .29D0*TC(I)*R/PC(I)*FAC
	        VTOT = VTOT + V*X(I)

	      endif

9	    enddo
	    Z(2) = P*VTOT/RT

	  endif

	endif
!-------------------------------------------------------------------------------------
!
!	Volume calculation
!
100	ICYC = ICYC+1	 	! icyc=1 vapor,icyc=2 liquido

	IER = 0
	if (IOPT == 1 .AND. IGUES > 1) ICYC = 2
	NTRIAL = 0
	NHELP = 1
	NIT = 0
	ROLD = P/Z(ICYC)/RT
	RO = ROLD
	zeta3 = phi3*pi*ro/6
	if (zeta3 > 0.9999999) then
	  
!  	  The initial volume is lower than the CS-covolume
	  if (iCyc == 1) then
	    ro = P / (RT + P*pi*phi3/6)
	  elseif (iCyc == 2) then
	    ro = 4.2/phi3/pi
	  endif

	endif
	if (NST > 2 .OR. NST == 2 .AND. (Delta(1,1) /= 0 .OR. Delta(2,2) /= 0)) then

!	  Initialize non-bonded fraction with 3 direct substitution steps
	  do i = 1, 3
	    do j = 1, NST

	      Xs_aux(j) = 1D0/(1D0 + sum(ro*sm(:NST)*Xs(iCyc,:NST)*Delta(:NST,j)) )

	    enddo
	    Xs(iCyc,:NST) = Xs_aux(:NST)
	  enddo

	endif
!
!-----------------------------------------------------------------------
!	
!	Volume calculation loop
!	
200	CONTINUE

	NIT = NIT+1
300	CONTINUE
!  	Total iteration counter (durect subs. + Newton).
	nroIter(iCyc) = nroIter(iCyc) + 1
	DNOM = RO*QT/RT
	DFV = 0.D0; DFVDV = 0.D0
	dFVas(iCyc) = 0.D0; dFVVas = 0.D0
!
!	Attractive contribution
	DO 35 K=1,NG
		            !Según paper Jørgensen, estas son algunas de las variables:
	  HEL2K = 0.D0	!H2(k) = Sum(theta(j)*tau(j,k)*Qt*g(k,j)*rho/TR, j=1,NG)
	  HEL4K = 0.D0	!H4(k) = Sum(theta(j)*tau(j,k), j=1,NG)
	  HEL5K = 0.D0	!H5(k) = Sum(theta(j)*tau(j,k)*(g(j,k)*Qt*rho/RT)*alfa(j,k)*Dg(j,k)*Qt*rho/RT
	  HEL15K = 0.D0
	  HEL12K = 0.D0
	  HEL7K = 0.D0	!H7(k) = SUM(
	  AKK = A(K,K)
	  DO 5 J=1,NG

	    ARGU = ALFA(J,K)*(A(J,K) - AKK)*DNOM
!	    write (*, *) A(J,K),AKK,DNOM
	    AAH = A(J,K)*DNOM
	    AA(J,K)=AAH
	    E(ICYC,J,K)=DEXP(ARGU)
	    IF (TETA(J) >= 1.D-16) then

	      TETAE = TETA(J)*E(ICYC,J,K)
	      HEL4K = HEL4K + TETAE
	      TETAA = TETAE*AAH
	      HEL2K = HEL2K + TETAA
	      TETAA = TETAA*ARGU
	      HEL5K = HEL5K + TETAA
	      TETAA = TETAA*ARGU
	      HEL15K = HEL15K + TETAA
	      TETAE = TETAE*ARGU
	      HEL12K = HEL12K+TETAE
	      TETAE = TETAE*ARGU
	      HEL7K = HEL7K + TETAE
	      
	    endif

5	  enddo
	  HEL2K = HEL2K/HEL4K
	  HEL5K = HEL5K/HEL4K
	  HEL7K = HEL7K/HEL4K
	  HEL12K = HEL12K/HEL4K
	  HEL15K = HEL15K/HEL4K
	  IF (RNYT(K) >= 1.D-16) then

	    DFV = DFV + RNYT(K)*Q(K)*(HEL5K + HEL2K - HEL2K*HEL12K)
	    DFVDV = DFVDV - RO**2*RNYT(K)*Q(K)*((2.D0*HEL2K*HEL12K**2 - HEL2K*(HEL7K + 2.D0*HEL12K) - 2.D0*HEL5K*HEL12K) + HEL15K + 2.D0*HEL5K)*ZZ/2.D0
	    
	  endif
	  HELP2(ICYC,K) = HEL2K              !h2k/h4k
	  HELP4(ICYC,K) = HEL4K              !h4k
	  HELP5(ICYC,K) = HEL5K              !h5k/h4k
	  HELP7(ICYC,K) = HEL7K
	  HELP12(ICYC,K)= HEL12K            !h6k/h4k
	  HELP15(ICYC,K)= HEL15K
	  
35	enddo
!
!     DFV y DFVDV son la primera y segunda derivada de Ar/RT con respec al volumen
!
	DFV = RO*DFV*ZZ/2		! (dFdV)att
	DFVDV = DFVDV - RO*DFV*2		! (d2FdV2)att
!
!	Contribución repulsiva
	PI6 = PI/6*XLAM3
	Y = 1.D0/(1 - PI6*RO)
	DYDV = -Y**2*PI6*RO**2
	DYDVV = 2.D0/Y*DYDV**2 - 2*RO*DYDV
!
!     DFVR es la contribución repulsiva de DFV
	DFVR = (3*XLAM1*XLAM2/XLAM3 + XLAM2**3/XLAM3**2*(2*Y - 1 - 1.D0/Y) & 
	                                                               + 1.D0/Y)*DYDV
! 	
	PREP(ICYC) = -RT*DFVR	!P_rep. It is then send into a COMMON block

!	Acá se agrega la contribución respulsiva.
	DFV = DFV + DFVR	      ! (dFdV)att+rep
!
	DFVDV = DFVDV + DFVR*DYDVV/DYDV + (XLAM2**3/XLAM3**2 * (2.D0 + 1.D0/Y**2) - 1.D0/Y**2)*DYDV**2		 ! (d2FdV2)att+rep
!
!	Association contribution
	if (NST > 0) then
!
	  if (NST == 1) then
	    
!	    Analytical solution for 1 self-associating site (1A)
	    if (Delta(1,1) > 0) then
	    
	      s_Delta_V = sm(1)*Delta(1,1)*ro
	      root(iCyc) = dsqrt(4D0*s_Delta_V + 1)
	
	      Xs(iCyc,1) = 2D0/(1D0 + root(iCyc))
! 	      dXs_aux(iCyc,1) = -Xs(iCyc,1)/(1D0 + root)/root  !auxiliary variable for derivative calculations
	      dXs_dV(iCyc,1) = -2D0*Xs(iCyc,1)**2 /root(iCyc)*s_Delta_V*ro

	    else
		
	      s_Delta_V = 0
	      root(iCyc) = 1	      
		
	    endif
	    sm_Xs(iCyc,1) = sm(1)*Xs(iCyc,1)
	    
	  elseif (NST == 2 .AND. Delta(1,1) <= 0 .AND. Delta(2,2) <= 0) then
	    
!         Two cross associating sites, with same or different mole amounts.
	    if (Delta(1,2) > 0) then
	
	      s1_Delta_V = sm(3-major)*Delta(1,2)*ro
	      s2_Delta_V = sm(major)*Delta(1,2)*ro
	      b_aux(iCyc) = 1 + s1_Delta_V - s2_Delta_V
	      root(iCyc) = dsqrt( b_aux(iCyc)*b_aux(iCyc) + 4D0*s2_Delta_V )
	
	      Xs(iCyc,major) = 2D0/(b_aux(iCyc) + root(iCyc))

! 	      dXs_aux(iCyc,major) = -Xs(iCyc,major)/(b_aux + root)   !auxiliary variable for derivative calculations	      
	      dXs_dV(iCyc,major) = -Xs(iCyc,major)**2 * (1D0 - b_aux(iCyc) + (b_aux(iCyc)*(1D0 - b_aux(iCyc)) - 2D0*s2_Delta_V)/root(iCyc))*ro/2D0
	      if (sm(major) == sm(3-major)) then
	  
	        Xs(iCyc,3-major) = Xs(iCyc,major)
	        dXs_dV(iCyc,3-major) = dXs_dV(iCyc,major)
	  
	      else
	        
	        Xs(iCyc,3-major) = 1D0/(1 + s2_Delta_V*Xs(iCyc,major))
! 	        Xs(iCyc,3-major) = (1d0 - Xs(iCyc,major))/Xs(iCyc,major)/s1_Delta_V	
	        dXs_dV(iCyc,3-major) = Xs(iCyc,3-major)**2 * s2_Delta_V*(dXs_dV(iCyc,major) - Xs(iCyc,3-major)*ro)
	      
	      endif
	      sm_Xs(iCyc,major) = sm(major)*Xs(iCyc,major)		
	      sm_Xs(iCyc,3-major) = sm(3-major)*Xs(iCyc,3-major)		
	      
	    else
		
	      b_aux(iCyc) = 1
	      root(iCyc) = 1
	      sm_Xs(iCyc,:NST) = sm(:NST)
		
	    endif
	  
	  else
	    
! 	    General case for NST sites with any associating scheme:
	    
!	    Iniciación de las fracciones no asociadas por continuacion. La extrapolación
!	    está acotada ente 0 y 1:
	    do i = 1,NST

	      auxXs = Xs(iCyc,i) + dXs_dV(iCyc,i)*(1.D0/RO - 1.D0/ROld)
	      if (auxXs > 1D0) then
         
	        Xs(iCyc,i) = 1
         
	      elseif (auxXS <= 0.D0) then
         
	        Xs(iCyc,i)= Xs(iCyc,i)/5d0
         
	      else
         
	        Xs(iCyc,i) = auxXs
         
	      endif
         
	    enddo
!	    Calling general calculation of Xs
	    Xs_aux = Xs(iCyc,:NST)
	    call OptiNewton (maxIt, NST, sm_aux, Delta_aux, RO, Xs_aux, H_aux, pivot_vector, in)
	    Xs(iCyc,:NST) = Xs_aux
! 	  
!  Cálculo de la derivada de la fracción no asociada respecto VOLUMEN
!
!	    Number of moles of non-bonded sites
	    sm_Xs(iCyc,:NST) = sm(:NST)*Xs(iCyc,:NST)
	    
!	    Vector -dg/dv provisionally stored in dXs_dV:
	    dXs_dV(iCyc,:NST) = -sm(:NST)*matmul( Delta(:NST,:NST), sm_Xs(iCyc,:NST) )*ro*ro

! 	    do i = 1, NST	    
! 	      sm_Xs(iCyc,i) = 0.D0
! 	      dXs_dV(iCyc,i) = 0.D0
! 	      if (sm(i) >= 1D-16) then
! 	  	
! 	        sm_Xs(iCyc,i) = sm(i)*Xs(iCyc,i)
! 	        do j = 1, NST	  
! 	          if (sm(j) >= 1D-16 .AND. Delta(i,j) > 1D-16) then
! 	        
! 	            dXs_dV(iCyc,i) = dXs_dV(iCyc,i) - sm(i)*sm(j)*Xs(iCyc,j)*Delta(i,j)
! 	          
! 	          endif
! 	        enddo
! 	        dXs_dV(iCyc,i) = dXs_dV(iCyc,i)*ro*ro
! 	        
! 	      endif
! 	    enddo
	    
! 	    call LUDcmp (H(iCyc,:NST,:NST), NST, NST, indx(iCyc,:NST), d1)
	    call LUBksb (H_aux, NST, NST, pivot_vector, dXs_dV(iCyc,:NST))
!
	  endif

!	  Association contribution to d(Ar/RT)/dV y d2(Ar/RT)/dV2:
	  dFVas(iCyc) = ro*(Sum_S - sum(sm_Xs(iCyc,:NST)))/2
	  dFVVas = -(dFVas(iCyc) + dot_product(sm(:NST), dXs_dV(iCyc,:NST))/2)*ro
!
	endif
!
!	Agregado de la contribución asociativa
	dFV = dFV + dFVas(iCyc)
	dFVdV = dFVdV + dFVVas
!-------------------------------------------------------------------------------------
!
	DPDV(ICYC) = RT*(1.D0 + DFVDV/RO/RO)	! -V**2 * dP/dV
	G = RO*R*T - R*T*DFV - P	            ! g(rho) = Pcalc - Pspec
	DGRO = R*T*(1.D0 + DFVDV/RO**2)
!
!  Cálculo de (ro)j+1 por el metodo de Newton
!
	if (dgro > 0.D0 .OR. NIt <= 1) then
	
	  if (NHelp <= 1) then
	    if (iCyc == 2) then
	      if (iGues <= 2 .AND. iGues >= 0) then
! 	        Liquid
	        if (NTrial /= 0) then
	          
	          gCrit = g*gOld
	          
	          if (g > gOld .AND. dgro < 0.D0) goto 400
	            
	          if (gCrit < 0) then
	          
	            NHelp = NHelp + 1
!	            Density interpolation
! 	            ro = -g*(ro-rold)/(g-gold) + ro
!
!	            This was modified to avoid ro-rOld = 0
	            rOld2 = rOld
	            rOld  = ro
	            ro = -g*(ro - rOld2)/(g - gOld) + ro
	            
	            goto 300
	            
	          endif
	          
	        endif
	        gOld = g
	        rOld = ro
	        if (iGues == 0 .AND. g >= 0D0) then
	        
	          sro = Z(1) - P/RT/(ro - g/dgro) !diff between V and L
	          
	          if (sro < 0D0) goto 400
	          
	        endif
	        gc = dabs(g)
	        if (GC/dabs(dgro) < EPS) then
              
	          goto 400 !exit Newton loop if g = 0, since the solution has been found...
              
	        elseif (gc > 1) then
	    
! 	          For errors greater than certain value, it's better to approximate successively
	          rOld = ro
	          ro = ro - 0.04D0*ro*g/gc
	          NTrial = NTrial + 1
	          
	          goto 300
	        
	        endif	        
	        
	      endif
	    elseif (dgro <= 0.D0) then
	    
	      NTrial = NTrial + 1
	      if (NTrial > 25) then
	        
	        IER = 100
	        
	        goto 400
	        
	      endif
	      src = ro - g/dgro
	      if (src < 0) goto 400
	      ro = 0.9D0*ro
	      
	      goto 300
	      
	    endif
	  endif
	  NHelp = NHelp + 1
	  dro = -g/dgro
	  DltRo = ro/20
	  if (iCyc == 1) DltRo = 4*DltRo
	  if (dro > DltRo) dro = DltRo
	  if (dro < -DltRo) dro = -DltRo
	  
	else
!	  At present iteration, dP/drho < 0, which means mechanical unstability. Rho will be restored and
!	  calculation will proceed with a half Newton step.
	  rOld = ro
	  RO   = RO - DRO
	  DRO  = .5D0*DRO

	endif
	RCR = dabs(dro)
	if (RCR >= eps .AND. NIt <= 20) then

!
!	  Almacenamiento del valor anterior de la densidad...   y posterior actualizacion
	  rOld = ro	  
	  ro = ro + dro
	  goto 200
	
	elseif (NIt > 20) then
	
	  IEr = 100 + iCyc
	  
	endif
 !
!	Culminó el calculo del volumen de la fase iCyc
400 	if (NST > 0) then

!       Store variables from actual phase.
	  H(iCyc,:NST,:NST) = H_aux
	  indx(iCyc,:NST) = pivot_vector

	endif	
 	if (iGues <= 1) then
    
	  Pcrit = 1D3*dabs(g)
	  if (Pcrit > P) iExt(iCyc) = 1
	  if (dgro < 0D0 .OR. iExt(iCyc) == 1) IEr = 100*(icyc - 1) + 99
 	
 	endif
 	IEr1(iCyc) = IEr
 	Z(iCyc) = P/ro/RT
	IF (Igues <= 1) THEN
	
	  IF (iOpt == 1 .AND. icyc == 1) GOTO 100
	  
	  IF (iopt /= -1 .OR. iCyc /= 2) THEN
	  
	    IF (iOpt == -1 .AND. IEr /= 0 ) GOTO 100
          
          IF (iOpt == 1) THEN
          
	      IF ( ier.NE.0 ) THEN
	      
	        Z(2) = Z(1)
	        GOTO 450
	        
	      ELSE
		 
	        IF (iEr1(1) /= 0) GOTO 500
	        
	        Zcrt = DABS(Z(1) - Z(2))
	        IF (Zcrt <= 1.D-10) THEN
	        
	          IEr1(2) = 199
	          Z(2) = Z(1)
	          
	        ENDIF
	        
	        GOTO 450
	        
	      ENDIF
	      
	    ELSEIF (IEr == 0) THEN
          
	      GOTO 500
	      
	    ENDIF
	  ENDIF
	  IEr = 99
	  Z(1) = Z(2)
	  
	  GOTO 500 	

450 	  if (IMin /= 0 .AND. IEr1(2) == 0) then
!	   
!	    Cálculo de la energia residual de Gibbs para ambas fases.
!	      r     r
!	     G     A
!	    --- = --- + Z - 1
!	    R T   R T
!
	    GIB1=0.D0
	    GIB2=0.D0
!
!	    Contribucion atractiva: A_att
	    DO K=1,NG
         
	      GIB1 = GIB1 + RNYT(K)*Q(K)*HELP2(1,K)
	      GIB2 = GIB2 + RNYT(K)*Q(K)*HELP2(2,K)

	    enddo
	    GIB1=-GIB1*ZZ/2.D0       ! (A)att vapor
	    GIB2=-GIB2*ZZ/2.D0       ! (A)att liq
!
!	    Contribución asociativa: A_assoc
	    if (NST > 0) then
          
!	      Vapor:
	      Gib1 = Gib1 + dot_product(sm(:NST), dlog(Xs(1,:NST))) + (Sum_S - sum(sm_Xs(1,:NST)))*.5D0	!(A)att + (A)ass vapor
!	      Lí­quido:
	      Gib2 = Gib2 + dot_product(sm(:NST), dlog(Xs(2,:NST))) + (Sum_S - sum(sm_Xs(2,:NST)))*.5D0	!(A)att + (A)ass lí­quido
          
	    endif
!
!	    Contribución repulsiva
	    Y1 = 1.D0/(1.D0 - PI6*(P/Z(1)/RT))
	    R1 = XLAM1*XLAM2/XLAM3*3.D0
	    R2 = XLAM2**3/XLAM3**2
!         
	    GIB1 = GIB1 + R1*(Y1 - 1.D0) + R2*(Y1**2 - Y1 - DLOG(Y1)) + DLOG(Y1) + Z(1) - DLOG(Z(1))	         ! Gres/RT+1=Ares/RT+Z-lnZ vapor
	    GIB2 = GIB2 + R1*(Y  - 1.D0) + R2*(Y**2  - Y  - DLOG(Y))  + DLOG(Y)  + Z(2) - DLOG(Z(2))	         ! Gres/RT+1=Ares/RT+Z-lnZ lí­quido
!         
	    IC = 1
	    IF (GIB1 < GIB2) IC = -1
	    IPHAS = 1 + (IC + 1)/2
	
	    GOTO 201
	    
	  endif
	endif


! 	IF (DGRO > 0.D0 .OR. NIT <= 1) GOTO 26
! 
! !	  At present iteration, dP/drho < 0, which means mechanical unstability. Rho will be restored and
! !	  calculation will proceed with a half Newton step.
! 	  rOld = ro
! 	  RO   = RO-DRO
! 	  DRO  = .5D0*DRO
! 	  GOTO 25
! 
! 26	IF (NHELP > 1) GOTO 20
! 
! !	  Si la fase es lí­quida salta
! 	  IF (ICYC == 2) GOTO 21
! 
! 	  IF (DGRO > 0.D0) GOTO 20
! 
! 	    NTRIAL = NTRIAL + 1
! 	    if (NTRIAL  >  25) IER=100
! 
! 	  IF (IER /= 0) GOTO 10
! 
! 	    SRC=RO-G/DGRO
! 
! 	  if (SRC < 0.D0) GOTO 10
! 
! 	    RO=RO*.90D0
! 	    GOTO 24
! 
! 21	if (IGUES > 2 .OR. IGUES < 0) GOTO 20
!    
!         if (NTRIAL == 0) GOTO 22
!         
!         Gcrit=G*GOld
!         
!         if (G > GOLD .AND. DGRO < 0.D0) GOTO 10
!         
!         if (GCRIT < 0.D0) GOTO 23
! 
! 22	GOLD=G
! 	ROLD=RO
! 	if (IGUES /= 0 .OR. G < 0.D0) GOTO 27
! 
! 	  SRO = Z(1) - P/RT/(RO - G/DGRO)
! 
!       if (SRO < 0.D0) GOTO 10
! 27	GC = DABS(G)
! 	if (GC < EPS) then
! 
! 	  goto 200 !exit Newton loop if g = 0, since the solution has been found...
! 
! 	endif
! 	rOld = ro
!       RO = RO - (.04D0*G/GC*RO) !rho - 0.04g/|g|·rho = rho·(1 - 0.04 g/|g|) = si g > 0, rho disminuye un 4%; si g < 0, rho aumenta un 4%.
!       NTRIAL=NTRIAL+1
!       
!       GOTO 24 !arriba hasta el comienzo, sin sumar iteración.
! 
!    23 NHELP=NHELP+1
! 	rOld2 = rOld
! 	rOld  = ro
! !  Extrapolation of rho:
! !	RO    = -G*(RO - ROLD)/(G - GOLD) + RO
! !
! !  Modifiqué esto para  que ro - rOld sea <> 0
! 	RO    = -G*(RO - ROLD2)/(G - GOLD) + RO
! !  Creo que debería ser de esta manera:	
! ! 	RO    = -Gold*(RO - ROLD2)/(G - GOLD) + RO	
!       GOTO 24
! 
!       
! 20	NHELP = NHELP + 1
! !  Paso del método de Newton tradicional:
! 	DRO=-G/DGRO
! !  Paso máximo permitido = 5% del valor de rho para el líquido:      
! 	DLTRO=.05D0*RO
! !  20% del valor de rho para el gas:      
! 	if (ICYC < 2) DLTRO=4.D0*DLTRO
! 	if (DRO > DLTRO) DRO=DLTRO
! 	if (DRO < -DLTRO) DRO=-DLTRO
! 
! 25	RCR=DABS(DRO)
! 
! !  Criterio de salida:
! ! 	if (RCR < EPS .OR. NIT > 20) GOTO 200
! 	if (RCR < EPS .OR. NIT > 20) GOTO 200
!
!	Almacenamiento del valor anterior de la densidad...   y posterior actualizacion
! 	rOld = ro
! 	RO=RO+DRO

! 	GOTO 8
! !
! !	Culminó el calculo del volumen de la fase iCyc
! 200 	if (NST > 0) then
! 
! !       Store variables from actual phase.
! 	  H(iCyc,:NST,:NST) = H_aux
! 	  indx(iCyc,:NST) = pivot_vector
! 
! 	endif
! !200	if (NIT < 20) GOTO 10
! 	if (NIT == 20) then! GOTO 10
! 	  IER=100*ICYC
! 	endif
! 10	if (IGUES <= 1) then !GOTO 110
! 
! 	  PCRIT=1.D3*DABS(G)
! 	  if (PCRIT > P) IEXT(ICYC)=1
! 	  if (DGRO < 0.D0 .OR. IEXT(ICYC) == 1) IER = 100*(ICYC - 1) + 99
! 	  
! 	endif
! 110	IER1(ICYC)=IER
! 	Z(ICYC)=P/RO/R/T
! 	if (IGUES > 1) GOTO 113
! 	if (IOPT == 1 .AND. ICYC == 1) GOTO 4
! 	if (IOPT == -1 .AND. ICYC == 2) GOTO 115
! 	if (IOPT == -1 .AND. IER /= 0) GOTO 4
! 	if (IOPT == 1) GOTO 112
! 	IF (IER == 0) GOTO 113
! 115	IER=99
! 	Z(1)=Z(2)
! 	GOTO 113
! 112	if (IER /= 0) GOTO 114
! 	if (IER1(1) /= 0) GOTO 113
! 	ZCRT = DABS(Z(1) - Z(2))
! !	if (ZCRT > 1.D-10) GOTO 210
! 	if (Zcrt > 1.D-8) goto 210
! 	IER1(2)=199
! 114	Z(2)=Z(1)
! 
! 210	IF (IMIN == 0 .OR. IER1(2) /= 0) GOTO 113
! !-----------------------------------------------------------------------
! !	Cálculo de la energia residual de Gibbs para ambas fases.
! !	  r     r
! !	 G     A
! !	--- = --- + Z - 1
! !	R T   R T
! !
!       GIB1=0.D0
!       GIB2=0.D0
! !
! !	Contribucion atractiva: A_att
! 	DO 202 K=1,NG
! 
! 	  GIB1=GIB1+RNYT(K)*Q(K)*HELP2(1,K)
! 	  GIB2=GIB2+RNYT(K)*Q(K)*HELP2(2,K)
! 
! 202	enddo
! 
! 	GIB1=-GIB1*ZZ/2.D0       ! (A)att vapor
! 	GIB2=-GIB2*ZZ/2.D0       ! (A)att liq
! !
! !	Contribución asociativa: A_assoc
! 	if (NST > 0) then
! 
! !	  Vapor:
! 	  Gib1 = Gib1 + dot_product(sm(:NST), dlog(Xs(1,:NST))) + (Sum_S - sum(sm_Xs(1,:NST)))*.5D0	!(A)att + (A)ass vapor
! !	  Lí­quido:
! 	  Gib2 = Gib2 + dot_product(sm(:NST), dlog(Xs(2,:NST))) + (Sum_S - sum(sm_Xs(2,:NST)))*.5D0	!(A)att + (A)ass lí­quido
! 
! 	endif
! !
! !	Contribución repulsiva
! 	Y1 = 1.D0/(1.D0 - PI6*(P/Z(1)/RT))
! 	R1 = XLAM1*XLAM2/XLAM3*3.D0
! 	R2 = XLAM2**3/XLAM3**2
! !
! 	GIB1 = GIB1 + R1*(Y1 - 1.D0) + R2*(Y1**2 - Y1 - DLOG(Y1)) + DLOG(Y1) + Z(1) - DLOG(Z(1))	         ! Gres/RT+1=Ares/RT+Z-lnZ vapor
! 	GIB2 = GIB2 + R1*(Y  - 1.D0) + R2*(Y**2  - Y  - DLOG(Y))  + DLOG(Y)  + Z(2) - DLOG(Z(2))	         ! Gres/RT+1=Ares/RT+Z-lnZ lí­quido
! !
! 	IC=1
! 	IF (GIB1 < GIB2) IC=-1
! 	
! 	IPHAS = 1 + (IC + 1)/2
! 	GOTO 201
! 113	CONTINUE
500	IC = 1
	if (IER1(1) /= 0 .AND. IOPT == -1) IC = -1
	if (IER1(2) /= 0 .AND. IOPT == 1) IC = -1
201	CONTINUE
! 	if (NST > 0) then
! 
! 	  deallocate (Xs_aux, sm_aux, Delta_aux, H_aux, pivot_vector)
! 
! 	endif
!----------------------------
	RETURN
ENDsubroutine ZMAX
!---------------------------------------------------------------------------------------------
!
!	Esta subrutina obtiene los coeficientes de actividad de una mezcla
!	(preferentemente lí­quida) de NC componentes.C
!	Se basa en la fórmula
!
!	                 ------
!	ln gamma(i) = ln phi(i) - SUM( x(i)*ln phi(i), i=1,NC )
!
!	donde:
!	       ------
!	       phi(i) = coeficiente de fugacidad del comp. i en la mezcla
!
!	       x(i)   = fracción molar de i
!
!	       phi(i) = coeficiente de fugacidad de i puro a la misma T y P
!
subroutine GAMMA (NC, NG, NST, iTyp, iGz, iTemp, n, T, P, lnGamma, dlnGammadT, Z)
!
	implicit real*8 (A-H, O-Z)
!
	parameter (NCM = 30, NGM = 30, NSM=24)
	real*8 lnGamma(NCM), n(NCM), Nt
	dimension phi(NCM),dLPhi(NCM,NCM),dLPhiT(NCM),dLPhiP(NCM),dlnGammadT(NCM),x(NCM),aux(NCM)
!
	iDer = 0
!	write(3,*)ityp,igz
!
!	Normalización
	Nt = sum(n(:NC))
!	x(:NC) = n(:NC)/Nt
!
!	Coeficientes de fugacidad de la mezcla:
	call GCEOS (NC, NG, NST, iDer, iTemp, T, P, n, phi, dLPhi, dLPhiT, DLPHIP, Z, IGz, iTyp, IC)
!
!	Iniciación para cálculo:
	lnGamma(:NC) = phi(:NC)
	dlnGammadT(:NC) = dLPhiT(:NC)
!
!	Lazo de cálculo
	do i = 1, NC
!
!	  Cálculo de los coeficientes de fugacidad de los puros:
	  x(:NC) = 1.D-20
	  x(i) = 1D0
	  call GCEOS (NC, NG, NST, iDer, iTemp, T, P, x, phi, dLPhi, dLPhiT, DLPHIP, Zpure, IGz, iTyp, IC)
	  lnGamma(i) = lnGamma(i) - phi(i)
	  dlnGammadT(i) = dlnGammadT(i) - dLPhiT(i)
!
	enddo
!
!
	return
endsubroutine GAMMA
!
!---------------------------------------------------------------------------------------------
!
! 
! 
!	La subrutina VIRIAL halla el 2º coeficiente virial de una mezcla
!
!
subroutine GCAvirial (NC, NG, NST, T, n, x, B, Bvl, Bat, Bas, calc)

	implicit real*8 (a-h, o-z)
	parameter (NCM = 30, NGM = 30, NSM = 24, PI = 3.1415926536D0)
	real*8 n(NC),nt,lambda1,lambda2,lambda3
	dimension x(NC),aux(NCM)
	logical calc
	integer sigma
!
!	Este COMMON provee de variables importantes al resto de las subrutinas
	COMMON/ZZZ/HELP2(2,NGM),HELP4(2,NGM),HELP5(2,NGM),xg(NGM),QT,      &
     &           theta(NGM),HELP6(2,NGM),HELP11(2,NGM),HELP12(2,NGM),    &
     &           E(2,NGM,NGM),PREP(2),DPDV(2),lambda1,lambda2,lambda3
!
	COMMON/COORD/zz
!
!	Variables especí­ficas del término dispersivo
	COMMON/GROUP2/Q(NGM),g(NGM,NGM),DgDT(NGM,NGM),ALFA(NGM,NGM),R,NY(NCM,NGM)
!
!
!	Propiedades moleculares
	COMMON/MOL/DC(NCM),D(NCM),DT(NCM),HA(NCM),HB(NCM)
	COMMON/CRIT/TC(NCM),PC(NCM), omega(NCM)
!
!	Versión asociativa del common "ZZZ"
	common/ZZZAs/Xs(2,NSM),sm(NSM),dXs_dV(2,NSM)
!
!	Versión asociativa del common "GROUP2"
	common/GrupAs2/sigma(NSM,NCM), major
	common/GrupAs4/Delta(NSM,NSM),dDeldT(NSM,NSM)
!
!	¿Se calcularon variables que aprovecha ZMAX?
	if (.NOT. calc) then
!
	  nt = sum(n(:NC))
	  x(:NC) = n(:NC)/nt
	  xg(:NG) = matmul(x(:NC),dfloat(ny(:NC,:NG)))
	  QT = dot_product(q(:NG),xg(:NG))
	  theta(:NG) = q(:NG)/QT*xg(:NG)

	  lambda1 = dot_product(x(:NC),d(:NC))
	  lambda2 = dot_product(x(:NC),d(:NC)**2)
	  lambda3 = dot_product(x(:NC),d(:NC)**3)

	  if (NST > 0) sm(:NST) = matmul(dfloat(sigma(:NST,:NC)),x(:NC))

	endif
!
!	Contribución atractiva:
	aux(:NG) = matmul(g(:NG,:NG),theta(:NG))
	Bat = -zz/2D0*QT**2*dot_product(theta(:NG),aux(:NG))/R/T
!
!	Contribución repulsiva
	Bvl = PI*(lambda1*lambda2/2D0 + lambda3/6D0)
!
!	Contribución asociativa
	Bas = 0D0
	if (NST > 0 .AND. sum(sm(:NST)) >= 1D-16) then

	  aux(:NST) = matmul(sm(:NST),Delta(:NST,:NST))
	  Bas = -dot_product(aux(:NST),sm(:NST))/2D0/sum(sm(:NST))

	endif
!
	B = Bvl + Bat + Bas
!	write(*,'(<NG>("| ",<NG>(F12.0,2X)," |",/))')((g(i,j),j=1,NG),i=1,NG)
!	write (*, *) QT
!	write (*, *) theta(:NG)
!	write (*, *) B,Bvl,Bat,Bas
!
endsubroutine GCAvirial
!
!---------------------------------------------------------------------------------------------
!
!
!	OptiNewton halla el valor máximo de la función Q :
!
!	                          Qmax = A/RT|assoc,eq.
! 
!---------------------Francisco, 14/07/2011------durante v-1.8.0---------
!	Cálculo de la fraccion no asociada, como se muestra en el libro
!	"Thermodynamic models: Fundamentals & computational aspects" de Mollerup & Michelsen y
!	en M.Michelsen Ind.Eng.Chem.Res. 45(2006)8449-8453.
!
!	Siendo Lmáx = Aas/RT|eq:
!		cuando grad(Lmax,X) = 0, se obtiene la formula general de cál-
!		culo de las fracciones no asociadas.
!
!	Siendo que se trata de una optimizacion de L iterando sobre X, se re-
!	suelve por Newton-Raphson.
!           H (X(u + 1) - X(u)) + g = 0
!
!	Pasos:
!     calcular/estimar X
!     Calcula g
!     calcula H
!     resolver DX: X = X + DX
!   1 si L no aumenta,
!           X = X + alfa·DX    #alfa = 0.5
!           ir a 1
!     fin
!
!	El cálculo de las derivadas con respecto a V de la fracción no aso-
!	ciada Xs, están basadas en la modificación de A.E.Andreatta (ver Tesis
!	PhD y T.M.Soria y col. Fluid Phase Equilib. 302(2011)1-9) al traba-
!	jo de S.P.Tan y col. Ind.Eng.Chem.Res. 43(2004)203-208.
subroutine OptiNewton (maxit, NTS, sm, Delta, rho, X, H, indx, in)
!
	implicit real*8 (a-h, o-z)
!
!  Número máximo de iteraciónes
!  parameter(maxit = 1000)
!
!  Parámetros de entrada de asociación
	dimension delta(NTS,NTS),sm(NTS)
!
!  Parámetros de salida e intermedios
	dimension X(NTS), Xnew(NTS), g(NTS), H(NTS,NTS), minXindx(1), indx(NTS)
!
!  Parámetros auxiliares
	dimension Xold(NTS), DX(NTS), A(NTS,NTS), gnew(NTS), Hnew(NTS,NTS)
!
!  Cálculo inicial.
	call Qfunction(NTS,X,sm,rho,Delta,Q,g,H)
	im = 0
!
!  Comienzo del lazo de convergencia
	do in = 1,maxit
!
!  Valores de variables auxiliares
	  im = 1 + im
! 	  A = H
	  DX = -g
!
!  Resolviendo sistema de ecuaciones lineales
	  call LUDcmp (H, NTS, NTS, indx, d)
	  call LUBksb (H, NTS, NTS, indx, DX)
!	  Reiniciando alfa
	  alpha = 1.D0
!
!	  Calculando nuevo X
1	  do i = 1,NTS
	    Xnew(i) = X(i) + alpha*DX(i)
!
!	    Si una fracción es negativa de disminuye pero no se hace 0
	    if (Xnew(i) <= 0D0) then
	      Xnew(i) = 2.D-1*X(i)
	    endif
	  enddo
!
!	  Salida del lazo de convergencia si el error es menos que una tolerancia.
	  if (maxval(dabs(alpha*DX)) <= 1.d-15)exit
!
!	  Nuevo valor de Q:
2	  call Qfunction (NTS, Xnew, sm, rho, Delta, Qnew, gnew, Hnew)
!
!	  Aumentó Q?
	  if (Qnew > Q*(1D0 + 1.D-14)) then

3	    Q = Qnew
	    X = Xnew
	    g = gnew
	    
	    H = Hnew

	  else

!	    Si no lo hizo, disminuye a la 1/3 del paso.
	    alpha = alpha/3.
	    im = im + 1
	    goto 1

	  endif
	enddo
!
	return
	stop
100	format(5(F10.8,2X))
	endsubroutine
!---------------------------------------------------------------------------------------------
!
!
!	La subrutina Qfunction calcula el valor de la función Q, su gra-
!	diente y Hessiano para el cálculo de segundo orden de la fracción
!	no asociada, según lo sugerido por M. L. Michelsen, IECR 2002, 45,
!	8449-8453.
!
subroutine Qfunction (NTS, X, sm, rho, Delta, Q, g, H)
!
!	Q = Sum(sm(k)·(ln X(k) - X(k) + 1) - 1/2·Sum(Sum(sm(k)·sm(l)·Delta(k,l)·X(l)/V),l=1,NST),k=1,NST)
!
!	  = Q1 + Q2
!
!	g(k) = sm(k)/X(k) - sm(k) - Sum(sm(k)*sm(l)*Delta(k,l)*X(l)/V),l=1,NST)
!
!	  = sm(k)/X(k) - sm(k) - SUMA
!
!	"H"(k,l) = -(sm(k) + SUMA)/X(k)*d(k,l) - sm(k)*sm(l)*Delta(k,l)/V
!
!	Donde: sm(k)      = es la cantidad de moles de sitios "k"
!	       X(k)       = es la fracción de sitio "k" no asociada
!	       Delta(k,l  = fuerza de asociación entre sitios "k" y "l"
!	       NST        = número de sitios totales
!	       d(k,l)	= delta de Kronocker, = 1 si "l" = "k", sino = 0.
!	       g(k)       = gradiente de Q en la dirección "k"
!	      "H"(k,l)   = hessiano k,l
!	nota1: tomo a rho = 1/V donde V = vol total. Asume que son iguales porque los moles están normalizados (creo).
!
!	Francisco, febrero de 2010.
!
!	---------------------------------------
!	Siendo que cada vez que se corta un cálculo, el programa corta casi siempre en esta subrutina.
!	Eso indica que esta es un cuello de botella, por eso, comprimí­ el cálculo de las variables
!	Q, g y H en la menor cantidad de lazos posible. v-1.9.18
!
	implicit real*8 (a-h,o-z)
!	implicit real*16 (a-h,o-z)
	dimension X(NTS), sm(NTS), Delta(NTS,NTS), g(NTS), H(NTS,NTS)
!
!	Iniciación
	Q1 = 0
	Q2 = 0
	g = 0
	H = 0
!
!	Cálculo de Q1, gradiente y hessiano:
	do k = 1,NTS
	  if (sm(k) > 0) then
	    
	    do l = 1,NTS

	      if (Delta(l,k) > 0 .AND. sm(l) > 0) then
	        H(l,k) = -sm(k)*sm(l)*Delta(l,k)*rho
!	        Nota: g(k) aún NO ES el gradiente_k. Fran 18/02/2010.
	        g(k) = g(k) - H(l,k)*X(l)
	        Q2 = Q2 + H(l,k)*X(k)*X(l)/2
	      endif
	      
	    enddo
	    H(k,k) = H(k,k) - (sm(k) + g(k))/X(k)
!	    Para evitar elementos nulos en la diagonal.
	    if (H(k,k) == 0.0D0) H(k,k) = 1.0D-20
	    Q1 = Q1 + sm(k)*(dlog(X(k)) - X(k) + 1)
!
!	    Ahora sí­, g es el gradiente_k:
	    g(k) = sm(k)/X(k) - sm(k) - g(k)
	    
	  else
	      
! 	    Trace site, deactivated from calculation although they have their
!	    own X at infinite dilution:
	    H(k,k) = 1.D0
	    X(k) = 	1D0/(1 + sum(rho*sm(:NTS)*X(:NTS)*Delta(:NTS,k)) )
	      
	  endif

	enddo
!
!	Cálculo de Q:
	Q = Q1 + Q2
!
	return
endsubroutine
!--------------------------------------------------------------------------------
!
!
!	Halla la solución del sistema
!	                                 U · x = b
!	                                 =   -   -
!
!	por sustitición hacia atrás, complementando la subrutina LUDcmp que
!	descompone a A en U*L.
!	indx es el vector de pivoteo para b.
!
SUBROUTINE LUBksb (A, n, np, indx, b)
	INTEGER n, np, indx(n)
	REAL*8 a(np,np),b(n)
	INTEGER i,ii,j,ll
	REAL*8 sum
	ii=0
	do i = 1, n

	  ll=indx(i)
	  sum=b(ll)
	  b(ll)=b(i)
	  if (ii /= 0) then

	    do j=ii,i-1

	      sum=sum-a(i,j)*b(j)

	    enddo

	  else if (sum /= 0.D0) then

	    ii=i

	  endif
	  b(i)=sum

	enddo
	
	do i = n, 1 , -1

	  sum=b(i)
	  do j=i+1,n
	    sum=sum-a(i,j)*b(j)
	  enddo
	  b(i)=sum/a(i,i)

	enddo
	return
ENDsubroutine
!--------------------------------------------------------------------------------
!
!	LUDcmp halla la matriz U tal que
!
!	                         A = U · L
!	                         =   =   =
!	Donde U es una matriz "diagonal superior" y L "diagonal inferior
!	con sus elementos diagonales unitarios".
!
SUBROUTINE LUDcmp(a, n, np, indx, d)
	INTEGER n, np, indx(n), NMAX
	REAL*8 d,a(np,np),TINY
	PARAMETER (NMAX = 500, TINY=1.0e-20)
	INTEGER i,imax,j,k
	REAL*8 aamax, dum, sum, vv(NMAX)
      
	d = 1.D0
	do i=1,n

	  aamax=0.D0
	  do j=1,n

	    if (abs(a(i,j)) > aamax) aamax = abs(a(i,j))

	  enddo
	  if (aamax == 0.D0) then

	    write(*,'(/, "Matriz singular en LUDcmp", /)')
	    !stop

	  endif
        vv(i) = 1.D0/aamax
	enddo
      do j = 1, n

	  do i = 1, j - 1

	    sum=a(i,j)
	    do k=1,i-1

	      sum=sum-a(i,k)*a(k,j)

	    enddo
	    a(i,j)=sum

	  enddo
	  aamax = 0.D0
	  do i = j, n

	    sum = a(i,j)
	    do k=1,j-1

	      sum = sum - a(i,k)*a(k,j)

	    enddo
	    a(i,j)=sum
	    dum=vv(i)*abs(sum)
	    if (dum >= aamax) then

	      imax=i
	      aamax=dum

	    endif

	  enddo
	  if (j /= imax) then

	    do k=1,n

	      dum=a(imax,k)
	      a(imax,k)=a(j,k)
	      a(j,k)=dum
	      
	    enddo
	    d=-d
	    vv(imax)=vv(j) 

	  endif
	  indx(j)=imax
	  if (a(j,j) == 0.D0)a(j,j)=TINY
	  if (j /= n) then

	    dum=1.D0/a(j,j)
	    do i=j+1,n

	      a(i,j)=a(i,j)*dum

	    enddo
        
	  endif

	enddo
!
      return
ENDsubroutine LUDcmp
!
!--------------------------------------------------------------------------------

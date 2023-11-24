subroutine SRKPR_version (CEq, version, thermo_name)

	integer, intent(in)        :: CEq
	character(20), intent(out) :: thermo_name, version

	version = "0.3.5.0"
	if (CEq == 0) then
	  thermo_name = "SRK-EOS"
	else
	  thermo_name = "PR-EOS"
	endif
	return
endsubroutine SRKPR_version
!---------------------------------------------------------------------------------------
!
!
!
!
!------------------------------------------------------------------------------------------
!
!  Subroutine ANEW calculates parameters a and b and their partial composition 
!  and temperature derivatives.
!
SUBROUTINE ANEW(X, A0, AT, ATT, B0, AD1, AD2, ADT, BD1)
!  X:      (I):      Normalized composition
!  A0:     (O):      Mixture "a" (a/RT)-PARAMETER
!  AT:     (O):      T-derivative of A0
!  ATT:    (O):      T-derivative of AT
!  B0:     (O):      Mixture b-parameter
!  AD1:    (O):      Compsotion derivative of A0
!  ADT:    (O):      T-derivative of AD1
!
	IMPLICIT REAL*8(A-H,O-Z)	
	integer, parameter ::        NCM = 30, NSM = 24
	
	integer ::                   i, j, NADR, NDER, nTemp, sigma
	
	real(8) ::                   FAC, P, S0, S1, T 
	real(8), dimension(NCM) ::   Ac0, Ac1, AD1, ADT, X, Ax0, Ax1, As1, AsT1, Bc, BD1
	
	real(8), dimension(NCM*(NCM+1)/2) ::   AD2, CK
	
	COMMON/STYR/                 N, NH, NDER, NTEMP, NST    
	COMMON/PAR/                  CK, BC, AC0, AC1, b_matrix(NCM,NCM)
	common /GrupAs2/             sigma(NSM,NCM), major	
	common /GrupAs4/             Delta_T(NSM,NSM), dDeldT_T(NSM,NSM) !Temperature dependent-only part of Delta	
	
	S0 = 0.D0
	S1 = 0.D0
	DO 10 I = 1,N
	
	  XI = X(I)
	  AS1(I) = 0.D0
	  AX0(I) = XI*AC0(I)		  
	  IF (NTEMP > 0)then
	  
	    AST1(I) = 0.D0			  
	    AX1(I) = XI*AC1(I)		  
	    S1 = S1 + AX1(I)			
	  endif
	  S0 = S0 + AX0(i)  
	  
10	enddo
	A0 = S0*S0	    
	NADR = 0	
15	DO 20 I = 1, N
	
	  FAC = 2*AC0(I)
	  AD1(I) = S0*FAC
	  
	  IF (NTEMP .NE. 0) ADT(I) = FAC*S1 + 2*S0*AC1(I)				 
	  
	  IF (NDER .GE. 2)then	!LT.2) GO TO 20	 
	    DO 25 J=I,N
	    
	      NADR = NADR+1			 
	      AD2(NADR) = FAC*AC0(J)	    
	
25	    enddo
	  endif
	  
20	Enddo
	NADR = 0
!	
!  START MAIN LOOP FOR A-PARAMETER
	DO 50 I = 1,N
	
	  IF(I < N)then	!EQ.N) GO TO 50	     
	  
	    K1 = I + 1				
	    NADR = NADR + 1			 
	    XI = X(I)			     
	    AX0I = AX0(I)			 
	    
	    IF (NTEMP .NE. 0) AX1I = AX1(I)  
	    
	    FAC = 2*AC0(I)			
	    DO 60 K = K1, N			
!	
!  START INNER LOOP		  
!
	      NADR = NADR + 1			 
	      EF = CK(NADR)			 
	      IF(EF /= 0)then
	      
	        AX0K = AX0(K)			 
	        AS1(I) = AS1(I) + AX0K*EF	   
	        AS1(K) = AS1(K) + AX0I*EF	   
	        IF (NTEMP .GT. 0)then
	        
	          AST1(I) = AST1(I) + AX1(K)*EF     
	          AST1(K) = AST1(K) + AX1I*EF	 
	          
	        endif
65	        IF (NDER .GE. 2)then
	          
	          FAC1 = AC0(K)*FAC*EF		
	          AD2(NADR) = AD2(NADR) - FAC1	
	          
	        endif
	      endif
60	    enddo

	  endif

50	enddo		    
	A0 = 0D0
	AT = 0.D0			 
	DO 70 I = 1,N			 
	  
	  FAC = 2*AC0(I)			
	  AD1(I) = AD1(I) - FAC*AS1(I)	
	  A0 = A0 + AD1(I)*X(I)/2D0	     
	  IF (NTEMP.GT.0)then
	  
	    ADT(I) = ADT(I) - 2*AC1(I)*AS1(I) - FAC*AST1(I)				 
	    AT = AT + X(I)*ADT(I)		 
	    
	  endif
70	enddo			    
	AT = AT/2D0
!					     
!     START B-CALCULATION	     
!					     
	B0 = 0.D0				 
	DO 120 I = 1,N			
	
	  B0 = B0 + X(I)*BC(I)		  
	  BD1(I) = BC(I)			
	  
120	enddo

! 	Complete the matrix of covolumes, in order of evaluate the association strenght. This may be also useful
! 	if an l_ij is introduced.
	do k = 1, N
	  do l = k, N
	    
	    b_matrix(l,k) = (bc(l) + bc(k))/2
	    b_matrix(k,l) = b_matrix(l,k)

      enddo
    enddo
    do k = 1, N
	  do i = 1, NST
	    if (sigma(i,k) > 0) then
		
            do l = k, N  
	        do j = i, NST	     		  
	          if (sigma(j,l) > 0) then
	    
	            if(Delta_T(i,j) > 0) then
	
! 	              kappa(j,i) = beta(j,i)*b_matrix(l,k)
! 	              kappa(i,j) =  kappa(j,i)
	              Delta_T(i,j) = Delta_T(i,j)*b_matrix(l,k)
	              Delta_T(j,i) = Delta_T(i,j)
	              if (nTemp == 1) then
          
	                dDeldT_T(i,j) = dDeldT_T(i,j)*b_matrix(l,k)
	                dDeldT_T(j,i) = dDeldT_T(i,j)
                  
	              endif
	            
	            endif

	          endif	          
	        enddo
	      enddo  
	      
	    endif		    
	  enddo
	enddo

	RETURN				
ENDsubroutine ANEW
!-----------------------------------------------------------------------------
!
!  CUBIC solves for the compressibility factor			     
!  MTYP=1: liquid; MTYP=-1: vapour; MTYP= 0: min. Gibbs energy
!
!
!     MTYP:   (I):      DESIRED PHASE TYPE 1=LIQ,-1=VAP,0=MIN. G
!     A:      (I):      EOS AP/T
!     B:      (I):      EOS BP/T
!     Z:      (O):      RETURNED COMPRESSIBILITY FACTOR
!				     
SUBROUTINE CUBIC (MTYP, A, B, ZZ, IC)  

	IMPLICIT REAL*8 (A-H, O-Z)
	
	integer, parameter    :: NSM = 24
	                      
	integer               :: iCyc
	real(8), dimension(2) :: Aassoc, Z
	
	COMMON /CUB/             C, C1, C2
	COMMON /STYR/            NC, NH, NDER, NTEMP, NST	
	common /ZZZAs/           Xs(2,NSM), sm(NSM), dXs_dV(2,NSM)
	common /GrupAs3/         sm_Xs(2,NSM), Sum_S, dFVas(2)
	common /phase_ID/        iCyc
	
	IC = 0
	iCyc = 1
	BC = B*C
	X2 = 1D0 - BC
	X1 = A - B*(1D0 + B + C + 2*BC)	     
	X0 = B*(A - BC*(1D0 + B))
	
	if (NST == 0) then      
	  
	  Zc = X2/3D0                          !Zc (1/3 for SRK, 0.31 for PR...)
	
!  Evaluate F at Z = Zc
	  F = ((Zc - X2)*Zc + X1)*Zc - X0      !Z - (Z/(Z - B) -A Z/(Z (Z + B) + C*B*(Z - B))	
	    

	  Z(1) = B                                !Liquid root
	  IF(F < 0. .OR. B > Zc) Z(1) = Z(1) + 1D0   !Vapor root of F < 0 (volume was too small) or if it may be supercritical.
!
!  If B is small, check for the possibility of a zero pressure solution
	  IF ( B < 1.D-5 ) THEN
	    DD = X1*X1 - 4.D0*X2*X0
	    IF ( DD > 0.D0 ) Z(1) = B
	  ENDIF
	  
	else
	  
	  B2 = 0.475D0*B
	  if (mTyp <= 0) then
	    Z(1) = 1
	  else
	    Z(1) = 1.1D0*B
	  endif
	  
	endif
      
!  Third order newton for cubic EoS
	Z0 = Z(1)
	call third_Order_Newton (B, B2, BC, X0, X1, X2, Z0, Z(1), DF1, DF2)
! 	DZ = 1
! 10	do while (DZ > 1.D-10)
! 
! 	  DF2 = 3*Z-X2
! 	  DF1 = Z*(DF2 - X2) + X1	  
! 	  F = ((Z - X2)*Z + X1)*Z - X0
! 	  DZ = F/DF1			    
! 	  DZ = DZ*(1D0 + DZ*DF2/DF1)
!           
! 	  Z = Z - DZ	
!       
! 	enddo
	
!  Converged; how many roots are desired?

!   
!     DF2 = d2F/dZ2. If DF2 < 0, we have found a liquid root. If else, we have a vapor root.
!	if (DF2 >= 0) then
!	  
!	  IC = -1
!	  
!	else
!	  
!	  IC = 1
!	  
!	endif
	IC = 1    
	IF (MTYP*DF2 >= 0) then

	  if (NST == 0) then	
	
! 
!  If the product mTyp*DF2 is greater than 0, it means we have found the "other" phase, not the
!  desired.
!  If mTyp = 0, we WANT TO check both phases.
	    E1 = Z(1) - X2			     
	    E0 = Z(1)*E1 + X1
	  
!  D < 0 means no more roots (supercritical?)  
	    D = E1*E1 - 4.D0*E0			
	    IF(D >= 0.) then 

!  Get remaining
	      Z(2) = (DABS(E1) + DSQRT(D))/2	
	      IF (E1 > 0.) Z(2) = -Z(2)
	      IF (Z(1) > Z(2)) Z(2) = E0/Z(2)
! 	    
!  Refine by a single newton step
	      DF2 = 3.D0*Z(2) - X2
	      DF1 = Z(2)*(DF2 - X2) + X1
	      F = ((Z(2) - X2)*Z(2) + X1)*Z(2) - X0
	      DF1R = 1.D0/DF1
	      DZ = F*DF1R
	      DZ = DZ*(1D0 + DZ*DF2*DF1R)
	      Z(2) = Z(2) - DZ
            
	    else

	      IC = -1
        
        endif
	  else
	    if (DF2 > 0 .AND. mTyp == 1) then
	  
	      D = -1 !asume supercritical
		
	    elseif (DF2 < 0 .AND. mTyp == 1) then
	  
	      D = -1 !asume supercritical
	      
	    else
	      if (mTyp == 0 .AND. DF2 > 0) then
	  
	        Z0 = 1.1D0*B ! now we need the liquid
            IC = -1
	      
	      elseif (mTyp == 0 .AND. DF2 < 0) then
	  
	        Z0 = 2 !we already found a liquid when we started with a gas (mType = 0)
	      
	      endif
	      iCyc = 2	    
	      call third_Order_Newton (B, B2, BC, X0, X1, X2, Z0, Z(iCyc), DF1, DF2)
	      if (dabs(Z(1) - Z(2)) < 1D-8) then
	  
	        D = -1
	        iCyc = 1
	        
	      else
	    
	        D = 1
	        
	      endif
	    endif
	  endif
	  if (D > 0) then
	    IF (Z(2) > B) then
	      IF (MTYP == 0) then

!	        Check Gibbs energy of both phases:
	        F = DLOG((Z(1) - B)/(Z(2) - B)) + A/B/(C2 - C1)*DLOG((Z(1) + C2*B)*(Z(2) + C1*B)/(Z(1) + C1*B)/(Z(2) + C2*B)) - (Z(1) - Z(2))
	        if (NST > 0) then
	  
	          Aassoc(1) = sum(sm(:NST)*dlog(Xs(1,:NST)) - sm_Xs(1,:NST)/2) + Sum_S/2
	          Aassoc(2) = sum(sm(:NST)*dlog(Xs(2,:NST)) - sm_Xs(2,:NST)/2) + Sum_S/2
	          F = F + Aassoc(2) - Aassoc(1)
	  
	        endif
	        IF (F < 0.D0) then 
		      iCyc = 2
	          if (Z(2) > Z(1)) then
                IC = -1
	          else
	            IC = 1
	          endif
              Z(1) = Z(2)
!	          IC = -IC
	          
	        elseif (F > 0) then
	          iCyc = 1
	          if (Z(1) > Z(2)) then
                IC = -1
	          else
	            IC = 1
	          endif

	        endif
	        
	      elseIF (MTYP == 1) then
	    
	        IF (Z(2) < Z(1)) then
	          Z(1) = Z(2)
	        else
	          iCyc = 1
	        endif
	        IC = MTYP	        
	    
	      else
	      
60	        IF (Z(2) > Z(1)) then
	          Z(1) = Z(2)
	        else
	          iCyc = 1
	        endif       
	        IC = MTYP
	  
	      endif

	    endif
	  endif
 
	endif
	  
	ZZ = Z(1)
	RETURN	
ENDsubroutine CUBIC
!----------------------------------------------------------------------------------------------
real(8) function function_Zcubic(X0, X1, X2, Z) result (F)

	implicit none
	real(8) :: X0, X1, X2, Z
	
	F = ((Z - X2)*Z + X1)*Z - X0
	
endfunction function_Zcubic
	
!----------------------------------------------------------------------------------------------
!
! THIRD_ORDER_NEWTON calculates the root of F(Z). This corresponds to the while loop
! in subrroutine CUBIC in CEOS version prior to 0.1.0. It was moved outside when adding association
!
! c = 0 (SRK)
! F(Z) = Z*(Z - B)*(Z + B) - Z*(Z + B) + A*(Z - B) - Zassoc*(Z - B)*(Z + B)
!
! c = 1 (PR)
! F(Z) = (Z - B)*(Z*(Z + B) + c*B(Z - B)) - Z*(Z + B) - c*B(Z - B) + A*(Z - B) - Zassoc*(Z - B)*[(Z + B) + c*(Z - B)*B/Z]
!
SUBROUTINE third_Order_Newton (B, B2, BC, X0, X1, X2, Z0, Z, DF1, DF2)
  
	implicit none
	
	integer, parameter                 :: max_iter = 20, n_screen = 10
	real(8), parameter                 :: Zmax = 2
	                                   
	integer                            :: i, iter, j, n_candidates = 0
	                                   
	real(8), intent(in)                :: B, B2, BC, X0, X1, X2, Z0
	real(8)                            :: alpha, aux, aux2, D_D1, DZ = 1D0, F, function_Zcubic, V_dLng_dV, & 
	                                      Z, Znew, DF1, DF2, Zassoc, dZassoc1, dZassoc2
	                       
	integer, dimension(:), allocatable :: candidate
	real(8), dimension(:), allocatable :: Zscreen, Fscreen
	
! 	variables from COMMON
	integer                            :: NC, NH, NDER, NTEMP, NST 
	COMMON /STYR/                         NC, NH, NDER, NTEMP, NST 

	Z = Z0
	iter = 0
!	Main loop. 
	
	do while (dabs(DZ) > 1.D-10)
	  
	  iter = iter + 1	  
	  if (iter == max_iter) then
	  
! 	    Lets make a screening of F varying Z. If some possible roots are found,  
!	    some iterations of Bibisection method will be performed. If more than
!	    one possible root if found, we should prefer the closer to Z0.

! 	    Start screening with logarithmic scaling, within the interval (B,Zmax).
	    allocate (Zscreen(n_screen+1), Fscreen(n_screen+1), candidate(n_screen+1))
	    n_candidates = 0
	    candidate(:n_screen+1) = 0
	    Zscreen(1) = B
	    aux = dlog(B)
	    Fscreen(1) = function_Zcubic(X0, X1, X2, Zscreen(1))
	    if (NST > 0) then
	  
!	      Association contribution    
	      call Z_assoc (B2, Zscreen(1), Zassoc, dZassoc1, dZassoc2)
	      Fscreen(1) = Fscreen(1) - Zassoc*(Zscreen(1) - B)*(Zscreen(1) + B + BC*(Zscreen(1) - B)/Zscreen(1))
	
	    endif	    
	    aux2 = dlog(Zmax)
	    D_D1 = (aux2 - aux)/dble(n_screen)
	    do i = 1, n_screen
	
	      Zscreen(i+1) = dexp(aux + D_D1*i)
	      Fscreen(i+1) = function_Zcubic(X0, X1, X2, Zscreen(i+1))
	      if (NST > 0) then
	    
!	        Association contribution    
	        call Z_assoc (B2, Zscreen(i+1), Zassoc, dZassoc1, dZassoc2)
	        Fscreen(i+1) = Fscreen(i+1) - Zassoc*(Zscreen(i+1) - B)*(Zscreen(i+1) + B + BC*(Zscreen(i+1) - B)/Zscreen(i+1))
	  
	      endif
! 	      Look for a candidate:
	      if (Fscreen(i)*Fscreen(i+1) < 0D0) then
	  
	        n_candidates = n_candidates + 1
	        candidate(n_candidates) = i
	 
	        if (n_candidates > 1) then
	          if (dabs(Z0 - Zscreen(i+1)) > dabs(Z0 - Zscreen(candidate(n_candidates-1))) ) then
	            n_candidates = n_candidates - 1
	            exit !we've started to move away from Z0
	          endif
	        endif
	        
	      endif
	      
	    enddo
	    if (n_candidates == 0) then
	
! 	      No candidates have been found. 
	      Z = Zmax !this one shouldn't fail...
	      
	    else
	
! 	      Screening finished, start Bibisection within the interval with candidate closer to Z0
	      i = candidate(n_candidates)
	      do j = 1, 3
	    
	        Z = (Zscreen(i) + Zscreen(i+1))/2
	        F = function_Zcubic(X0, X1, X2, Z)
	        if (NST > 0) then
	    
!               Association contribution    
	          call Z_assoc (B2, Z, Zassoc, dZassoc1, dZassoc2)
	          F = F - Zassoc*(Z - B)*(Z + B + BC*(Z - B)/Z)
		    
	        endif
	        if (F*Fscreen(i+1) < 0) then	    
	          Zscreen(i) = Z
	          Fscreen(i) = F
	        else   
	          Zscreen(i+1) = Z
	          Fscreen(i+1) = F	          
	        endif	  
		  
	      enddo
	    
	    endif
	    deallocate (Zscreen, Fscreen)	    
	  endif
! 	  Z_old = Z
	  DF2 = 3*Z-X2	! 1/2 * d2F/dZ2
	  DF1 = Z*(DF2 - X2) + X1	  
	  F = ((Z - X2)*Z + X1)*Z - X0
	  if (NST > 0) then
	    
!         Association contribution    
	    call Z_assoc (B2, Z, Zassoc, dZassoc1, dZassoc2)
	    aux2 = (Z + B + BC*(Z - B)/Z)
	    aux = (Z - B)*aux2
	    F = F - Zassoc*aux
! 	    aux2 = 2*Z - B**2*BC/Z/Z +BC
	    DF1 = DF1 - dZassoc1*aux - Zassoc*(2*Z - B**2*BC/Z/Z +BC)
	    DF2 = DF2 + (-dZassoc2*aux - 2*dZassoc1*(aux2+(Z-B)*(1D0+BC/Z-BC*(Z-B)/Z**2))    &
	              - Zassoc*((2*(1D0+(BC)/(Z)-(BC*(Z-B))/(Z*Z))+(Z-B)* (-(2*BC)/(Z*Z)+(2*BC* (Z-B))/(Z*Z*Z)))))/2
	  
	  endif
	  
	  DZ = F/DF1			    
! 	  DZ = DZ*(1D0 + DZ*DF2/DF1)	  
	  
!  See pages 93 to 94 of "Thermodynamic Models: Fundamentals & 
!  Computational Aspects", M.L. Michelsen, J.M. Mollerup, 2nd Ed.	 
! 
!  (In page 94, Eq. (134) and (135) are probably wrong, there 
!  should be Delta and not Delta1).
!
	  if (DZ*DF2/DF1 < -5D-1) then
	    
	    D_D1 = -5D-1  
	    
	  elseif (dabs(DZ*DF2/DF1) > 1D0) then
	    
	    D_D1 = 0D0
	    
	  else
	    
	    D_D1 = DZ*DF2/DF1
	    
	  endif
	  DZ = DZ*(1D0 + D_D1)
   
!  To avoid negative values of the volume:	  
	  alpha = 1D0
	  do
	    
	    Znew = Z - DZ*alpha
	    if (Znew < B) then
		
	      alpha = alpha/3.1622776602D0
	      
	    else
		
	      exit
	      
	    endif
	  enddo
	  Z = Znew	
      
	enddo
	
	DZ = 1D0
	return
endsubroutine third_Order_Newton
!----------------------------------------------------------------------------------------------
!
!  Evaluate the association contribution to the compressibility factor and its derivatives with
!  respect the... compressibility factor.
! 
subroutine Z_assoc (B2, Z, ZassocOut, dZassoc1, dZassoc2)
	
	implicit none
	
	integer, parameter   :: maxIt = 100, NCM = 30, NSM = 24
	real(8), parameter   :: R = 82.05D0, tol = 1D-14
	
	real(8), intent(in)  :: B2, Z
	real(8)              :: ZassocOut, dZassoc1, dZassoc2
! 	Auxiliary variables
	integer              :: i, j, k, l, N_iter_Xs = 0
	real(8)              :: aux = 0D0, rho, s1_dDV_V, s1_Delta_V, s2_Delta_V, sum2, V_dLng_dV,    &
	                        d2Xs_dV2(2,NSM)
	real(8), save        :: rhOld(2), XsOld(2,NSM), dXs_dVOld(2,NSM)
	
	integer, dimension(:),   allocatable :: pivot_vector
	real(8), dimension(:,:), allocatable :: Delta_aux, H_aux
	real(8), dimension(:),   allocatable :: Xs_aux, sm_aux	                        
	
! 	Variables in COMMON vector:
	integer              :: iCyc, indx, major, NC, NH, NDER, NTEMP, NST, sigma
	real(8)              :: AC0, AC1, b_aux, b_matrix, BC, beta, CK, d2DeldV2, dDeldT, dDeldT_T, dDeldV,&
	                        Delta, Delta_T, dFVas, dXs_dV, d2Zassoc_dZ2, dZassoc_dZ, eps_R, H, kappa, P,&
	                        root, sm, sm_Xs, Sum_S, sum1, T, x, Xs, Zassoc
	                        
!
	COMMON /STYR/           NC, NH, NDER, NTEMP, NST
	common /phase_ID/       iCyc	
	COMMON /PAR/            CK(NCM*(NCM+1)/2), BC(NCM), AC0(NCM), AC1(NCM), b_matrix(NCM,NCM)
!
	common /ZZZAs/          Xs(2,NSM), sm(NSM), dXs_dV(2,NSM)
	common /ZZZAs2/         Zassoc(2), dZassoc_dZ(2), d2Zassoc_dZ2(2), sum1(2)
!
	common /GrupAs1/        kappa(NSM,NSM), eps_R(NSM,NSM)
	common /GrupAs1bis/     beta(NSM,NSM)
	common /GrupAs2/        sigma(NSM,NCM), major
	common /GrupAs4/        Delta_T(NSM,NSM), dDeldT_T(NSM,NSM) !Temperature dependent-only part of Delta
	common /GrupAs4bis/     Delta(NSM,NSM), dDeldT(NSM,NSM), dDeldV(NSM,NSM), d2DeldV2(NSM,NSM)

	common /GrupAs3/        sm_Xs(2,NSM), Sum_S, dFVas(2)
	common /GrupAs5/        H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)
	
	common /PTx/            P, T, x(NCM)

	sm_Xs(iCyc,:NST) = 0 !initializate values.

!	Allocate internal variables which are passed here in COMMON vectors, but they are also input variables
!	in the subroutine OptiNewton.
	allocate (sm_aux(NST), Xs_aux(NST), Delta_aux(NST,NST), H_aux(NST,NST), pivot_vector(NST))
	sm(:NST) = matmul(dfloat(sigma(:NST,:NC)),x(:NC))
	sm_aux = sm(:NST)
	  
	
	rho = P/Z/T !rho*R [atm/K], same units as b_matrix
	v_dLng_dv = -B2/(Z - B2)	
	sm(:NST) = matmul(dfloat(sigma(:NST,:NC)), x(:NC))
	Sum_S = sum(sm(:NST))
	
! 	Evaluate complete Delta matrix
	do i = 1, NST 
	  do j = i, NST	          
	 
	    if(Delta_T(j,i) > 0) then	  
	      Delta(j,i) = Delta_T(j,i)*Z/(Z - B2)
	      Delta(i,j) = Delta(j,i)
	      dDeldV(j,i) = Delta(j,i)*v_dLng_dv/Z*P/T
	      dDeldV(i,j) = dDeldV(j,i)
	      d2DeldV2(j,i) = Delta(j,i)*((v_dLng_dv/Z)**2 + (2*Z - B2)*B2/(Z - B2)**2 /Z**2 )*P*P/T/T
	      d2DeldV2(i,j) = d2DeldV2(j,i)	
	    endif
	        
	  enddo
	enddo
	Delta_aux = Delta(:NST,:NST)

	if (NST == 2) then
	  if (sm(2) > sm(1)) then
	
	    major = 2	!this is only for the analytical solution of 2 cross associating sites
	
	  else
    
	    major = 1
	    
	  endif
	endif
	Sum_S = sum(sm(:NST))
!
	do i = 1, NST
	  indx(iCyc,i) = i
	enddo
! 	if (NST == 1) then
! 	  
! !	  Analytical solution for 1 self-associating site (1A)
! 	  if (Delta(1,1) > 0) then
! 	  
! 	    s1_Delta_V = sm(1)*Delta(1,1)*rho
! 	    s1_dDV_V = sm(1)*dDeldV(1,1)*rho
! 	    root(iCyc) = dsqrt(4D0*s1_Delta_V + 1D0)
!     
! 	    Xs(iCyc,1) = 2D0/(1D0 + root(iCyc))
! 	    
! 	    dXs_dV(iCyc,1) = Xs(iCyc,1)**2 /root(iCyc)*(s1_Delta_V*rho - s1_dDV_V)
! 	    
! 	    !esta hay que revisarla... por alguna razón me da mal.
! 	    d2Xs_dV2(iCyc,1) = (2*dXs_dV(iCyc,1)**2 * (Xs(iCyc,1)*root(iCyc) + 1D0)/Xs(iCyc,1)**2            & 
! 	    
! 	                       - Xs(iCyc,1)**2 *rho*(sm(1)*d2DeldV2(1,1) - 2*s1_dDV_V + 2*s1_Delta_V*rho))/root(iCyc)
! 
! 	  else
! 	    
! 	    s1_Delta_V = 0
! 	    s1_dDV_V = 0
! 	    root(iCyc) = 1	    
! 	    
! 	  endif
! 	  sm_Xs(iCyc,1) = sm(1)*Xs(iCyc,1)
! 	  
! 	elseif (NST == 2 .AND. Delta(1,1) <= 0 .AND. Delta(2,2) <= 0) then
! 	  
! !         Two cross associating sites, with same or different mole amounts.
! 	  if (Delta(1,2) > 0 .AND. sm(major) == sm(3-major) ) then
!     
! 	    s1_Delta_V = sm(3-major)*Delta(1,2)*rho
! 	    s2_Delta_V = sm(major)*Delta(1,2)*rho
! 	    b_aux(iCyc) = 1D0 + s1_Delta_V - s2_Delta_V
! 	    root(iCyc) = dsqrt( b_aux(iCyc)*b_aux(iCyc) + 4D0*s2_Delta_V )
!     
! 	    Xs(iCyc,major) = 2D0/(b_aux(iCyc) + root(iCyc))
!     
! 	    
! 	    aux = b_aux(iCyc)*(1D0 - b_aux(iCyc)) + 2*sm(major)*dDeldV(1,2) - 2D0*s2_Delta_V
! 	    dXs_dV(iCyc,major) = -Xs(iCyc,major)**2 * (1D0 - b_aux(iCyc) + aux/root(iCyc))*rho/2D0
! 	    d2Xs_dV2(iCyc,major) = Xs(iCyc,major)**2 * (2D0*dXs_dV(iCyc,major)**2/Xs(iCyc,major)**3                    &
! 	                                                - ( 2*(b_aux(iCyc) - 1D0) - aux**2/root(iCyc)**3               &
! 	                                                    + ((b_aux(iCyc) - 1D0)**2 - 2*aux)/root(iCyc) )*rho*rho/2  &
! 	                                                - sm(major)*d2DeldV2(1,2)*rho/root(iCyc))
! ! 	    if (sm(major) == sm(3-major)) then
! 	      
! 	      Xs(iCyc,3-major) = Xs(iCyc,major)
! 	      dXs_dV(iCyc,3-major) = dXs_dV(iCyc,major)
! 	      d2Xs_dV2(iCyc,3-major) = d2Xs_dV2(iCyc,major)
! 	
! ! 	For any strange reason, this calculation fails for different sites
! 		
! ! 	    else	
! ! 	    
! ! 	      Xs(iCyc,3-major) = 1D0/(1D0 + s2_Delta_V*Xs(iCyc,major))
! ! 	      dXs_dV(iCyc,3-major) = Xs(iCyc,3-major)**2 * ( s2_Delta_V*(Xs(iCyc,major)*rho - dXs_dV(iCyc,major)) - sm(major)*Xs(iCyc,major)*dDeldV(major,3-major)  )
! ! 	      d2Xs_dV2(iCyc,3-major) = 2*dXs_dV(iCyc,3-major)**2/Xs(iCyc,3-major)  &
! !                                      + Xs(iCyc,3-major)**2 *sm(major)*rho*(-d2Xs_dV2(iCyc,major)*Delta(major,3-major)      &
! !                                                                            +2*dXs_dV(iCyc,major)*dDeldV(major,3-major)     &
! !                                                                            -2*dXs_dV(iCYc,major)*Delta(major,3-major)*rho  &
! !                                                                            +Xs(iCyc,major)*d2DeldV2(major,3-major)         &
! !                                                                            -2*Xs(iCyc,major)*dDeldV(major,3-major)*rho   )
! ! 	    
! ! 	    endif
! 	    sm_Xs(iCyc,major) = sm(major)*Xs(iCyc,major)		
! 	    sm_Xs(iCyc,3-major) = sm(3-major)*Xs(iCyc,3-major)		
! 	    
! 	  else
! 	    
! 	    b_aux(iCyc) = 1
! 	    root(iCyc) = 1
! 	    sm_Xs(iCyc,:NST) = sm(:NST)	    
! 	    
! 	  endif
! 	
! 	else
	  
! 	  General case for NST sites with any associating scheme:
	  
!	  Iniciación de las fracciones no asociadas por continuacion. La extrapolación
!	  está acotada ente 0 y 1:
	  if (rhOld(iCYc) > 0) then
	    do i = 1,NST
        
	      aux = XsOld(iCyc,i) + dXs_dVOld(iCyc,i)*(1.D0/rho - 1.D0/rhOld(iCyc))
	      if (aux > 1D0) then
	    
	        Xs(iCyc,i) = 1D0
	    
	      elseif (aux <= 0.D0) then
	    
	        Xs(iCyc,i)= Xs(iCyc,i)/5d0
	    
	      else
	    
	        Xs(iCyc,i) = aux
	    
	      endif
	    enddo
	  else
!  Inicio de las fracciones no asociadas por 3 iteraciones de sustición directa
	    Xs(iCyc,:NST) = 1
	    do k = 1,3
	      do j = 1,NST
        
	        Xs(iCyc,j) = 1D0/(1D0 + sum(rho*sm(:NST)*Xs(iCyc,:NST)*Delta(:NST,j)) )
        
	      enddo
	    enddo	  
	  endif
!  Llamado a calculo en para obtener fracciones no asociadas.
	  Xs_aux = Xs(iCyc,:NST)
	  call OptiNewton (maxIt, NST, sm_aux, Delta_aux, rho, Xs_aux, H_aux, pivot_vector, N_iter_Xs)
	  Xs(iCyc,:NST) = Xs_aux	  
! 	  call Optim (NST, iAct, Xs(iCyc,:NST), SAFTQfunction, Qsp, g, H, 1, tol)
! 	  
!  Cálculo de la derivada de la fracción no asociada respecto VOLUMEN
!
!	  Number of moles of non-bonded sites
	  sm_Xs(iCyc,:NST) = sm(:NST)*Xs(iCyc,:NST)
	  
!	  Vector -dg/dv provisionally stored in dXs_dV:
	  dXs_dV(iCyc,:NST) = -sm(:NST)*matmul( Delta(:NST,:NST)*rho - dDeldV(:NST,:NST), sm_Xs(iCyc,:NST) )*rho
	  call LUBksb (H_aux, NST, NST, pivot_vector, dXs_dV(iCyc,:NST))

	  do i = 1, NST
        
!       Vector -h provisionally stored in d2Xs_dV2 
	    d2Xs_dV2(iCyc,i) = sm(i)*rho*sum(sm_Xs(iCyc,:NST)*(d2DeldV2(:NST,i) - 2*dDeldV(:NST,i)*rho + 2*Delta(:NST,i)*rho**2 )  &

	                                     + 2*sm(:NST)*dXs_dV(iCyc,:NST)*(dDeldV(:NST,i) - Delta(:NST,i)*rho) )               &

	                       - 2*sm(i)*dXs_dV(iCyc,i)**2/Xs(iCyc,i)**3
	  enddo
	  call LUBksb (H_aux, NST, NST, pivot_vector, d2Xs_dV2(iCyc,:NST))   
	  H(iCyc,:NST,:NST) = H_aux
	  indx(iCyc,:NST) = pivot_vector

!
! 	endif

	sum1(iCyc) = Sum_S - sum(sm_Xs(iCyc,:NST))
	sum2 = sum(sm(:NST)*dXs_dV(iCyc,:NST))*T/P

	aux = 1D0 - V_dLng_dV
	Zassoc(iCyc) = -aux * sum1(iCyc) /2
! 	dFVas(iCyc) = -Zassoc*rho
	dZassoc_dZ(iCyc) = -v_dLng_dv*sum1(iCyc)/(Z - B2)/2 + aux*sum2/2
	d2Zassoc_dZ2(iCyc) =  (V_dLng_dV/(Z - B2)*sum1(iCyc) + V_dLng_dV*sum2 )/(Z - B2) + aux*sum(sm(:NST)*d2Xs_dV2(iCyc,:NST))*T**2/P**2/2
	
! 	Store variables
	rhOld(iCyc) = rho
	XsOld(iCyc,:NST) = Xs(iCyc,:NST)
	dXs_dVOld(iCyc,:NST) = dXs_dV(iCyc,:NST)
	
! 	Output:
	ZassocOut = Zassoc(iCyc)
	dZassoc1 = dZassoc_dZ(iCyc)
	dZassoc2 = d2Zassoc_dZ2(iCyc)
	
	return
endsubroutine Z_assoc

!----------------------------------------------------------------------------------------------
!
!  INDAT reads in components and selects proper parameters from parameter tables
!  To change from 15 to K components:  
!  REPLACE DIM. 15 BY K  AND DIM. 120 BY  K(K+1)/2
!
SUBROUTINE INDAT (iCEq, NC, inputFile, outputFile)

	implicit none
	
	integer, parameter            :: NCData = 15, NCM = 30, NSM = 24
	real(8), parameter            :: eps = epsilon(1.D0)
	
	integer, intent(in)           :: iCEq, NC, inputFile, outputFile
	integer                       :: NONHY, major, NOWN_AIJ
	integer                       :: i, id_site1, id_site2, iEntalp, iOut, iPrt, j, k, L, L1, L2, N, N_usr_CKV, NAD, &
&	                                 nDer, NH, NST, nTemp, M
	integer, dimension(NCM)       :: IV
	integer, dimension(NSM,NCM)   :: sigma

	real(8)                       :: AE, AV, C, C1, C2, cK, conA, conB, EL, OM
	real(8), dimension(NCM)       :: ac, ac0, ac1, acen, bc, o, Pcr, q, Tcr, TSqr, PM, Peneloux, HHA, HHB, HHC, HHD,& 
&	                                 HHE, HHF, HHG
	real(8), dimension(NCData)    :: omega, P, T, HA, HB, HC, HD, HE, HF, HG, MW
	real(8), dimension(NCData,4)  :: cKv
	real(8), dimension(NCM,NCM)   :: b_matrix, usr_CKV = 0.D0
	real(8), dimension(NSM,NSM)   :: beta, eps_R, kappa
	CHARACTER(10)                 :: NAME(NCData), NEW(NCM)
! 	character(20)                    versCEOS, thermoName
! 
! 	common /versSUB/                 versCEOS, thermoName
	
	COMMON /CUB/                     C, C1, C2		
	COMMON /CRIT/                    TCR, PCR, ACEN
	COMMON /STYR/                    N, NH, NDER, NTEMP, NST   
	COMMON /OVER/                    AC, Q, TSQR
	COMMON /PAR/                     CK(NCM*(NCM+1)/2), BC, AC0, AC1, b_matrix
	COMMON /COUT/                    IOUT

	common /GrupAs1/                 kappa, eps_R
	common /GrupAs1bis/              beta
	common /GrupAs2/                 sigma, major	
	
	COMMON /NAME/                    NEW
	common /GCPROM/                  PM, Peneloux, HHA, HHB, HHC, HHD, HHE, HHF, HHG
	common /ENTALP/                  iEntalp
!
!  Components list
!     
	DATA NAME/'CH4       ','C2        ','C3        ','iC4       ','nC4       ','iC5       ','nC5       ','nC6       ', &
&	          'nC7       ','nC8       ','nC9       ','H2O       ','N2        ','CO2'      ,'H2S       '/	    
!
!  Critical temperatures and pressures					 
!
	DATA T,P/190.6,305.4,369.8,408.1,425.2,460.4,469.6,507.4,540.2,568.8,594.6,647.3,126.2,304.2,373.2,                &
&	45.4,48.2,41.9,36.,37.5,33.4,33.3,29.3,27.,24.5,22.8,218.,33.5,72.8,88.2/		
!
!  Acentric factors		  
!
	DATA OMEGA/.008,.098,.152,.176,.193,.227,.251,.296,.351,.394,.440,.344,.04,.225,.1/   
!
!  Binary interaction coefficients					     
!
	DATA CKV  /2*.45,.53,2*.52,6*.5,4*0.,.02,.06,10*.08,3*0.,.12,11*.15,2*0.,.12,.08,2*.07,4*.06,.05,2*.04,2*.03,0.,.12,0./	  
!  
!  Passur-Danner coefficient for ideal gas enthalpy
! 	
	data HA /-5.58114, -0.76005, -1.22301, 13.28660, 29.11502, 27.62342, 27.17183, 32.0356, 30.70117, 29.50114, 28.5654,&
&	        -2.46322, -0.68925, 4.77805, -0.61782/
	data HB /0.5648334, 0.273088, 0.179733, 0.036637, 0.00204, -.031504, -0.002795, -.023096, -.023143, -.022402, -.021654,&
&	        0.477392, 0.253664, 0.114433, 0.238575/
	data HC /-.282973, -.042956, 0.066458, 0.349631, 0.434879, 0.469884, 0.440073, 0.461333, 0.460981, 0.459712, 0.458518, &
&	        -.052512, -.014549, 0.101132, -.024457/
	data HD /0.417399, 0.312815, 0.250998, 0.005361, -.08181, -.098283, -.086288, -.097402, -.098074, -.098062, -.097973, &
&	        0.064594, 0.012544, -.026494, 0.041067/
	data HE /-1.525576, -1.38989, -1.247461, -.298111, 0.072349, 0.102985, 0.081764, 0.103368, 0.104752, 0.104754, 0.404654,&
&	        -.202759, -.017106, 0.034706, -.130126/
	data HF /1.958857, 2.007023, 1.893509, 0.538662, -.01456, -.029485, -.019715, -.030643, -.03134, -.031355, -.031318, &
&	        0.236310, -.008239, -.01314, 0.144852/
	data HG /-.623373, 0.045543, 0.178189, 0.609350, 0.829122, 0.871908, 0.736161, 0.767792, 0.711737, 0.664632, 0.626960, &
&	        -.339830, 0.050052, 0.343357, -.045932/

	data MW /16.04246, 30.06904, 44.09562, 58.1222, 58.1222, 72.14878, 72.14878, 86.17536, 100.20194, 114.22852, 128.2551, &
&	        18.01528, 28.0134, 44.0095, 34.08088/

! 	versCEOS = '0.3.1.1'
! 
!
!  Specify component list, printing option and equation  
	READ (inputFile, *) (IV(I), I = 1, NC), IPRT
	C = dble(ICEQ)
	iOut = outputFile
!					     
!     EQUATION OF STATES USED IS:   
!					     
!     P = RT/(v - b) - a/(v (v + b) + c b (v - b))
!					     
!     (C=0) CORRESPONDS TO SRK. FOR PR, USE C=1				 
!		
!     (LibreOffice Math code: P = {R T} over (v - b) - a over {v (v + b) + c b (v - b)} 
!	
	if(iCEq == 1)then
	  
	  WRITE (*, 6)
! 	  thermoName = "PR-EOS"
	  IF (outputFile >= 1) WRITE (outputFile, 6)	
6	  FORMAT(//,' Peng-Robinson equation',//)		  

	elseif(iCEq == 0)then
	  
	  WRITE (*, 8)
! 	  thermoName = "SRK-EOS"
	  IF(outputFile >= 1)WRITE (outputFile, 8)	 
8	  FORMAT(//,' Soave-Redlich-Kwong Equation',//)		
	 
	endif
! 
!  Get constats for parameters "a" and "b", depending of the CEOS used. This could
!  be useful in the future, for different CEoS than PR and SRK.
! 
9	CALL CONST(C, C1, C2, CONA, CONB) 
	N = NC	
1	FORMAT(16I5)
	NH = 0
	DO 10 I = 1, N			 
	
	  L=IV(I)
	  IF(L <= NCData)then
	  
	    IF (NH == 0 .AND. L > 11) NH = I
	    NEW(I)=NAME(L)    
	    TCR(I) = T(L)
	    PCR(I) = P(L)
! 	    OM = OMEGA(L)
	    acen(i) = omega(L)
	    PM(i) = MW(L)
	    HHA(i) = HA(L)
	    HHB(i) = HB(L)
	    HHC(i) = HC(L)*1e-3
	    HHD(i) = HD(L)*1e-6
	    HHE(i) = HE(L)*1e-10
	    HHF(i) = HF(L)*1e-14
	    HHG(i) = HG(L)
	    Peneloux(i) = 0
	
	  else
	    
2	    read (inputFile, *) NEW(I)
	    read (inputFile, *) TCR(I), PCR(I), ACEN(I)
	    read (inputFile, *) PM(I), Peneloux(i), HHA(I), HHB(I), HHC(I), HHD(I), HHE(i), HHF(i), HHG(i)
! 	    OM = ACEN(I)
3	    FORMAT(A10, 3F10.4)		 

	  endif
4	  TSQR(I) = 1.D0/DSQRT(TCR(I))	
	  BC(I) = CONB*TCR(I)/PCR(I)           ! BC = b/R
	  AC(I) = CONA*TCR(I)/DSQRT(PCR(I))    ! AC = sqrt(ac/R)
!
!  Temperature dependence selected according to SRK			
!
	  Q(I) = 0.48 + acen(i)*(1.574 - .176*acen(i))   
! 	  q(i) =  0.48508 + acen(i)*(1.55171 - 0.15613*acen(i))   
	  IF (ICEQ == 1) Q(I) = .37464 + acen(i)*(1.54226 - .26992*acen(i))		     
	  
10	enddo
	IF (IPRT >= 0) then
	
	  WRITE (*, 11)			 
	  write (outputFile, '("  ",100("-"))')
	  IF(outputFile >= 1) WRITE (outputFile, 11)
11	  FORMAT(X,'Pure component properties',//,3X,'Compound',6X,'Tc(K)   Pc(atm)   omega  ac(atm cm6/mol2)  b(cm3/mol)   M-value    M(g/mol)  c(cm3/mol)'/)

	  DO 111 I = 1, N	
	  
	    WRITE (*, 12) NEW(I), TCR(I), PCR(I), acen(i), (82.05*ac(i))**2, 82.05*bc(i), Q(I), PM(i), Peneloux(i)
	    IF (outputFile >= 1) WRITE (outputFile, 12) NEW(I), TCR(I), PCR(I), acen(i), (82.05*ac(i))**2, 82.05*bc(i), Q(I), PM(i), Peneloux(i)
	    
111	  enddo
	iEntalp = NC
	do i = 1,NC
	  if ((HHA(i) == 0) .AND. (HHB(i) == 0) .AND. (HHC(i) == 0) .AND. (HHD(i) == 0) .AND. (HHE(i) == 0) .AND. (HHF(i) == 0)) then

	    iEntalp = iEntalp - 1

	  endif
	enddo
	if (iEntalp /= 0) then

	  write (outputFile, '(/," Passut-Danner constant for ideal gas enthalpy calculation:",/,"   Compound   A(BTU/lb)  B(BTU/(lb ºR))  ", &
     &             "C(BTU/(lb ºR2)·1E3  D(BTU/(lb ºR3))·1E6  E(BTU/(lb ºR4))·1E10  F(BTU/(lb ºR5))·1E14     G(BTU/(lb ºR))",/,135("-"))')
	  do i = 1,NC

	    write (outputFile, 14)New(i),HHA(i),HHB(i),1.D3*HHC(i),1.D6*HHD(i),1D10*HHE(i),1D14*HHF(i),HHG(i)

	  enddo

	endif

	endif
12	FORMAT(3X, A10, X, 2(F8.2,2X), F6.4, 4X, G12.4, F12.1, 5X, F7.4, 4X, 2(F7.2,4X))
14	format(3X,A10,X,F9.6,3X,F9.6,10X,F9.6,12X,F9.6,13X,F9.6,13X,F9.6,12X,F9.6)
! 
!  Read user specified kij:
	read (inputFile, *) N_usr_CKV
	if (N_usr_CKV > 0) then
	  do k = 1, N_usr_CKV
	  
	    read (inputFile, *)i, j, EL
	    do l = 1, N
	      if(i == IV(l))then
	        do m = l, N
	          if(j == IV(m))then

	            usr_CKV(l,m) = EL
	            usr_CKV(m,l) = EL
	            
	          endif
	        enddo
	      endif
	    enddo
	    
	  enddo
	endif

15	NAD = 0
	DO I = 1, N
	  DO J = I, N
	
	    NAD = NAD + 1	
	    EL = 0.D0
	    L1 = IV(I)  
	    L2 = IV(J)
	    IF(L1 <= NCData .AND. L2 <= NCData)then
	
!  Compounds included in the internal database
	      IF (NH /= 0 .AND. J >= NH) EL = CKV(L1,L2-11)

	    endif
	    if(N_usr_CKV > 0)then
!  User specified parameters, now allows overwrite a non-zero value 
!  in database with zero from user

	      EL = usr_CKV(i,j)
	
	    endif
20	    CK(NAD)=EL
	  
	  enddo
	enddo

	IF ((IPRT > 0 .AND. NH /= 0) .OR. N_usr_CKV /= 0) then

!  Database interaction parameter printing
	  WRITE (*, 30) (NEW(K), K = 1, N)    
	  IF(outputFile >= 1) WRITE (outputFile, 30) (NEW(K), K = 1, N)
30	  FORMAT(/, '  Binary interaction coefficients ', //, 13X, <N>(A10, 2X))
!
	  do i = 1, N
	  
	    L1 = 1 + (i - 1)*N - ((i*i - i)/2 - i + 1)
	    L2 = i*N - (((i + 1)**2 - i - 1)/2 - i)
	    M = i - 1
	    IF (outputFile >= 1) WRITE (outputFile, 40) new(i), CK(L1:L2)
	    write (*, 40) new(i), CK(L1:L2)
	    
	  enddo
40	  FORMAT(1X, A10, 2X, <12*(i-1)>X, <N>(F6.4, 6X))
	
	endif
	
! 	Read the association contribution, if any.
 
	read (inputFile, *) NST, NOWN_AIJ !   Total number of associating sites in the mixture, and number of non-zero values
	if (NOWN_AIJ == 0) then
	  
	  NST = 0    !By default, all interactions are zero. If the user doesn't specify any, it's like there's no association.
	  
	endif
	if (NST > 0) then
	  eps_R(:NST,:NST) = 0.D0
	  beta(:NST,:NST) = 0.D0
! 	  kappa(:NST,:NST) = 0.D0
	  do i = 1, NOWN_AIJ
	    
	    read (inputFile, *) id_site1, id_site2, AE, AV
	    eps_R(id_site1,id_site2) = AE
	    eps_R(id_site2,id_site1) = AE
	    beta(id_site1,id_site2) = AV
	    beta(id_site2,id_site1) = AV
	    
	  enddo
! 	  Association energy: 
	  write (outputFile, '(//, " Association energy matrix, [epsilon/R] (K)", /)')
	  do i = 1, NST
	  
	    IF (outputFile >= 1) WRITE (outputFile, 41) i, eps_R(i,i:NST)
	    write (*, 41) i, eps_R(i,i:NST)
	    
	  enddo
41	  FORMAT(1X, I2, 2X, <12*(i-1)>X, <NST>(F6.1, 6X))	  
! 	  Association volume:  
	  write (outputFile, '(//, " Association volume matrix, [beta]", /)')
	  do i = 1, NST
	  
	    IF (outputFile >= 1) WRITE (outputFile, 42) i, beta(i,i:NST)
	    write (*, 42) i, beta(i,i:NST)
	    
	  enddo	  
42	  FORMAT(1X, I2, 2X, <12*(i-1)>X, <NST>(F6.4, 6X))
	  write (outputFile, '(//, " Site configuration of each compound:", /)')
	  do i = 1, N
	    
!	    Read "site contribution" of every compound:        
	    read (inputFile, *) sigma(:NST,i)
	    write (outputFile, '(X, A10, X, <NST>(I2, 4X))') new(i), sigma(:NST,i)
	    
	  enddo	
	endif
	

	RETURN				
	
ENDsubroutine INDAT				   
!					     
!     CALCULATE EQUATION CONSTANTS WHICH REPRODUCE CORRECT TC AND PC    
!					     
SUBROUTINE CONST (C, C1, C2, AC, BC)

	IMPLICIT REAL*8(A-H,O-Z)
	
	U = 1D0 + C
	W = -C
	D = U*U - 4D0*W
	C2 = (U + DSQRT(D))/2D0
	C1 = W/C2	     
	S1 = 1D0 + ((U + W)*(U + 3D0) - W)/2D0
	S2 = 1.D0 + U + W			    
	R = DSQRT(S1*S1 - S2*S2*S2)	 
	A = (S1+R)**(1.D0/3D0)
	B = S2/A
	X = A + B + 1D0
	
	BC = 1.D0/(3D0*X + U - 1D0)		   
	AC = BC*DSQRT(X*(X*X - 3D0*W) - U*W)  
	
	RETURN
	
ENDsubroutine CONST
!-------------------------------------------------------------------------------
! 
!  Temset calculates the temperature dependent part of the EoS PURE COMPONENT
!  parameters, here SQRT(A/RT)
!
!  The results, the vectors AC0, AC1, AND AC2 represent value, 1st derivative
!  and 2nd derivative of SQRT(a/RT)
! 
!  Since version 0.1.0, here it is calculated the temperature dependent association
!  strenght
!
!     
SUBROUTINE TEMSET (T)

!
	IMPLICIT REAL*8 (A-H, O-Z)
	
	integer, parameter :: NCM = 30, NSM = 24
	
	COMMON /CRIT/         TCR(NCM), PCR(NCM), omega(NCM)  
	COMMON /STYR/         N, NH, NDER, NTEMP, NST    
	COMMON /OVER/         AC(NCM), Q(NCM), TSQR(NCM)					   
	COMMON /PAR/          CK(NCM*(NCM+1)/2), BC(NCM), AC0(NCM), AC1(NCM), b_matrix(NCM,NCM)

	integer            :: sigma
	real(8)            :: kappa

	common /GrupAs1/      kappa(NSM,NSM), eps_R(NSM,NSM)
	common /GrupAs1bis/   beta(NSM,NSM)
	common /GrupAs2/      sigma(NSM,NCM), major	
	common /GrupAs4/      Delta_T(NSM,NSM), dDeldT_T(NSM,NSM) !Temperature dependent-only part of Delta
	common /GrupAs4bis/   Delta(NSM,NSM), dDeldT(NSM,NSM)	
	
	DATA TOLD/0.D0/		   
	
	IF (T /= TOLD) then !RETURN	   
	
	  TOLD=T				
	  SQT = DSQRT(T)			
	  T2R = 1.D0/T/2.D0
	  DO 10 I = 1, N			 
	
	    ALF = TSQR(I)			 
    !   X = 1.D0 - ALF*SQT			 
	    Q1 = AC(I)*(1.D0 + Q(I))/SQT	   
	    AC0(I) = Q1 - AC(I)*Q(I)*ALF  ! sqrt( a(i)/T ) / R
	    AC1(I) = -Q1*T2R              ! d(AC0(i))/dT
	  
10    enddo
    endif

	if (NST > 0) then
	  
	  Delta_T(:NST,:NST) = 0.D0
	  Delta(:NST,:NST) = 0.D0
	  do i = 1,NST	  
	 
	    do j = i,NST
        
	      if (beta(i,j) > 0 .AND. eps_R(i,j) > 0) then
	
	        exponential = dexp(eps_R(i,j)/T)
	        Delta_T(i,j) = beta(i,j)*(exponential - 1.D0)
	        Delta_T(j,i) = Delta_T(i,j)
	        if (nTemp == 1) then
          
	          dDeldT_T(i,j) = -beta(i,j)*exponential*eps_R(i,j)/T**2 
	          dDeldT_T(j,i) = dDeldT_T(i,j)
          
	        endif
	        
	      endif
        
	    enddo
	  enddo
	endif

	RETURN
	
ENDsubroutine TEMSET
!-------------------------------------------------------------------------------------
! 
!  CUBGEN calculates fugacity coefficients and their derivatives, using mixture 
!  parameter routine anew.
!  Other modifications are possible by replacing ANEW
!
SUBROUTINE CUBGEN (MT, IC, T, P, Z, XX, FG, FT, FP, FX, aux)
!
!  PARAMETERS:
!
!  MT:     (I):      Phase type desired
!                    1 = LIQ, -1 = VAP, 0 = MIN G PHASE
!  IC:     (O):      Indicator for phase type; ic returns
!                    1:  If called with MT = 1/-1 and phase is found
!                   -1:  If called with MT = 1/-1 and phase is not found
!                    2:  If called with MT = 0 and min G phase is liquid
!                   -2:  If called with MT = 0 and min G phase is vapour
!  T:      (I):      Temperature (K)
!  P:      (I):      Pressure (ATM)
!  Z:      (O):      Calculated compressibility factor
!  XX:     (I):      Mixture total composition
!  FG:     (O):      Log fugacity coefficient vector
!  FT:     (O):      T-derivative of FG
!  FP:     (O):      P-derivative of FG
!  FX:     (O):      Scaled composition derivative of FG
!  AUX:    (O):      various residual properties

	IMPLICIT REAL*8 (A-H, O-Z)
	
	integer, parameter :: NCM = 30, NSM = 24
	
	DIMENSION             XX(NCM), X(NCM), FG(NCM), FT(NCM), FP(NCM), FX(NCM,NCM)		 
	DIMENSION             PD(NCM), AD1(NCM), BD1(NCM), ADT(NCM), AD2(NCM*(NCM+1)/2)
	
	real(8)            :: AUX(6), d2FasdTdn(N), d2Lng_dndm(N,N), d2Lng_dndV(N), dLng_dn(N), &
	                      dPdnas(N), dXs_dn(NSM,N), Ntot
	
! 	Variables in COMMON vector:
	integer            :: sigma
	real(8)            :: dPdn(NCM), d2Zassoc_dZ2, dZassoc_dZ, Pspec, Sum_SmsX, Tspec,  &
	                      x_spec(NCM), Zassoc
	
	COMMON /STYR/         N, NH, NDER, NTEMP, NST    
	COMMON /CUB/          C, C1, C2
	
	common /PTx/          Pspec, Tspec, x_spec
	
	common /ZZZAs/        Xs(2,NSM), sm(NSM), dXs_dV(2,NSM)
	common /ZZZAs2/       Zassoc(2), dZassoc_dZ(2), d2Zassoc_dZ2(2), Sum_SmsX(2)
	common /GrupAs2/      sigma(NSM,NCM), major	
	common /GrupAs3/      sm_Xs(2,NSM), Sum_S, dFVas(2)
	common /GrupAs4/      Delta_T(NSM,NSM), dDeldT_T(NSM,NSM) !Temperature dependent-only part of Delta
	common /GrupAs4bis/   Delta(NSM,NSM), dDeldT(NSM,NSM), dDeldV(NSM,NSM),              &
	                      d2DeldV2(NSM,NSM)
	
	common /GrupAs5/      H(2,NSM,NSM), root(2), b_aux(2), indx(2,NSM)	
	common /phase_ID/     iCyc
	
	common /ScndDer/      dPdV, dPDT, dPdn	
	
	CALL TEMSET(T)		    

!  Calculate molar fractions (x) from number of moles (xx)
	AUX = 0
	nTot = sum(xx(:N))
	x(:N) = xx(:N)/Ntot
	
	Pspec = P; Tspec = T; x_spec(:N) = x(:N)
	
	CALL ANEW (X, A, AT, ATT, B, AD1, AD2, ADT, BD1)
	
	APT = A*P/T			   
	BPT = B*P/T
!  Calculate compressibility factor
	CALL CUBIC (MT, APT, BPT, Z, IC)
	
! 	IC = 1
! 	
! 	TL = 1.D0 - BPT*C - 3.D0*Z
! 	IF (MT*TL < 0.) IC = -1	
!  V/R:	
	V = Z*T/P                           ! V/R
	S1 = 1.D0/(V + C1*B)		    
	S2 = 1.D0/(V + C2*B)		    
	P1 = P/T + A*S1*S2                  !Prep, if there isn't association
	if (NST > 0) then
	  
	  P1 = P1 - Zassoc(iCyc)/V
	  
	endif
	PA = -S1*S2                         !Patt/a = dP/da
	PN = P1                            
	FAC = C1*S1 + C2*S2		   
	P2 = A*PA                           !Patt
	PB = P1*P1 - FAC*P2                 !probably, this is dP/db
	
!  DPDV = d(P/T) / d(V/R)
	DPDV = -P1*P1 - P2*(S1 + S2)
!  
!  Association contribution:
	if (NST > 0) then
	  
	  DPDV = DPDV + P*dZassoc_dZ(iCyc)/T/V - Zassoc(iCyc)/V/V
	  dLng_dn(:N) = 0.475D0*BD1(:N)/(V - 0.475D0*B)
	  if (nDer > 1 .OR. nTemp > 0) then

!  Calculation of the number of associating site moles derivatives of the non-bonded
!  fraction. This is required for the calculation of d2[ln(phi_k)]/[dni dnj], and
!  for d[ln(phi_k)]/dT and dP
	    
	    do i = 1, N
	
	      d2Lng_dndm(:N,i) = dLng_dn(:N)*dLng_dn(i)
	      d2Lng_dndV(i) = -dLng_dn(i)/(V - 0.475D0*B)
	      
!
!	      d(Passoc/T)/dn(i) = - d2F/(dn(i)dV)  Sum_SmsX(iCyc) comes from a COMMON vector.
	      dPdnas(i) = -sum(sigma(:NST,i)*dXs_dV(iCyc,:NST)/Xs(iCyc,:NST))               &
	                  + (d2Lng_dndV(i)*Sum_SmsX(iCyc) - dLng_dn(i)*sum(sm(:NST)*dXs_dV(iCyc,:NST)))/2		
	
! 	      Vector -dg/dn(i) provisionally stored in dXs_dn(i):
! 	      Contribution from (dg/dn)_X=cnt, where "g" is the _gradient_ of Q at stationary point.
	      dXs_dn(:NST,i) = sm(:NST)*dLng_dn(i)*(1D0 - Xs(iCyc,:NST))/Xs(iCyc,:NST)               
! 	      Contribution from (dg/ds)_n=cnt:
	      do j = 1, NST
	
	        if (sm(j) >= 1D-16) dXs_dn(j,i) = dXs_dn(j,i)                                & 
	                            + sm(j)*dot_product(sigma(:NST,i), Xs(iCyc,:NST)*Delta(:NST,j))/V
	  
	      enddo
! 	      Solve for the value of dXs/dn(i):
	      call LUBksb (H(iCyc,:NST,:NST), NST, NST, indx(iCyc,:NST), dXs_dn(:NST,i))
		
	    enddo
	    
	  endif
	  if (nTemp > 0) then

! 	    
!	    dPdTas = Zassoc(iCyc)/V
	    dFasdT = 0D0
	    d2FasdTdV = 0D0
	    d2FasdTdn(:N) = 0.D0
 	    do i = 1, NST  
!         
	      do j = 1, NST

	        Hij = sm_Xs(iCyc,i)*sm_Xs(iCyc,j)*dDeldT_T(i,j)/(V-0.475D0*B)/2 
	        dFasdT = dFasdT - Hij
	        if (Delta(i,j) > 0) d2FasdTdV = d2FasdTdV - Hij*(2*dXs_dV(iCyc,j)/Xs(iCyc,j) + dDeldV(i,j)/Delta(i,j) - 1D0/V)
	        if (sm(j) > 0) d2FasdTdn(:N) = d2FasdTdn(:N) - Hij*(2*dble(sigma(i,:N))/sm(i) + 2*dXs_dn(i,:N)/Xs(iCyc,i) + dLng_dn(:N))
! 	        dPdTas = dPdTas + sm_Xs(iCyc,i)*SUM( sm_Xs(iCyc,:NST)*dDeldT(i,:NST)*(dXs_dV(iCyc,:NST)/Xs(iCyc,:NST) - 5D-1/V - dDeldV(i,:NST)/Delta(i,:NST))/2 )/V
	      enddo
 	     
 	    enddo	    
	    
	  endif
	  
	endif
	AUX(2) = DPDV*T*v/P
	
	XL1 = DLOG(V/(V - B))              !dF
	
!  Coefficients for fugacity coefficient calculation	
	FN = XL1                           !dF/dn keeping every constant (remember, Frep is multiplied by "n")
	                                   !See Eq.(66) in page 88.
	XL2 = DLOG(S1/S2)/(C2 - C1)        !
	FA = -XL2/B                        !-dFatt/da
	F2 = -A*FA                         !Fatt
	FF = FN - F2			    
	GB = -V*PA			    
	F2B = (A*GB - F2)/B		   
	FB = P1 - F2B			   
	FNB = P1
	FAB = -F2B/A			  
	GBB = -GB*FAC	 
	F2BB = (A*GBB - 2.D0*F2B)/B	    
	FBB = P1*P1 - F2BB		    
	XLZ = DLOG(Z)		           !ln(Z)
	FNP = FN - XLZ	
	
!  Fugacity coefficients:
	FG(:N) = FNP + FA*AD1(:N) + FB*BD1(:N) !Eq.(66) in page 88, plus ln(Z)
	                                       ! because the reference here is T and P
	if (NST > 0) then
	  
	  do i = 1, N
	    
	    FG(i) = FG(i) - dLng_dn(i)*Sum_SmsX(iCyc)/2 !Sum_SmsX(iCyc) comes from a COMMON
	    if (maxval(sigma(:NST,i)) > 0) then
	
	      FG(i) = FG(i) + sum(dfloat(sigma(:NST,i))*dlog(Xs(iCyc,:NST)))
 
	    endif
	    
	  enddo
	
	endif

	if(nDer >= 2 .OR. nTemp >= 1)then

! 	  nt·d(lnPhi(i))/dn(j) = nt·d2F/dn(i)dn(j) + 1 + n/RT·(dP/dn(i)·dP/dn(j))/ (dP/dV) [here, R = 1]
!
!                                     _ 
!	  dP/dn(i)/(dP/dV) = dV/dn(i) = V(i)  
! 	  PD(:N) = (PN + PA*AD1(:N) + PB*BD1(:N))/DPDV
	  do i = 1, N
	    
	    dPdn(i) = PN + PA*AD1(i) + PB*BD1(i)
	    if (NST > 0) dPdn(i) = dPdn(i) + dPdnas(i)
	    PD(i) = dPdn(i)/dPdV
	    
	  enddo
	    
	  if (nDer >= 2) then
	  
	    NADR = 0
	    DO 30 I = 1, N			 
	      
! 	      PDI = PD(I)			   
	      CC = 1.D0 + FNB*BD1(I) + PN*PD(I)		  
	      CA = FAB*BD1(I) + PA*PD(I)		    
	      CB = FNB + FAB*AD1(I) + FBB*BD1(I) + PB*PD(I)	 
	      DO 301 J = I, N
	      
	        NADR = NADR + 1			 
	        FX(I,J) = CC + CA*AD1(J) + CB*BD1(J) + FA*AD2(NADR)
	        if (NST > 0) then
	    
	          FX(i,j) = FX(i,j) + dPdnas(j)*PD(i) !contribution not summed before, because of how terms are added.
	          
	          FX(i,j) = FX(i,j) - d2Lng_dndm(i,j)*Sum_SmsX(iCyc)/2
	          if (maxval(sigma(:NST,i)) > 0) then
	
	            FX(i,j) = FX(i,j) + sum(sigma(:NST,i)*dXs_dn(:NST,j)/Xs(iCyc,:NST))

	          endif
	          aux0 = 0D0
	          if (maxval(sigma(:NST,j)) > 0) then
	
	            aux0 = sum(sigma(:NST,j)*(1D0 - Xs(iCyc,:NST)))
	
	          endif
	          FX(i,j) = FX(i,j) - dLng_dn(i)*(aux0 - sum(sm(:NST)*dXs_dn(:NST,j)))/2
	          
	        endif
	        FX(J,I) = FX(I,J)
		  
301	      enddo
30	    enddo

50	  endif   
	  if (nTemp > 0) then
	  
	    DFT = FA*AT
	    if (NST > 0) DFT = DFT + dFasdT
	    
	    HR = -T*DFT + Z - 1.D0
!  Residual enthalpy/R stored in AUX(3)
	    HR = HR*T
	    SR = -T*DFT - FF + XLZ
!  Residual entropy/R stored in AUX(4)
	    AUX(3) = HR
	    AUX(4) = SR  
	    
	    PTR = PA*AT
	    FTT = FA*ATT
!
!       (dP/dT)_Vn = P/T - RT·(d2F/dT/dV)_n  .  See page 63.
	    DPDT = P/T + T*PTR
        if (NST > 0) dPdT = dPdT + T*d2FasdTdV

	    CP = -T*(T*FTT + 2.D0*DFT) - DPDT**2/DPDV - 1.D0
!  Residual heat capacity/R in AUX(5)
	    AUX(5) = CP
	   
	    DVDT = -DPDT/DPDV/T
! AUX(6) is pressure derivative of resid. entropy	    
	    aux(6) = dvdT/v
	    
	    CC = 1.D0/T
	    CB = FAB*AT
	    DO 70 I = 1, N
	
!                          _
!         d(lnPhi(i))/dP = V(i)/RT - 1/P (See Eq. (15)). The association contribution has been already added.
	      FP(I) = -PD(I)/T - 1.D0/P	     
!                                              _
!         d(lnPhi(i))/dP = d2F/dTdn(i) + 1/T - V(i)/RT·dP/dT
	      FT(I) = CC + CB*BD1(I) + FA*ADT(I) + DPDT*PD(I)/T
          if (NST > 0) FT(i) = FT(i) + d2FasdTdn(i)
	    
70	    enddo

	  endif
	
	endif
	RETURN				
ENDsubroutine CUBGEN
!-------------------------------------------------------------------------------------
!
!  Subroutine CubicVirial calculates the predicted second virial coefficient 
!  of a mixture from SRK or PR equations
SUBROUTINE CubicVirial(NC, n, T, Bcoef, dBcdT)

	implicit none
	
	integer, parameter      :: NCM = 30
	real(8), parameter      :: R = 82.05D0
	
	integer                 :: NC, NDER, NH, NTEMP, N1, NST
	
	real(8)                 :: Bcoef, dBcdT, T
	
	real(8)                 :: a, aT, aTT, b, AD2(NCM*(NCM+1)/2)
	real(8), dimension(NCM) :: bd, ad, adT, n

	COMMON /STYR/              N1, NH, NDER, NTEMP, NST  
	
	N1 = NC
	CALL ANEW(n, a, aT, aTT, B, AD, AD2, ADT, bd)
	
	Bcoef = b*R - a/T
	dBcdT = a/T/T - aT/T
	
	return
	
ENDsubroutine CubicVirial

! 
!  This package calculates associating contribution to Helmholtz free energy and its derivatives
!  based on SAFT. 
!
!  Associating parameters are red by subroutine SAFTPar
!
!  Derivative functions are calculated depending of integer variables:
!  
!    v_deriv: derivatives wrt the total volume (0, 1, 2)
!
!    n_deriv: derivatives wrt the number of moles of every component (0, 1, 2)
! 
!    t_deriv: derivatives wrt the temperature of the system (0, 1, 2)
!
!  
! module SAFT !it should be a module in the future

!   implicit none
!   real(8), parameter                                :: R = 8.314D0
!   
!   integer, save                                     :: NST
!   integer                                           :: inputFile, outputFile, n_deriv, t_deriv, v_deriv
!   integer, dimension(:,:),     allocatable, save    :: sigma
!   
!   real(8), dimension(:),       allocatable          :: d2Xs_dT2, d2Xs_dV2, dXs_dT, dXs_dV, sm, Xs
!   real(8), dimension(:,:),     allocatable          :: dXs_dn, dXs_ds, kappa, eps_R
!                                
!                                
!   integer,                                  private :: i, j, k, l, m, n, m1, m2, maxSiteN, id_group1, id_group2, id_site1, id_site2
!   integer, dimension(:),       allocatable, private :: pivot, Massoc
!   integer, dimension(:,:),     allocatable, private :: nyAss
!                                
!   real(8),                                  private :: Q, Sum_S
!   real(8), dimension(:),       allocatable, private :: g, dg_dT, dg_dV, sm_Xs
!   real(8), dimension(:,:),     allocatable, private :: C, d2CdV2, d2DeldTV, d2DeldV2, dCdT, dCdV, dDeldT, dDeldV, Delta, H
!   real(8), dimension(:,:,:),   allocatable, private :: d2DeldTn, d2DeldVn, dDeldn
!   real(8), dimension(:,:,:,:), allocatable, private :: assocEnergy, assocVol, d2Deldnm
!   
!   save
!                                             public  :: SAFT_Helmholtz, SAFTPar, Strenght
!                                                 
!   contains
!
!
!--------------------------------------------------------------------------
!
!
!
!
!
subroutine SAFT_Helmholtz (model, NC, NST, n_deriv, v_deriv, t_deriv, T, V, n, &
                           Fas, dFV, dFVV, dFTV, dFVn, dFn, dFnm, dFTn, dFT, dFTT)

	implicit none
	
	integer,                   intent(in)          :: n_deriv, NC, NST, model, t_deriv, v_deriv
	real(8),                   intent(in)          :: T, V
	real(8), dimension(NC)   , intent(in)          :: n
	
	real(8),                   intent(out)         :: dFT, dFTT, dFTV, dFV, dFVV, Fas
	real(8), dimension(NC),    intent(out)         :: dFn, dFTn, dFVn
	real(8), dimension(NC,NC), intent(out)         :: dFnm 
	 

! 	Internal variables:
	integer                                        :: i, j, k, l
	
	integer, dimension(NST,NC)                     :: sigma
	
	real(8)                                           help, invV, s2X2_2V
	real(8), dimension(NST)                        :: dXs_dT, dXs_dV, sm, Xs
	real(8), dimension(NC,NST)                     :: dXs_dn
	real(8), dimension(NST,NST)                    :: d2DeldT2, d2DeldTV, d2DeldV2, dDeldT, dDeldV, Delta
	real(8), dimension(NC,NST,NST)                 :: d2DeldTn, d2DeldVn, dDeldn
	real(8), dimension(NC,NC,NST,NST)              :: d2Deldnm
	
!   	Get site non-bonded fractions:
! 	call NonBondedFraction (model, NC, NST, v_deriv, t_deriv, n_deriv, T, V, n, &
! 	                        sigma, sm, Delta, dDeldV, d2DeldV2, dDeldT, d2DeldT2, &
! 	                        dDeldn, d2DeldVn, d2DeldTn, d2Deldnm, Xs, dXs_dV, dXs_dT, dXs_dn)
	call NonBondedFraction (model, NC, NST, v_deriv, t_deriv, n_deriv, T, V, n, &
	                        sigma, sm, Delta, dDeldV, d2DeldV2, dDeldT, d2DeldT2, d2DeldTV, &
	                        dDeldn, d2DeldVn, d2DeldTn, d2Deldnm, Xs, dXs_dV, dXs_dT, dXs_dn)	
	Fas = 0
	dFV = 0
	dFVV = 0
	dFT = 0
	dFTV = 0
	dFTT = 0
	dFn(:NC) = 0
	dFTn(:NC) = 0
	dFVn(:NC) = 0
	dFnm(:NC,:NC) = 0
	invV = 1D0/V
	do k = 1, NST
	  if (sm(k) > 0) then
	    
	    Fas = Fas + sm(k)*(dlog(Xs(k)) - Xs(k)/2) + sm(k)/2
	    do l = 1, NST
	      if (sm(l) > 0) then
		  
	        s2X2_2V = sm(l)*Xs(l) * sm(k)*Xs(k) * invV/2
	        if (v_deriv > 0) then
	          
	          help = s2X2_2V*(dDeldV(l,k) - Delta(l,k)*invV)
	          dFV = dFV - help
	          if (v_deriv > 1) dFVV = dFVV + 2*help*invV - 2*help*dXs_dV(l)/Xs(l) - s2X2_2V*d2DeldV2(l,k)
	          
	        endif
	        if (t_deriv > 0) then
	        
	          help = s2X2_2V*dDeldT(l,k)
	          dFT = dFT - help
	          dFTV = dFTV + help*invV - 2*help*dXs_dV(l)/Xs(l) - s2X2_2V*d2DeldTV(l,k)
	          if (t_deriv > 1) dFTT = dFTT - 2*help*dXs_dT(l)/Xs(l) - s2X2_2V*d2DeldT2(l,k)
	        
	        endif
	        if (n_deriv > 0) then
	        
	          do i = 1, NC
	        
	            help = s2X2_2V*dDeldn(i,l,k)
	            dFn(i) = dFn(i) - help
	            dFVn(i) = dFVn(i) + help*invV - 2*help*dXs_dV(l)/Xs(l) - s2X2_2V*d2DeldVn(i,l,k)
	            if (t_deriv > 0) dFTn(i) = dFTn(i) - s2X2_2V*(d2DeldTn(i,l,k) + 2*dXs_dT(l)/Xs(l)*dDeldn(i,l,k))
	            if (n_deriv > 1) then
	              
	              dFnm(:NC,i) = dFnm(:NC,i) - 2*help*dXs_dn(:NC,l)/Xs(l) - s2X2_2V*d2Deldnm(:NC,i,l,k)

	            endif
			
	          enddo
	          
	        endif
	        
	      endif
	    enddo
	  endif
	  if (n_deriv > 0) then
	    do i = 1, NC
	  
! 	      if (sigma(k,i) > 0) then
	  
	        dFn(i) = dFn(i) + sigma(k,i)*dlog(Xs(k))
		  dFVn(i) = dFVn(i) + sigma(k,i)*dXs_dV(k)/Xs(k)
	        if (t_deriv > 0) dFTn(i) = dFTn(i) + sigma(k,i)*dXs_dT(k)/Xs(k)
	        if (n_deriv > 1) then
	
	          dFnm(:NC,i) = dFnm(:NC,i) + sigma(k,i)*(dXs_dn(:NC,k)/Xs(k)) - sigma(k,:NC)*Xs(k)*invV*sum(sm(:NST)*Xs(:NST)*dDeldn(i,:NST,k) )
	
	        endif
	  
! 	      endif
	  
	    enddo
	  endif
	enddo
	return
endsubroutine SAFT_Helmholtz
!-------------------------------------------------------------------------------------------------
!
!  Subroutine Strenght calculates different parameters required for the model, and their temperature
!  dependence
!
subroutine Strenght (model, NC, NST, t_deriv, v_deriv, n_deriv, T, V, n,&
	               Delta, dDeldT, d2DeldT2, d2DeldTV, d2DeldTn,&
	               dDeldV, d2DeldV2, d2DeldVN, dDeldn, d2Deldnm)

	implicit none
	integer, parameter                             :: NCM = 30, NSM = 24
	
	integer, intent(IN)                            :: model, NC, NST, n_deriv, t_deriv, v_deriv
	real(8), intent(IN)                            :: T, V, n(NC)
	                                                                                            
	integer                                        :: i_start
	real(8), dimension(NST,NST), intent(out)       :: d2DeldT2, d2DeldTV, d2DeldV2, dDeldT, dDeldV, Delta 
	real(8), dimension(NC,NST,NST)   , intent(out) :: dDeldn, d2DeldTn, d2DeldVN
	real(8), dimension(NC,NC,NST,NST), intent(out) :: d2Deldnm

	real(8), dimension(NSM,NSM)                    :: common_eps_R, common_kappa
	common /GrupAs1/                                  common_kappa, common_eps_R
	common /GrupAs2/                                  common_sigma, common_major, common_Massoc	                                                
    
	integer                                        :: A, B, B_start, i, j, k, l, m, common_major
	integer, dimension(NCM)                        :: common_Massoc
	integer, dimension(NSM,NCM)                    :: common_sigma
	
	real(8)                                           invT, T_old
	real(8), dimension(NC,NC)                    :: RDF, dRDF_T, dRDF_T2, dRDF_TV, dRDF_V, dRDF_V2, dRDF_V3
	real(8), dimension(NC,NC,NC)                 :: dRDF_n, dRDF_Tn, dRDF_Vn
	real(8), dimension(NC,NC,NC,NC)              :: dRDF_nm
	DATA T_OLD /0.D0/
  
!	Evaluate the term kappa*(exp(e/RT) - 1) only if the temperature has changed from previous call.
! 	if (T /= T_old) then	  
! 	  
! 	  T_old = T
	  invT = 1D0/T
	  do k = 1, NST
	    do l = 1, NST
	    
	      Delta(l,k) = common_kappa(l,k)*(dexp(common_eps_R(l,k)*invT) - 1)
	      if (t_deriv > 0) then
	    
	        dDeldT(l,k) = -(Delta(l,k) + common_kappa(l,k))*common_eps_R(l,k)*invT*invT
	        if (t_deriv > 1) then
	          
	          d2DeldT2(l,k) = -dDeldT(l,k)*(common_eps_R(l,k) + 2*T)*invT*invT
	          
	        endif
	    
	      endif
	    
	    enddo
	  enddo
	  
! 	endif

!     EoS isn't GCA, so the RDF is a function of mixture volume and number of moles (and temperature if d(i) = f(T)):
	call subRDF (model, NC, NST, n_deriv, t_deriv, v_deriv, T, RDF, dRDF_T, dRDF_T2,&
	          dRDF_TV, dRDF_V, dRDF_V2, dRDF_n, dRDF_Vn, dRDF_Tn, dRDF_nm)
      
	k = 0
	do i = 1, NC
	  do A = 1, common_Massoc(i)
	
	    k = k + 1       !Defino actual sitio "A" como sitio "k"
	    l = 0 !k
	    do j = 1, NC      
! 	      if (common_Massoc(i) == common_Massoc(j)) then
! 		  
! 	        B_start = A !simmetrical submatrix (2B + 2B, 2B + 4B, etc)
! 	        
! 	      else
		  
	        B_start = 1 !unsimmetrical submatrix (2B + 3B...)
	        
! 	      endif  
	      do B = B_start, common_Massoc(j)
	
	        l = l + 1 !defino actual sitio "B" como sitio "l"
	        dDeldV(l,k) = Delta(l,k)*dRDF_V(j,i)
	        d2DeldV2(l,k) = Delta(l,k)*dRDF_V2(j,i)
	        if (t_deriv > 0) then
	          
	          d2DeldT2(l,k) = d2DeldT2(l,k)*RDF(j,i) + 2*dDeldT(l,k)*dRDF_T(j,i) + &
	                          Delta(l,k)*dRDF_T2(j,i) 
	          d2DeldTV(l,k) = dDeldT(l,k)*dRDF_V(j,i) + Delta(l,k)*dRDF_TV(j,i)
	       
	        endif
	        if (n_deriv > 0) then
	       
	          do m = 1, NC
	       
	            dDeldn(m,l,k) = Delta(l,k)*dRDF_n(m,j,i)
	            d2DeldVn(m,l,k) = Delta(l,k)*dRDF_Vn(m,j,i)
	            if (t_deriv > 0) then
	       
	              d2DeldTn(m,l,k) = dDeldT(l,k)*dRDF_n(m,j,i) + Delta(l,k)*dRDF_Tn(m,j,i)
	            
	            endif
	            if (n_deriv > 1) then
		 
	              d2Deldnm(:NC,m,l,k) = Delta(l,k)*dRDF_nm(:NC,m,j,i)			
		 
	            endif
	 
	          enddo  
	    
	        endif    
	        if (t_deriv > 0) dDeldT(l,k) = dDeldT(l,k)*RDF(j,i) + Delta(l,k)*dRDF_T(j,i)
	        Delta(l,k) = Delta(l,k)*RDF(j,i)

	      enddo  
	    enddo
	
	  enddo
	enddo
	

    return
    
  endsubroutine Strenght
!-------------------------------------------------------------------------------------------------

!-------------------------------------------------------------------------------------------------
!
!  Subroutine NonBondedFraction calculates, naturally, the associating site non-bonded fraction at specified
!  T, V and s (number of site moles) of the mixture.
!
!  Derivatives are calculated in a separate subroutine. 
!
subroutine NonBondedFraction (model, NC, NST, v_deriv, t_deriv, n_deriv, T, V, n, &
	                        sigma, sm, Delta, dDeldV, d2DeldV2, dDeldT, d2DeldT2, d2DeldTV, &
	                        dDeldn, d2DeldVn, d2DeldTn, d2Deldnm, Xs, dXs_dV, dXs_dT, dXs_dn)


	integer, parameter                      :: maxIt = 100, NCM = 30, NSM = 24
	real(8), parameter                      :: tol_x = 1.D-14
	                                        
	integer, intent(IN)                     :: n_deriv, NC, NST, model, t_deriv, v_deriv	                                        
	real(8), intent(IN)                     :: T, V, n(NC)
	
	integer, dimension(NST,NC), intent(out) :: sigma	                                        
	real(8), dimension(NST), intent(out)    :: dXs_dT, dXs_dV, sm, Xs
	real(8), dimension(NC,NST), intent(out) :: dXs_dn
	
	common /GrupAs2/                           common_sigma(NSM,NCM), common_major, common_Massoc(NCM)
	
	integer                                 :: common_sigma, common_major, common_Massoc, i, k, l
	integer, dimension(NST)                 :: pivot_vector
	
	real(8)                                 :: invV, rho
	real(8), dimension(NST,NST), intent(out):: d2DeldT2, d2DeldTV, d2DeldV2, dDeldT, dDeldV, &
	                                           Delta
	real(8), dimension(NC,NST,NST), intent(out) :: d2DeldTn, d2DeldVn, dDeldn
	real(8), dimension(NC,NC,NST,NST), intent(out) :: d2Deldnm	
	
	real(8), dimension(NST,NST)             :: H, dg_ds	
	   
	sigma = common_sigma(:NST,:NC)
	sm = matmul(sigma, n) !	Moles of sites
	invV = 1D0/V
! 	rho = sum(n(:NC))*invV
	call Strenght (model, NC, NST, t_deriv, v_deriv, n_deriv, T, V, n, &
	               Delta, dDeldT, d2DeldT2, d2DeldTV, d2DeldTn,&
	               dDeldV, d2DeldV2, d2DeldVN, dDeldn, d2Deldnm)


! 	Initialize non-bonded fraction by 3 direct substitution iterations
	Xs = 1
	do i = 1, 3
	  do k = 1, NST

	    Xs(k) = 1/(1 + sum(sm(:NST)*Xs(:NST)*Delta(:NST,k))*invV)

	  enddo
	enddo

!	Calling second order minimization to get non-bonded site fraction
!	call Optim (iAct, Xs, SAFTQfunction, Q, g, H, Restr, tol_x)
	call OptiNewton (maxIt, NST, sm, Delta, invV, Xs, H, pivot_vector, i)    
!	Number of moles of not-bonded sites:
! 	sm_Xs = sm*Xs
	if (v_deriv > 0) then
	
!     Volume derivatives. Hessian matrix is helpfuly already diagonalized
	  do k = 1, NST
	    
	    !This is -dg/dv
	    if (sm(k) > 0) then
		
	      dXs_dV(k) = sm(k)*sum(sm(:NST)*Xs(:NST)*(dDeldV(:NST,k) - Delta(:NST,k)*invV))*invV
	      
	    else
		
	      dXs_dV(k) = 0
	      
	    endif

	  enddo
!	  First derivative:
	  call LUBksb (H, NST, NST, pivot_vector, dXs_dV)
	
	endif
	
!     Number of moles derivatives:      
	if (n_deriv > 0) then
	         
!       General case, any quantity of sites
	  do k = 1, NST
      
! 	    dg/dn = (dg/ds)n=cte * ds/dn + (dg/dn)s=cte
!	    ds/dn = sigma
	    do i = 1, NC
	
! 	      Store -dg/dn: (note opossite sign)
	      dg_ds(:NST,k) = sm(k)*Xs(:NST)*Delta(:NST,k)*invV
	      dXs_dn(i,k) = dot_product(dg_ds(:NST,k),sigma(:NST,i))
	      if (sm(k) > 0) dXs_dn(i,k) = dXs_dn(i,k) + sm(k)*sum(sm(:NST)*Xs(:NST)*dDeldn(i,:NST,k)*invV)
	
	    enddo
	  
	  enddo
	  do i = 1, NC
	  
	    call LUBksb (H(:NST,:NST), NST, NST, pivot_vector(:NST), dXs_dn(i,:NST))
	
	  enddo
	endif
!     
!	Temperature derivatives  
	
	if (t_deriv > 0) then

	  do k = 1, NST
	    
	    !This is -dg/dT
	    if (sm(k) > 0) then
		
	      dXs_dT(k) = sm(k)*sum(sm(:NST)*Xs(:NST)*dDeldT(:NST,k))*invV
	      
	    else
		
	      dXs_dT(k) = 0
		
	    endif
		
	  enddo	  
	  call LUBksb (H, NST, NST, pivot_vector, dXs_dT)
	    
	endif
	
	return

endsubroutine NonBondedFraction
!---------------------------------------------------------------------------------------------
!
!  SAFTPar reads associating parameters for the SAFT contribution to EoS's. It is common for
!  most people who work with SAFT models to specify assoc. interaction param. in a 4 index way. 
!  That is: site k in component A with site l in component B. However, this isn't necessary from
!  the calculations perspective and these 4 index matrixs are converted here into single 2 index 
!  matrixs.
!  In case of GCA, the subroutine reads interaction param. between sites in groups/segments
!
subroutine SAFTPar (inputFile, outputFile, model, NC, NST, CName, d03)

	implicit none
  
	integer, parameter  :: NCM = 30, NSM = 24
	
	integer, intent(in) :: inputFile, model, NC, NST, outputFile
	real(8), intent(in) :: d03(NC,NC)
	character(10), dimension(NC), intent(in) :: CName
	
	
	integer             :: A, B, i, j, k, l, major, id_site1, id_site2, NOWN_AIJ
	real(8)             :: AE, AV

	common /GrupAs1/       kappa, eps_R
	common /GrupAs2/       sigma, major, Massoc	
	integer             :: Massoc(NCM), sigma(NSM,NCM)
	real(8)             :: eps_R(NSM,NSM), kappa(NSM,NSM)
 
!	User specified values
	read (inputFile, *) NOWN_AIJ     
!
!	Read association Energy/R (or energy/k_b) [K] and volume [-]	     
	eps_R(:NST,:NST) = 0.D0
	kappa(:NST,:NST) = 0.D0
	do i = 1, NOWN_AIJ
	  
	  read (inputFile, *) id_site1, id_site2, AE, AV
	  eps_R(id_site1,id_site2) = AE
	  eps_R(id_site2,id_site1) = AE
	  kappa(id_site1,id_site2) = AV
	  kappa(id_site2,id_site1) = AV
	  
	enddo
! 	Association energy: 
	write (outputFile, '(//, " Association energy matrix, [epsilonA/k] (K)", /)')
	do i = 1, NST
	
	  IF (outputFile >= 1) WRITE (outputFile, 41) i, eps_R(i,i:NST)
	  write (*, 41) i, eps_R(i,i:NST)
	  
	enddo
41	FORMAT(1X, I2, 2X, <12*(i-1)>X, <NST>(F6.1, 6X))	  
! 	Association volume:  
	write (outputFile, '(//, " Association volume matrix, [kappa]", /)')
	do i = 1, NST
	
	  IF (outputFile >= 1) WRITE (outputFile, 42) i, kappa(i,i:NST)
	  write (*, 42) i, kappa(i,i:NST)
	  
	enddo	  
42	FORMAT(1X, I2, 2X, <12*(i-1)>X, <NST>(F6.4, 6X))
	write (outputFile, '(//, " Site configuration of each compound:", /)')	
	do i = 1, NC 
	  
! 	  Read "site contribution" of every compound:        
	  read (inputFile, *) sigma(:NST,i)
	  Massoc(i) = 0
	  do k = 1, NST
	  
	    if (sigma(k,i) > 0) Massoc(i) = Massoc(i) +1 !sum(sigma(:NST,i))
		
	  enddo
	  write (outputFile, '(X, A10, X, <NST>(I2, 4X))') CName(i), sigma(:NST,i)	  
	  
	enddo
	
! 	Store of kappa(k,l)*d0(i,j)^3 in kappa(k,l) 
	k = 0
	do i = 1, NC
	  do A = 1, Massoc(i)
	    
	    k = k + 1
	    l = 0	
	    do j = 1, NC
	      do B = 1, Massoc(j)
		  
	        l = l + 1
	        kappa(l,k) = kappa(l,k)*d03(j,i)
	        
	      enddo
	    enddo
	    
	  enddo
	enddo
	
	return
  
endsubroutine SAFTPar
!---------------------------------------------------------------------------------------------
!
!
!  Subroutine SAFTQfunctionSAFT calculates the value of the non-bonded fraction of the associating
!  sites in the mixture.
!
!  Calculations follows the work of Michelsen and Hendrix [9], Michelsen [10] and Mollerup and
!  Michelsen [11]. This is based on the optimization of a function named Q, wrt the fraction of
!  non-bonded sites as independent variables.
!
!  QfuncSAFT then is fed into an optimization subroutine (Optim).
!
!func(NDim, iAct, X, Q, g, H)
! subroutine SAFTQfunction (NST, iAct, Xs, Q, g, H)
!
!	Q = Sum(sm(k)·(ln X(k) - X(k) + 1) - 1/2·Sum(Sum(sm(k)·sm(l)·Delta(k,l)·X(l)/V),l=1,NST),k=1,NST)
!
!	  = Q1 + Q2
!
!	g(k) = sm(k)/X(k) - sm(k) - Sum(sm(k)*sm(l)*Delta(k,l)*X(l)/V),l=1,NST)
!
!	  = sm(k)/X(k) - sm(k) - SUMA
!
!	H(k,l) = -(sm(k) + SUMA)/X(k)*d(k,l) - sm(k)*sm(l)*Delta(k,l)/V
!
!	Donde: sm(k)      = es la cantidad de moles de sitios "k"
!	       X(k)       = es la fracción de sitio "k" no asociada
!	       Delta(k,l) = fuerza de asociación entre sitios "k" y "l"
!	       NST        = número de sitios totales
!	       d(k,l)     = delta de Kronocker, = 1 si "l" = "k", sino = 0.
!	       g(k)       = gradiente de Q en la dirección "k"
!	       H(k,l)     = hessiano k,l
!	nota1: tomo a rho = 1/V donde V = vol total. Asume que son iguales porque los moles están normalizados (creo).
!
!
! 	real(8)  :: Q1 = 0.D0, Q2 = 0.D0
! !
! !  Initialization:
! 	g = 0.0D0
! 	H = 0.0D0
! !
! !	Cálculo de Q1, gradiente y hessiano:
! 	do k = 1, NST
! 	
! 	  do l = 1, NST
! 	
! 	    if (Delta(k,l) > 0D0) then
! 	      H(k,l) = -sm(k)*sm(l)*Delta(k,l)*rho
! !  	    Nota: g(k) aún NO ES el gradiente_k. Fran 18/02/2010.
! 	      g(k) = g(k) - H(k,l)*Xs(l)
! 	      Q2 = Q2 + 5D-1*H(k,l)*Xs(k)*Xs(l)
! 	    endif
! 	
! 	  enddo
! 	  H(k,k) = H(k,k) - (sm(k) + g(k))/Xs(k)
! !	Para evitar elementos nulos en la diagonal.
! 	  if (H(k,k) == 0.0D0)H(k,k) = 1.0D-20
! 	  Q1 = Q1 + sm(k)*(dLOG(Xs(k)) - Xs(k) + 1.D0)
! 	
! !	Ahora sí­, g es el gradiente_k:
! 	  g(k) = sm(k)/X(k) - sm(k) - g(k)
! 	
! 	enddo
! 	
! !  Cálculo de Q:
! 	Q = Q1 + Q2
! !
! 	return
! endsubroutine SAFTQfunction
! !
!---------------------------------------------------------------------------------------------
subroutine subRDF (model, NC, NST, n_deriv, t_deriv, v_deriv, T, &
	             RDF, dRDF_T, dRDF_T2, dRDF_TV, dRDF_V, dRDF_V2, dRDF_n, dRDF_Vn, dRDF_Tn, dRDF_nm)
	
	implicit none
	integer, parameter                           :: NCM = 30, NSM = 24
	                                             
	integer, intent(in)                          :: model, n_deriv, NC, NST, t_deriv, v_deriv
	real(8), intent(in)                          :: T
	
	real(8), intent(out), dimension(NC,NC)       :: RDF, dRDF_T, dRDF_T2, dRDF_TV, dRDF_V, dRDF_V2
	real(8), intent(out), dimension(NC,NC,NC)    :: dRDF_n, dRDF_Vn, dRDF_Tn
	real(8), intent(out), dimension(NC,NC,NC,NC) :: dRDF_nm

	integer                                      :: common_sigma(NSM,NCM), common_major, common_Massoc(NCM)
	real(8)                                      :: common_invV, common_Y, common_GZ, common_GZZ, common_zeta, &
	                                                common_zetaT, common_zetaTT, common_zetaTV, common_zetaV,  &
	                                                common_zetaVV, common_ZM2
	real(8), dimension(NCM)                      :: common_d, common_zetaN, common_zetaTN, common_zetaVN
	real(8), dimension(0:3,NCM)                  :: common_DD, common_DDT, common_DTt
	real(8), dimension(0:3)                      :: common_S, common_ST, common_ST2
	
	common /RDF_sPCSAFT/                            common_GZ, common_GZZ, common_ZM2, common_zeta, common_zetaT, &
	                                                common_zetaV
	
	common /PCSAFT_FV/                              common_d, common_DD, common_DDT, common_DTT, common_S, &
	                                                common_ST, common_ST2, common_invV, common_Y
	                                                
 	common /GrupAs2/                                common_sigma, common_major, common_Massoc	                                                
	 
! 	Internal variables
	integer                                      :: i, j, k, l
	real(8)                                      :: dave, dave2, dRDF_zeta2, dRDF_zeta22, DRDF_ZETA2_ZETA3,    & 
	                                                dRDF_zeta3, dRDF_zeta32, invRDF, help2, help3, invT, invT2,  &
	                                                invV2, Y2, Y3, zeta2, zeta2T, zeta2T2,zeta2TV, zeta2V, &
	                                                zeta2V2, zeta3,zeta3T, zeta3T2, zeta3TV,  zeta3V, zeta3V2
	                                                
	real(8)                                         daveI, daveII, daveIJ, daveJ, daveJJ, daveT, daveTT, di, dj, &
	                                                DRDF_DAVE, DRDF_DAVE2, DRDF_ZETA2_DAVE, DRDF_ZETA3_DAVE, &
	                                                inv_didj, inv_didj2, help4, help5, didj, help6, help7
	
	real(8), dimension(NC)                       :: dD_dT, d2D_dT2, minv, zeta2N, zeta2TN, zeta2VN, zeta3N, &
	                                                zeta3TN, zeta3VN

	
 	RDF = 0
 	dRDF_T = 0; dRDF_T2 = 0; dRDF_TV = 0
 	dRDF_V = 0; dRDF_V2 = 0
 	dRDF_n = 0; dRDF_Tn = 0; dRDF_Vn = 0; dRDF_nm = 0	
	
	if (model == 0) then	  

! 	  Simplified PC-SAFT (CS repulsive term, instead of Mansoori-Leeland expression)
	  do i = 1, NC
	    if (common_Massoc(i) > 0) then
	      do j = 1, NC
	        if (common_Massoc(j) > 0) then
		     
	          RDF(j,i) = 0.5D0*common_ZM2 *common_Y**3
	          invRDF = 1D0/RDF(j,i)
	          dRDF_zeta3 = RDF(j,i)*common_GZ                                   !d(g)/d(zeta3)
	          dRDF_zeta32 = common_GZZ*RDF(j,i) + dRDF_zeta3*dRDF_zeta3*invRDF  !d2(g)/d(zeta3)2
	          dRDF_V(j,i) = dRDF_zeta3*common_zetaV                             !d(g)/d(V) = d(g)/d(zeta3) * d(zeta3)/d(V)
	          dRDF_V2(j,i) = dRDF_zeta3*common_zetaVV + dRDF_zeta32*common_zetaV*common_zetaV!etc...
	          if (t_deriv > 0) then
	          
	            dRDF_T(j,i) = dRDF_zeta3*common_zetaT
	            dRDF_T2(j,i) = dRDF_zeta3*common_zetaTT + dRDF_zeta32*common_zetaT*common_zetaT
	            dRDF_TV(j,i) = dRDF_zeta3*common_zetaTV + dRDF_zeta32*common_zetaT*common_zetaV
	            
	          endif
	          if (n_deriv > 0) then
	          
	            do k = 1, NC
	              
	              zeta3N(k) = common_invV*common_DD(3,k)
	              zeta3VN(k) = -zeta3N(k)*common_invV
	              dRDF_n(k,j,i) = dRDF_zeta3*zeta3N(k)
	              dRDF_Vn(k,j,i) = dRDF_zeta3*zeta3VN(k) + dRDF_zeta32*zeta3N(k)*common_zetaV
	              if (t_deriv > 0) then 
	                
	                zeta3TN(k) = common_DDT(3,k)*common_invV
	                dRDF_Tn(k,j,i) = dRDF_zeta3*zeta3TN(k) + dRDF_zeta32*zeta3N(k)*common_zetaT
	                
	              endif 
	              if (n_deriv > 1) then
	      
	                dRDF_nm(:NC,k,j,i) = dRDF_zeta32*zeta3N(k)*zeta3N(l)
	      
	              endif
	              
	            enddo
	            
	          endif
	          
	        endif
	      enddo
	    endif
	  enddo
	  
	elseif (model == 2) then
	
! 	  Original SAFT term (complete Mansoori-Leeland expression for mixtures of HS

! 	  Auxiliary variables:
	  invV2 = common_invV*common_invV
	  zeta2 = common_S(2)*common_invV
	  Y2 = common_Y*common_Y
	  Y3 = Y2*common_Y
	  zeta2V  = -common_S(2)*invV2 !dzeta_2/dV
	  zeta2V2 = -2*zeta2V*common_invV 
	  zeta3V = -common_S(3)*invV2
	  zeta3V2 = -2*zeta3V*common_invV  	  
	  if (t_deriv > 0) then
	    
! 	    The derivatives from common are wrt (1/T) instead of T
	    invT = 1D0/T
	    invT2 = invT*invT
	    
	    zeta2T = -common_ST(2)*common_invV*invT2
	    zeta2T2 = -2*invT*zeta2T + invT2*invT2*common_ST2(2)*common_invV
	    zeta2TV = -zeta2T*common_invV
	    zeta3T = -common_ST(3)*common_invV*invT2
	    zeta3T2 = -2*invT*zeta3T + invT2*invT2*common_ST2(3)*common_invV !common_ST2(3)*common_invV
	    zeta3TV = -zeta3T*common_invV

	    minv(:NC) = common_d(:NC)/common_dd(1,:NC) !1/m(i) 
	    dD_dT(:NC) = -common_ddT(1,:NC)*minv(:NC)*invT2
	    d2D_dT2(:NC) = -2*invT*dD_dT(:NC) + invT2*invT2*common_DTT(1,:NC)*minv(:NC)
	    
	  endif
	  do i = 1, NC
	    if (common_Massoc(i) > 0) then
	      
	      di = common_d(i)
	      do j = 1, NC
		     
	        if (common_Massoc(j) > 0) then !Only compounds able to associate are of importance
	      
	          dj = common_d(j)
	          didj = di*dj
	          inv_didj = 1D0/(di + dj)
		    inv_didj2 = inv_didj*inv_didj
	          dave = didj*inv_didj ! This is the problem with the "pair interaction"  
	          dave2 = dave*dave                                          ! problem, which makes this term so messy.
	          help2 = 3*dave*zeta2*Y2
	          help3 = 2*dave2*zeta2*zeta2*Y3
	          
	          RDF(j,i) = common_Y + help2 + help3
	          dRDF_zeta2 = 3*dave*Y2 + 4*dave2*zeta2*Y3  !dg/zeta_2
	          dRDF_zeta3 = Y2 + 2*help2*common_Y + 3*help3*common_Y !dg/dzeta_3
	          
	          dRDF_zeta22 = 4*dave2*Y3
	          dRDF_zeta32 = 2*Y3 + 6*help2*Y2 + 12*help3*Y2
	          dRDF_zeta2_zeta3 = 6*dave*Y3 + 12*dave2*zeta2*Y3*common_Y
	          
! 	          Volume derivatives
	          dRDF_V(j,i) = dRDF_zeta2*zeta2V + dRDF_zeta3*zeta3V
	          help2 = dRDF_zeta22*zeta2V + dRDF_zeta2_zeta3*zeta3V
	          help3 = dRDF_zeta2_zeta3*zeta2V + dRDF_zeta32*zeta3V
	          dRDF_V2(j,i) = help2*zeta2V + dRDF_zeta2*zeta2V2 + help3*zeta3V + dRDF_zeta3*zeta3V2
	                         
! 	          Temperature derivatives               
	          if (t_deriv > 0) then
	
	            daveI = dj*inv_didj - dave*inv_didj !didj*inv_didj2
	            daveJ = di*inv_didj - dave*inv_didj !didj*inv_didj2
	            daveT = daveI*dD_dT(i) + daveJ*dD_dT(j)
	            dRDF_dave2 = 4*zeta2*zeta2*Y3
	            dRDF_dave = 3*zeta2*Y2 + dave*dRDF_dave2
			
	            dRDF_T(j,i) = dRDF_zeta2*zeta2T + dRDF_zeta3*zeta3T + dRDF_dave*daveT
	            dRDF_zeta2_dave = 3*Y2 + 8*dave*zeta2*Y3
	            dRDF_zeta3_dave = 6*zeta2*Y3 + 12*dave*zeta2*zeta2*Y3*common_Y
	            dRDF_TV(j,i) = help2*zeta2T + dRDF_zeta2*zeta2TV + help3*zeta3T + dRDF_zeta3*zeta3TV &
	                           + (dRDF_zeta2_dave*zeta2V + dRDF_zeta3_dave*zeta3V)*daveT
	            help4 = dRDF_zeta22*zeta2T + dRDF_zeta2_zeta3*zeta3T
	            help5 = dRDF_zeta2_zeta3*zeta2T + dRDF_zeta32*zeta3T
	            if (t_deriv > 1) then
			  
	              daveII = 2*(didj*inv_didj - dj)*inv_didj2
	              daveJJ = 2*(didj*inv_didj - di)*inv_didj2
	              daveIJ = 2*dave*inv_didj2
	              daveTT = (daveII*dD_dT(i) + daveIJ*dD_dT(j))*dD_dT(i) + daveI*d2D_dT2(i) + &
	                       (daveIJ*dD_dT(i) + daveJJ*dD_dT(j))*dD_dT(J) + daveJ*d2D_dT2(j) 
	              
	              dRDF_T2(j,i) = (help4 + dRDF_zeta2_dave*daveT)*zeta2T + dRDF_zeta2*zeta2T2 + &
	                             (help5 + dRDF_zeta3_dave*daveT)*zeta3T + dRDF_zeta3*zeta3T2 + &
	                             (dRDF_dave2*daveT + dRDF_zeta2_dave*zeta2T + dRDF_zeta3_dave*zeta3T)*daveT + &
	                             dRDF_dave*daveTT	                             
	  
	            endif
		      
	          endif
	          if (n_deriv > 0) then
	            
	            do k = 1, NC
	      
	              zeta2N(k) = common_DD(2,k)*common_invV
	              zeta2VN(k) = -zeta2N(k)*common_invV	              
	              zeta3N(k) = common_DD(3,k)*common_invV
	              zeta3VN(k) = -zeta3N(k)*common_invV
	              dRDF_n(k,j,i) = dRDF_zeta2*zeta2N(k) + dRDF_zeta3*zeta3N(k)
	              dRDF_Vn(k,j,i) = help2*zeta2N(k) + dRDF_zeta2*zeta2VN(k) &
	                                             +help3*zeta3N(k) + dRDF_zeta3*zeta3VN(k)
	              if (t_Deriv > 0) then
			    
	                
	                zeta2TN(k) = -common_DDT(2,k)*common_invV*invT2
	                zeta3TN(k) = -common_DDT(3,k)*common_invV*invT2
	                dRDF_Tn(k,j,i) = help4*zeta2N(k) + dRDF_zeta2*zeta2TN(k) + &
	                                 help5*zeta3N(k) + dRDF_zeta3*zeta3TN(k) + &
	                                 (dRDF_zeta2_dave*zeta2N(k) + dRDF_zeta3_dave*zeta3n(k))*daveT	   
                
	              endif
	              if (n_deriv > 1) then
	                do l = 1, k
	                   	                
	                  help6 = dRDF_zeta22*zeta2N(l) + dRDF_zeta2_zeta3*zeta3N(l)
	                  help7 = dRDF_zeta2_zeta3*zeta2N(l) + dRDF_zeta32*zeta3N(l)
	                  
! 	                  "zeta" variables are linear combinations, so second derivatives wrt n are zero.        
	                  dRDF_nm(l,k,j,i) = help6*zeta2N(k) + help7*zeta3N(k)
	                  dRDF_nm(k,l,j,i) = dRDF_nm(l,k,j,i)
	      
	                enddo
	              endif  
	            enddo
	          endif
	          
	        endif
	      enddo
	    endif
	  enddo
	elseif (model == 1) then
	
! 	  Simplified expression used in CPA and other SAFT derived equations
	
	endif
	
	return
endsubroutine subRDF

! endmodule SAFT

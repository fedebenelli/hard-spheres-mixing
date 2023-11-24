subroutine fugacity_from_gamma_phi (model, NC, phase_type, t_deriv, p_deriv, n_deriv, T, P, n, lnPhi, n_dlnPhidn, dlnPhidT, dlnPhidP, Z, phase_out)

	implicit none
	
	integer, parameter                       :: NCM = 30
	
	integer, intent(in)                      :: n_deriv, NC, model, p_deriv, phase_type, t_deriv
	real(8), intent(in)                      :: n(NCM), P, T
	                                         
	integer, intent(out)                     :: phase_out
	real(8), intent(out)                     :: Z
	real(8), dimension(NCM), intent(out)     :: dlnPhidP, dlnPhidT, lnPhi
	real(8), dimension(NCM,NCM), intent(out) :: n_dlnPhidn
	
	integer                                  :: mixture = 1
	real(8)                                  :: GL, GV, ln_Ntot, Ntot, ZL, ZV
	real(8), dimension(NCM)                  :: dlnPhidPL, dlnPhidPV, dlnPhidTL, dlnPhidTV, lnPhiL, lnPhiV
	real(8), dimension(NCM,NCM)              :: n_dlnPhidnL, n_dlnPhidnV

	GL = 0
	GV = 0
	Ntot = sum(n(:NC))
	ln_Ntot = dlog(Ntot)
	n_dlnPhidn(:NC,:NC) = 0
	dlnPhidT(:NC) = 0
	dlnPhidP(:NC) = 0	
	if (phase_type <= 0) then
	  
! 	  Vapor phase fugacity
	  call vapor_fugacity (model, mixture, NC, t_deriv, p_deriv, n_deriv, T, P, n, lnPhiV, n_dlnPhidnV, dlnPhidTV, dlnPhidPV, ZV)
	  
	endif
	if (phase_type >= 0) then
	  
! 	  Liquid phase fugacity
	  call fugacity_from_activity_model (model, NC, t_deriv, p_deriv, n_deriv, T, P, n, lnPhiL, n_dlnPhidnL, dlnPhidTL, dlnPhidPL, ZL)
	  
	endif
	if (phase_type == 0) then
	  
	  GV = sum (n(:NC) * (dlog(n(:NC)) + lnPhiV - ln_Ntot) )
	  GL = sum (n(:NC) * (dlog(n(:NC)) + lnPhiL - ln_Ntot) )
	  
	endif
	phase_out = 1	  
	if (phase_type > 0 .OR. GL < GV) then
	  
	  lnPhi(:NC) = lnPhiL(:NC)
	  n_dlnPhidn(:NC,:NC) = n_dlnPhidnL(:NC,:NC)
	  dlnPhidT(:NC) = dlnPhidTL(:NC)
	  dlnPhidP(:NC) = dlnPhidPL(:NC)
	  Z = ZL
	  
	else
	  
	  lnPhi(:NC) = lnPhiV(:NC)
	  n_dlnPhidn(:NC,:NC) = n_dlnPhidnV(:NC,:NC)
	  dlnPhidT(:NC) = dlnPhidTV(:NC)
	  dlnPhidP(:NC) = dlnPhidPV(:NC)	  
	  Z = ZV
	  if (phase_type == 0) phase_out = -1
	  
	endif
	return
endsubroutine fugacity_from_gamma_phi
! ---------------------------------------------------------------------------------------------------------------------
!
! Fugacity coeficient predicted by GE model
!
!
subroutine fugacity_from_activity_model (model, NC, t_deriv, p_deriv, n_deriv, T, P, n, lnPhi, n_dlnPhidn, dlnPhidT, dlnPhidP, Z)

	implicit none
	
	integer, parameter                       :: NCM = 30
	real(8), parameter                       :: R = 0.08205D0 !atm m3/(kmol K)
	
	integer, intent(in)                      :: model, n_deriv, NC, p_deriv, t_deriv
	real(8), intent(in)                      :: n(NCM), P, T
	
	real(8), intent(out)                     :: Z
	real(8), dimension(NCM), intent(out)     :: dlnPhidP, dlnPhidT, lnPhi
	real(8), dimension(NCM,NCM), intent(out) :: n_dlnPhidn
	
	integer                                  :: i, NDIF
	real(8), dimension(NCM)                  :: dlnFpuredT, dlnFpuredP, lnFpure, lnGamma, dlnGammadT, vpure, x
	real(8), dimension(NCM,NCM)              :: n_dlnGammadn
	
	integer                                  :: common_NG, common_NGA, common_NST
	common /GCA/                                common_NG, common_NGA, common_NST
    

	if (model == 0) then
	  
	  lnGamma(:NC) = 0
	  n_dlnGammadn(:NC,:NC) = 0
	  dlnGammadT(:NC) = 0
	  
	elseif (model == -1 .OR. model == -2) then

! 	  Package at 16/03/2020 needs to be normalized... 
	  Z = sum(n(:NC))
	  x(:NC) = n(:NC)/Z	  
	  if (n_deriv > 0) NDIF = 1
	  if (t_deriv > 0) NDIF = 2
	  call UNIQUACFAC (NDIF, 0, NC, common_NG, T, x, lnGamma, n_dlnGammadn, dlnGammadT)	
	  if (t_deriv > 0) dlnGammadT(:NC) = dlnGammadT(:NC)/lnGamma(:NC)
	  if (n_deriv > 1) then
	    do i = 1, NC
		
	      n_dlnGammadn(:NC,i) = n_dlnGammadn(:NC,i)/lnGamma(:NC)
		
	    enddo
	  endif
	  lnGamma(:NC) = dlog(lnGamma(:NC))
	  
	endif
	call pureFugacity (model, NC, t_deriv, p_deriv, T, P, lnFpure, dlnFpuredT, dlnFpuredP, vpure)
	lnPhi(:NC) = lnGamma(:NC) + lnFpure(:NC) - dlog(P)
	n_dlnPhidn(:NC,:NC) = n_dlnGammadn(:NC,:NC)
	if (t_deriv > 0) dlnPhidT(:NC) = dlnGammadT(:NC) + dlnFpuredT(:NC)
	if (p_deriv > 0) dlnPhidP(:NC) = dlnFpuredP(:NC) - 1D0/P
	Z = dot_product(n(:NC), vpure(:NC))*P/sum(n(:NC))/R/T
	
	return
endsubroutine fugacity_from_activity_model
! ---------------------------------------------------------------------------------------------------------------------
!
	
subroutine vapor_fugacity (model, mixture, NC, t_deriv, p_deriv, n_deriv, T, P, n, lnPhi, n_dlnPhidn, dlnPhidT, dlnPhidP, Z)
	
	implicit none
	integer, parameter                       :: NCM = 30
	real(8), parameter                       :: R = 82.05D0 !atm cm3/(mol K)
	
	integer, intent(in)                      :: mixture, model, n_deriv, NC, p_deriv, t_deriv
	real(8), intent(in)                      :: n(NCM), P, T
	
	real(8), intent(out)                     :: Z
	real(8), dimension(NCM), intent(out)     :: dlnPhidP, dlnPhidT, lnPhi
	real(8), dimension(NCM,NCM), intent(out) :: n_dlnPhidn
	
	integer                                  :: common_MDIM(NCM,NCM), common_NDIM, i, j       
	real(8)                                  :: aux, Btot, RT, x(NCM)
	real(8), dimension(NCM,NCM)              :: B, BF
	logical                                  :: dimerisation = .FALSE.
	
	COMMON /MPHU/                               common_MDIM, common_NDIM	
	
	if (model == 0) then
	  
	  lnPhi(:NC) = 0
	  dlnPhidP(:NC) = 0
	  dlnPhidT(:NC) = 0
	  n_dlnPhidn(:NC,:NC) = 0	
	  Z = 1
	  
	else
	  
	  if (common_NDIM > 0) dimerisation = .TRUE.
	  x = n(:NC)/sum(n(:NC))
	  RT = R*T
	  call SVIP (NC, mixture, T, B, BF)
	  if (dimerisation) then
	  
	    call PHUB (NC, mixture, P, T, x, B, BF, lnPhi, Z)
! 	    Code for fugacity derivatives in the chem theor. haven't been done yet.
! 	    All gas phase derivatives will be deprecated
	    dlnPhidT(:NC) = 0
	    dlnPhidP(:NC) = 0
	    n_dlnPhidn(:NC,:NC) = 0
	    
	  else
	    
! 	    When subroutine PHUB is updated to cover derivatives, this section should be erased and all 
!	    of this calculations should be done inside PHUB
	    Btot = dot_product(x(:NC),matmul(B(:NC,:NC),x(:NC)))
	    
	    if (Btot > 0 .OR. (P < -RT/4/Btot)) then
	    
	      do i = 1, NC
	        
	        aux = dot_product(x(:NC),B(:NC,i))
	        dlnPhidP(i) = (2*aux - Btot)/RT
	        lnPhi(i) = dlnPhidP(i)*P
	        dlnPhidT(i) = 0
	        if (t_deriv > 0) dlnPhidT(i) = -lnPhi(i)/T !here, we are deprecating completely the temp. dependence of B
	        n_dlnPhidn(:NC,i) = 0
	        if (p_deriv == 0) dlnPhidP(i) = 0
            
	      enddo
	      if (n_deriv > 1) then
	        do i = 1, NC	  
	          do j = 1, NC
	      	
	            n_dlnPhidn(j,i) =  2*B(j,i)*P/RT - lnPhi(j) - lnPhi(i)
	      	
	          enddo
	        enddo
	      endif
	      Z = 1 + Btot*P/RT
	      
	    else
	    
	      lnPhi(:NC) = 0
	    
	    endif
	    
	  endif
	  
	endif
	
	return
ENDsubroutine vapor_fugacity
! ---------------------------------------------------------------------------------------------------------------------
!
subroutine pureFugacity (model, NC, t_deriv, p_deriv, T, P, lnFpure, dlnFpuredT, dlnFpuredP, vpure)
  
	implicit none
	
	integer, parameter                  :: NCM = 30
	real(8), parameter                  :: R = .08205D0, lnPa = dlog(101325D0) !atm m3/(kmol K), Pa->atm
	
	integer, intent(in)                 :: model, NC, p_deriv, t_deriv
	real(8), intent(in)                 :: P, T

	real(8), dimension(NCM), intent(out):: lnFpure, dlnFpuredP, dlnFpuredT, vpure
	
	integer                             :: i, mixture
	real(8)                             :: B2, D2, dlnPOYdT, dVpuredT, lnPOY, RT, aux(NCM,NCM), ZVsat
	real(8), dimension(NCM)             :: dlnPhisatdP = 0, dlnPhisatdT = 0, lnPsat, lnPhisat = 0, n, psat
	
! 	Common variables:
	real(8), dimension(NCM)             :: Avol, Bvol, Cvol, Dvol, Avp, Bvp, Bvp2, Cvp, Dvp, Evp
	
	common /PureParameters/                Avol, Bvol, Cvol, Dvol, Avp, Bvp, Bvp2, Cvp, Dvp, Evp

	RT = R*T
	
	mixture = 0	
	do i = 1, NC
	  
! 	  Saturation volume
	  B2 = (1 - min(.9999999999D0,T/Cvol(i)))
	  D2 = B2**Dvol(i)
	  Vpure(i) = (Bvol(i)**(1 + D2)) /Avol(i)
	  dVpuredT = -D2*dlog(Bvol(i))/Cvol(i)/B2 * Vpure(i)
	  
! 	  Saturation pressure
	  B2 = Bvp(i)/(T + Bvp2(i))
	  D2 = Dvp(i)*T**Evp(i)
	  lnPsat(i) = Avp(i) + B2 + Cvp(i)*dlog(T) + D2 - lnPa
	  psat(i) = dexp(lnPsat(i))
	  
! 	  Fugacity coefficient at saturation point
	  n(:NC) = 0
	  n(i) = 1
	  call vapor_fugacity (model, mixture, NC, t_deriv, p_deriv, 0, T, psat(i), n, lnPhisat, aux, dlnPhisatdT, dlnPhisatdP, ZVsat)
	
! 	  Pointing correction
	  lnPOY = vpure(i)*(P - psat(i))/RT
	  dlnPOYdT = dVpuredT*(P - psat(i))/RT - lnPOY/T
	  
	  lnFpure(i) = lnPsat(i) + lnPhisat(i) + lnPOY
	  if (p_deriv > 0) dlnFpuredP(i) = dlnPhisatdP(i) + vpure(i)/RT
	  if (t_deriv > 0) dlnFpuredT(i) = -B2/(T + Bvp2(i)) + Cvp(i)/T + D2*Evp(i)/T + dlnPhisatdT(i) + dlnPOYdT
	  
	enddo
	
	return
endsubroutine pureFugacity
! ---------------------------------------------------------------------------------------------------------------------
!
subroutine readParametersGammaPhi (model, NC, inputFile, outputFile)
	implicit none

	integer, parameter          :: NCM = 30
	                            
	integer, intent(in)         :: inputFile, model, NC, outputFile
	                           
	real(8)                     :: etaMax, T
	                           
	integer                     :: i, j, NG, unimodel
	real(8), dimension(NCM)     :: Avol, Avp, Bvol, Bvp, Bvp2, Cvol, Cvp, dipole_moment, Dvol, Dvp, Evp, HHA, HHB, HHC, HHD, HHE, &
&	                               HHF, HHG, MW, omega, Pc, Peneloux, radius_gyr, Tc

	real(8), dimension(NCM,NCM) :: eta = 0
	common /PureParameters/        Avol, Bvol, Cvol, Dvol, Avp, Bvp, Bvp2, Cvp, Dvp, Evp
	                            
	                            
	CHARACTER(10)               :: compoundName(NCM)
	COMMON /NAME/                  compoundName
	                           
	integer                     :: common_MDIM(NCM,NCM), common_NG, common_NGA, common_NST, common_NDIM
	common /GCA/                   common_NG, common_NGA, common_NST	

	common /GCPROM/                MW, Peneloux, HHA, HHB, HHC, HHD, HHE, HHF, HHG
	                           
	COMMON /CRIT/                  TC, PC, omega
	COMMON /VIRDAT/                radius_gyr, dipole_moment, ETA
	COMMON /MPHU/                  common_MDIM, common_NDIM
	
	write (outputFile, '("Pure compound properties", /)')
	write (outputFile, '("Compound", 5X, "Tc(K)", 4X, "Pc(atm)", 5X, "omega", 3X, "R`(Å)", 3X, "mu(D)", 2X, "eta")')
	
	common_MDIM = 0; common_NDIM = 0
	do i = 1, NC

	  read (inputFile, *) compoundName(i)
	  read (inputFile, *) Avol(i), Bvol(i), Cvol(i), Dvol(i)
	  read (inputFile, *) Avp(i), Bvp(i), Bvp2(i), Cvp(i), Dvp(i), Evp(i)
	  read (inputFile, *) MW(I), Peneloux(i), HHA(I), HHB(I), HHC(I), HHD(I), HHE(i), HHF(i), HHG(i)	  
	  Tc(I) = Cvol(i)
	  T = Tc(i)
	  Pc(i) = exp(Avp(i) + Bvp(i)/(T + Bvp2(i)) + Cvp(i)*log(T) + Dvp(i)* T**Evp(i))/101325
	  T = .7*T
	  omega(i) = -1 - log10(exp(Avp(i) + Bvp(i)/(T + Bvp2(i)) + Cvp(i)*log(T) + Dvp(i)* T**Evp(i))/101325/Pc(i))
	  if (model == 0) then
	    radius_gyr(I) = 0
	    dipole_moment(I) = 0
	    eta(I,I) = 0
	  else
	    read (inputFile, *) radius_gyr(I), dipole_moment(I), eta(I,I)
	  endif
	  write (outputFile, '(A10, X, F7.2, 4X, G10.3, X, F6.4, 2X, F6.3, 2X, F6.3, 2X, F4.2)') compoundName(i), Tc(i), Pc(i), omega(i), radius_gyr(I), dipole_moment(I), eta(I,I)
        
	enddo
	if (model == -1 .OR. model == -2) then
	  
! 	  -1: UNIFAC -> 0
! 	  -2: UNIQUAC-> 1
	  unimodel = -1 - model
	  call PARINL (NC, NG, unimodel, inputFile, outputFile)
	  common_NG = NG
	
	endif
	
	etaMax = maxval(eta(:NC,:NC))
	if (etaMax > 3.5) then
	  
	  common_NDIM = 1
	  do i = 1, NC
	    do j = i, NC
		
	      if (eta(i,i) > 3.5 .AND. abs(eta(i,i) - eta(j,j)) < 1) then
		  
	        eta(j,i) = (eta(i,i) + eta(j,j))/2
	        eta(i,j) = eta(j,i)
	        common_MDIM(j,i) = 1 !cross association for polar compounds
	        common_MDIM(i,j) = 1
	        
		endif
		  
	    enddo
	  enddo
	  
	endif
	  
	return 
endsubroutine

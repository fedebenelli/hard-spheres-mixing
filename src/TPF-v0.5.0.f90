!  This program aims to calculate a flash at specified vapor fraction (V/F = beta) and pressure
!  or temperature. Iteration variables are T or P (depending the specification) and distribution
!  factors, K.
!  These variables are first scaled logarithmicaly, in order to avoid negative values:
!
!    xi(i) = ln(K(i)), i = 1, NC
!    xi(NC+1) = ln(T) or ln(P).
!
!  The system of NC equations (NC isofugacity conditions + molefraction sumations equal unity) are
!  solved though a Newton-Raphson procedure. Details of the equations may be found in [1].



!  References
!  [1] ML Michelsen, JM Mollerup, "Thermodynamic Models: Fundamentals & Computational Aspects"
!      2nd Ed. Tie-Line Publications, 2007, H�lte, pp. 302:303.

!definition of intel compiler directives
! verbose - output intensity of intermediate results
! 0: No output
! 1: Only initialization status
! 2: Detailed logs for flash calculations
!DIR$ DEFINE verbose = 0


module vaporFraction_flash

	use, intrinsic :: iso_fortran_env, only : ERROR_UNIT, INPUT_UNIT, OUTPUT_UNIT
	
	implicit none

	integer, parameter                   :: inputFile = 35, L = 2, logFile = 33, NCM = 30, outputFile = 31,thermoOut = 32, V = 1

	integer                              :: model, NC, number_heavy_comp, number_perm_gases
	integer, dimension(:), allocatable   :: heavyCompIndex, permGasIndex

	character(len=256)                   :: outputFilename

	logical, dimension(:), allocatable   :: gas, heavy
	logical                              :: GCAinitialized =.false.
	
	
	contains
!
!
!------------------------------------------------------------------------------------------------------------------------
!
!
	subroutine gcaInit(filename, flashDataStr)
  
    implicit none
  
    character(len=*), intent(in)                                           :: filename
    character(len=512), dimension(:), allocatable, optional, intent(inout) :: flashDataStr
                                                                          
    integer, parameter                                                     :: readInc = 100
    integer                                                                :: i, j, readCounter, input_status
    logical                                                                :: EoF = .false.
    character(10)                                                          :: name = "__________"
    character(len=512), dimension(:), allocatable                          :: saveStr1, saveStr2
    character(len=512)                                                     :: locData
  
!------------------------------------------------------------------------------
  
    open (unit = inputFile, file = filename, status = "old")

    read (inputFile, '(G)') name
  
  !DIR$ IF (verbose >= 1)
    write (*, '(/, " Output file = ", A, /)') trim(name)
  !DIR$ ENDIF
    outputFilename = 'output.'//trim(name)//'.txt'

    open (unit = thermoOut, file = 'thermoOut.'//trim(name)//'.txt')
    read (inputFile, *) model, NC
    call ReadParameters (inputFile, thermoOut, model, NC)
    close (unit = thermoOut)
  
    allocate (gas(NC), heavy(NC))
    gas = .FALSE.
    heavy = .FALSE.

  ! This algorithm allows taking into account non-soluble gases and non-volatile compounds.

    read (inputFile, *) number_perm_gases
    if (number_perm_gases > 0) then

      allocate (permGasIndex(number_perm_gases))
      read (inputFile, *) permGasIndex
      do j = 1, number_perm_gases
        do i = 1, NC
          if (i == permGasIndex(j)) then

            gas(i) = .TRUE.
  !           j = j + 1
            exit

          endif
        enddo
      enddo

    endif
    read (inputFile, *) number_heavy_comp
    if (number_heavy_comp > 0) then

      allocate (heavyCompIndex(number_heavy_comp))
      read (inputFile, *) heavyCompIndex
      do j = 1, number_heavy_comp
        do i = 1, NC
          if (i == heavyCompIndex(j)) then

            if (gas(i)) then

              gas(i) = .FALSE.
  !             j = j + 1
              exit

            else

              heavy(i) = .TRUE.
  !             j = j + 1
              exit

            endif

          endif
        enddo
      enddo
    endif

  ! optionally reading input data for flash executions until an error occures:
  ! * Error while reading (iostat <> 0)
    if (present(flashDataStr)) then
	
      allocate (saveStr1(readInc))
      readCounter = 0
      do while (.not. EoF)
	  
        do i = 1, readInc
	    
          read (inputFile, '(A)', iostat = input_status) locData
          if (input_status /= 0 ) then     
		
            EoF=.true.
            exit !exit loop and return.
            
          elseif(len(trim(locData))>0) then
		
            readCounter=readCounter+1
            saveStr1(readCounter)=locData
            
          endif
          
        enddo
        if (.not. EoF) then
	    
          if(allocated(saveStr2)) deallocate(saveStr2)
	
          allocate (saveStr2(readCounter))
          saveStr2(:) = saveStr1(1:readCounter)
          deallocate(saveStr1)
          allocate(saveStr1(readCounter+readInc))
          saveStr1(1:readCounter) = saveStr2
          deallocate (saveStr2)
	    
        end if
        
      enddo
      
      if (allocated(flashDataStr)) deallocate (flashDataStr)
      allocate (flashDataStr(readCounter))
      flashDataStr = saveStr1(1:readCounter)
	
    endif
  
!     close (unit = inputFile)
    GCAinitialized = .true.
  end subroutine gcaInit


  !  -----------------------------------------------------------------------------
  !
  !  This subroutine performs the flash calculation at specified vapor fraction and pressure or
  !  temperature. It uses scaled variables (xi) with for the solution procedure based on Newton-
  !  Raphson method with "step correction".

  !  The simple step correction simply checks 3 things:
  !  * The thermodynamic package has succefully found the asked phase (L or V)
  !  * The variable error is lower thant 2*||F||_2. ||F||_2 is a minimum at solution, but this procedure
  !    doesn't follow a minimization approach. However, such an increase in the error may be doe to an
  !    overjump in some variables different that T or P.
  !  * Maximum jump is delta(xi[NDim]) = 0.5, in order to avoid overjumps in pressure or temperature.
  !

	subroutine Flash (model, NC, calc_type, T, P, feed, trace, gas, heavy, beta, error_thermo, K, w, Z, var_spec, dXi_dS, lnPhi_out)

	implicit none

	integer, parameter                             :: L = 2, maxIt = 50, max_alpha_iter = 10, NCM = 30, V = 1
	real(8), parameter                             :: tolerance = 1.D-8, delta_xi_NDim_max = 2.5D-1
                                                    
                                                    
	integer, intent(in)                            :: model, NC
                                                    
	integer                                        :: alpha_iteration, guess = 0, iteration, i, j, NDim
	integer, dimension(2)                          :: outCond
	integer, dimension(:), allocatable             :: pivot

	real(8), intent(in)                            :: beta
	real(8), intent(inout)                         :: P, T
	real(8), dimension(:), intent(in)              :: feed !dimension(NCM)
	real(8), dimension(2), intent(inout)           :: Z
	real(8), dimension(:), intent(inout)           :: K!dimension(NCM)
	real(8), dimension(:,:), intent(out)           :: w !dimension(NCM,2)
	real(8), dimension(:,:), optional, intent(out) :: lnPhi_out !dimension(NCM,2)
	real(8), dimension(:,:), optional, intent(out) :: dXi_dS !dimension(NC+1,NC+1)

	real(8)                                        :: alpha, alpha_old, d, error, error_new
	real(8), dimension(:), allocatable             :: delta_xi, delta_xi_Old, F, xi, xi_new
	real(8), dimension(:,:), allocatable           :: jacobian, lnPhi
	real(8), dimension(:,:,:), allocatable         :: n_dLnPhi_dn
                                                    
                                                    
	character, intent(in)                          :: calc_type
                                                    
	logical, dimension(NC), intent(in)             :: gas, heavy, trace
	logical, dimension(:), optional, intent(in)    :: var_spec
	logical, dimension(:), allocatable             :: local_var_spec
	logical, intent(out)                           :: error_thermo


	guess = 2
	alpha_old = 1
	error_thermo=.false.
!	Initialization of K:
!	call initialization (model, NC, calc_type, T, P, feed, trace, gas, heavy, beta, K, w, Z, error_init)

!	Size of the problem
	NDim = NC + 1
	j = NC
	do i = 1, NC

	  if (gas(j) .OR. heavy(j)) then

!         Reduce the size of the problem for gas or heavy comp. introduced at the end, since
!         we aren't interested in their isofugacity.
	    NDim = NDim - 1
	    if (gas(j)) then

		K(j) = 1D100

	    elseif (heavy(j)) then

		K(j) = 1D-100

	    endif

	  else
!         The first comp. that distributes between V and L phases appeared. The size can't be reduced any more in this way.
	    exit

	  endif
	  j = j - 1

	enddo
	allocate (xi(NDim), xi_new(NDim), delta_xi(NDim), delta_xi_old(NDim), F(NDim), jacobian(NDim,NDim+1), pivot(NDim),  &
	          n_dLnPhi_dn(NCM,NCM,2), lnPhi(NCM,2))
	delta_xi = 0
	do i = 1, NDim - 1

	  if (gas(i)) then

	    K(i) = 1D100
	    xi(i) = dlog(K(i))

	  elseif (heavy(i)) then

	    K(i) = 1D-100
	    xi(i) = dlog(K(i))

	  else

	    xi(i) = dlog(K(i))

	  endif

	enddo
	if (calc_type == 'P') then

	  xi(NDim) = dlog(T)

	else

	  xi(NDim) = dlog(P)

	endif
!  	First error evaluation before start convergence loop
	call Func (model, NC, NDim, calc_type, T, P, feed(:NCM), trace, gas, heavy, beta, K, xi,&
	           error_thermo, w, Z, guess, F, jacobian, lnPhi, n_dLnPhi_dn, error)    
	if(error_thermo .OR. error > 1.0d0) then
	
! 	  Recalculate with an EOS intrinsic volume guess 
	  guess = 0
	  call Func (model, NC, NDim, calc_type, T, P, feed(:NCM), trace, gas, heavy, beta, K, xi,&
	             error_thermo, w, Z, guess, F, jacobian, lnPhi, n_dLnPhi_dn, error)
	  
	else
	  
! 	  Use previous results of compressibility factors as initial values for next iteration
	  guess = 2
	  
	end if
	!DIR$ IF(verbose>=2)
	call writeLog (0, NDim, NC, 0, 1D0, T, P, K(:NC), F, error)
	!DIR$ ENDIF
!
!
!	Convergence loop:
!
	do iteration = 1, maxIt
!
	  delta_xi_old = delta_xi
	  
!	  Solve Newton-Raphson equations: J�DXi = -f
	  call LUDcmp (jacobian(:NDim,:NDim), NDim, NDim, pivot, d)
	  if (error < tolerance) then

	    exit

	  endif
	  delta_xi = -F
	  call LUBksb (jacobian(:NDim,:NDim), NDim, NDim, pivot, delta_xi)

!  	  Cut step to avoid big jumps in pressure or temperature (which are the variables fed to
!  	  the thermo package):
	  if ( dabs(delta_xi(NDim)) < delta_xi_NDim_max ) then

	    alpha = 1.D0

	  else

	    alpha = delta_xi_NDim_max/dabs(delta_xi(NDim))

	  endif
!	  Trying to prevent oscilations...
	  d = norm2(delta_xi*alpha + delta_xi_old*alpha_old, NDim) !Cancelation of 2 consecutive steps.
	  d = d/error
	  if (d < 1d-1) then

	    alpha = alpha/2

	  endif
! 	  Step may be cutted even more if there's a problem with the convergence, which is performed
! 	  within this loop up to correct evaluation or
	  alpha_iteration = 0
	  do

	    alpha_iteration = alpha_iteration + 1
	    do i = 1, NDim - 1

	      if ((.NOT. gas(i)) .AND. (.NOT. heavy(i))) then
            
	        if (trace(i)) then
            
	          xi_new(i) = lnPhi(i, L) - lnPhi(i, V)
            
	        else
            
	          xi_new(i) = xi(i) + delta_xi(i)*alpha
            
	        endif
	        K(i) = dexp(xi_new(i))

	      endif

	    enddo
	    xi_new(NDim) = xi(NDim) + delta_xi(NDim)*alpha
	    if (calc_type == 'T') then

	      d = P  !store old pressure to update Z(L) if T is specified
	      P = dexp(xi_new(NDim))

	    else

	      T = dexp(xi_new(NDim))

	    endif
! 	    if ((iteration > 5 .OR. error < 1) .AND. (.NOT. error_thermo)) then
! 	      guess = 2 ! Using previous values of Z as initial values for next step, but only
	                ! after 5 iterations or the error is already low.
	    if ((.NOT. error_thermo) .AND. guess == 2) then		
	      if (calc_type == 'T') Z(L) = Z(L) * P / d
	    endif

	    call Func (model, NC, NDim, calc_type, T, P, feed, trace, gas, heavy, beta, K, xi_new, error_thermo, w, Z, guess, F, jacobian, lnPhi, n_dLnPhi_dn, error_new)
	   !DIR$ IF(verbose>=2)
	    call writeLog (iteration, NDim, NC, alpha_iteration, alpha, T, P, K(:NC), F, error_new)
	   !DIR$ ENDIF
	    if (.NOT. error_thermo .AND. error_new < 10*error) then

!	      No problems have been found
	      xi = xi_new
	      error = error_new
	      alpha_old = alpha
	      exit

	    else

	      if (alpha_iteration > max_alpha_iter) then

	        error_thermo = .TRUE.
	        return ! I was not possible to find a way to calculate fugacity coefficients

	      endif
!	      Cut step in this situation. The decrease factor reduces one order of magnitude delta_xi
!  	      after two alpha-steps
	      alpha = alpha/3.1622777D0
!	      Tell the thermo package to initialize ZL and ZV only if an problem has occurred.
	      if (error_thermo) guess = 0

	    endif

	  enddo

	enddo
	
	if (error > tolerance) then
	  
	  error_thermo = .TRUE.
	  
	end if
	
	if (present(dXi_dS)) then
!	  Calculation of the derivatives of the calculated variables wrt the specified ones.
!	  call sensitivity_vector (NDim, NC, jacobian(:NDim,:NDim), pivot, jacobian(:NDim,NDim+1), feed(:NC), w(:NC,L), w(:NC,V), &
!                         n_dLnPhi_dn(:NC,:NC,L), n_dLnPhi_dn(:NC,:NC,V), gas(:NC), heavy(:NC), dXi_dS)
	  allocate (local_var_spec(NC+2))
	  if ( present(var_spec)) then
	    
	    local_var_spec = var_spec
	    
	  else
	    
	    local_var_spec = .TRUE. !if var_spec hasn't been specified, we assume the user wants *everithing*
	    
	  endif
	  call sensitivity_vector2(NDim, NC, jacobian(:NDim,:NDim), pivot, jacobian(:NDim,NDim+1), beta, K(:NC), &
	                           w(:NC,L), w(:NC,V), n_dLnPhi_dn(:NC,:NC,L), n_dLnPhi_dn(:NC,:NC,V), gas(:NC), &
	                           heavy(:NC), local_var_spec, dXi_dS)

	end if
	if (present(lnPhi_out)) then
	  
	  lnPhi_out(1:NC,:) = lnPhi(1:NC,:)
	  
	end if
	

	return
	endsubroutine
!  -----------------------------------------------------------------------------
!
!  Evaluation of the isofugacity condition and mole fraction sums, together with the jacobian,
!  molefractions and other variables.
!
  
	subroutine Func (model, NC, NDim, calc_type, T, P, feed, trace, gas, heavy, beta, K, xi, error_thermo, w, Z, guess, F,    &
	                 jacobian, lnPhi, n_dLnPhi_dn, error)

	implicit none

	integer, parameter              :: calc = 6, NCM = 30

	integer                         :: guess, i, j, L = 2, model, NC, NDim, V = 1
	integer, dimension(2)           :: outCond

	real(8)                         :: B, beta, Cp, error, H, P, T

	real(8), dimension(2)           :: factor, Z
	real(8), dimension(:)           :: feed !dimension(NCM)
	real(8), dimension(:)           :: K !dimension(NCM)
	real(8), dimension(:,:)         :: lnPhi, w !dimension(NCM,2)
	real(8), dimension(:,:,:)       :: n_dLnPhi_dn !dimension(NCM,NCM,2)
	real(8), dimension(:,:)         :: jacobian !dimension(NDim,NDim+1)
	
	real(8), dimension(NCM,2)       :: dH_dn, dLnPhi_dP, dLnPhi_dT, w_loc

	real(8), dimension(NC)          :: denom, dxdxi, dydxi

	real(8), dimension(NDim)        :: F, xi

	character                       :: calc_type
	logical                         :: error_thermo
	logical, dimension(NC)          :: gas, heavy, trace
	
	w_loc=0.0d0
	w_loc(:NC,:)=w(:NC,:)

	call molar_fractions (NC, trace, gas, heavy, beta, K(:NC), feed(:NC), w_loc(:NC, :2), denom)
	call THERMO(model, NC, L-1, calc, guess, T, P, w_loc(:NCM, L), Z(L), lnPhi(:NCM, L), n_DLnPhi_dN(:NCM, :NCM, L), DLnPhi_dT(:NCM, L), DLnPhi_dP(:NCM, L), H, DH_dN(:NCM, L), Cp, B, outCond(L))
	call THERMO(model, NC, V-2, calc, guess, T, P, w_loc(:NCM, V), Z(V), lnPhi(:NCM, V), n_DLnPhi_dN(:NCM, :NCM, V), DLnPhi_dT(:NCM, V), DLnPhi_dP(:NCM, V), H, DH_dN(:NCM, V), Cp, B, outCond(V))
	
	w(:NC,:)=w_loc(:NC,:)
    
	IF (outCond(L) /= 1 .AND. outCond(V) /= 1 .OR. abs(Z(V) - Z(L)) < 1.E-6) then

!  	  Error during calculation of any phase, generally due to supercritical conditions.
	  error_thermo = .true.
	  error=1.0d16
	  return

	else

	  error_thermo = .false.

	endif
	do j = 1, NDim - 1

	  if ((.NOT. gas(j)) .AND. (.NOT. heavy(j))) then

	    dxdxi(j) = -w(j, V)*beta/denom(j)
	    dydxi(j) = w(j, V) + dxdxi(j)*K(j)
	    do i = 1, NDim - 1

	      if ((.NOT. gas(i)) .AND. (.NOT. heavy(i))) then

!             Jacobian based on number of moles derivatives (L*x(i) and V*y(i))
	        jacobian(i, j) = n_dLnPhi_dn(i, j, V)*dydxi(j) - n_dLnPhi_dn(i, j, L)*dxdxi(j)
          
	      else
          
	        jacobian(i, j) = 0.D0
          
	      endif

	    enddo
	    jacobian(j, j) = jacobian(j, j) + 1.D0

!	    F(i) is evaluated only if the compound is present
	    if (trace(j)) then

	      F(j) = 0.0D0

	    else

	      F(j) = xi(j) + lnPhi(j, V) - lnPhi(j, L) !isofugacity

	    endif
	    if (calc_type == 'T') then

	      jacobian(j, NDim) = (dLnPhi_dP(j, V) - dLnPhi_dP(j, L))*P    !used in the system of eq. to solve present point
	      jacobian(j, NDim+1) = (dLnPhi_dT(j, V) - dLnPhi_dT(j, L))*T  !used only for the sensitivity vector

	    else

	      jacobian(j, NDim) = (dLnPhi_dT(j, V) - dLnPhi_dT(j, L))*T    !used in the system of eq. to solve present point
	      jacobian(j, NDim+1) = (dLnPhi_dP(j, V) - dLnPhi_dP(j, L))*P  !used only for the sensitivity vector

	    endif
	    jacobian(NDim, j) = dydxi(j) - dxdxi(j)

	  else

    !  All elements of F and J refering to this component are null
	    F(j) = 0.D0
	    do i = 1, NDim

	      if (i == j) then
            
	        jacobian(i, j) = 1.D0
            
	      else
            
	        jacobian(i, j) = 0.D0
            
	      endif

	    enddo
	    jacobian(j, NDim) = 0.D0
	    jacobian(j, NDim+1) = 0.D0

	  endif

	enddo
	F(NDim) = sum( w(:NC, V) - w(:NC, L) ) ! mole fraction sums
	jacobian(NDim, NDim) = 0.D0
	jacobian(NDim, NDim+1) = 0.D0
	error = norm2(F, NDim)

	return
  endsubroutine Func
! -------------------------------------------------------------------------------------------------------
!
! Initialization of K searching for both, L and V at P, T and feed. Only if both phases shows errors
! the program is stoped.
!
	subroutine initialization (model, NC, calc_type, T, P, feed, trace, gas, heavy, beta, K, w, Z, error_init)

	implicit none

	integer, parameter            :: guess = 0, L = 2, NCM = 30, V = 1

	integer                       :: calc = 1, i, iteration, model, NC, outCondL, outCondV

	real(8)                       :: B, beta, Cp, F, H, P, step, SumL, SumV, T

	real(8), dimension(2)         :: Z
	real(8), dimension(NC,2)      :: w
	real(8), dimension(NCM,2)     :: w_loc

	real(8), dimension(NC)        :: denom, K, omega, Pc, Tc, feed
	real(8), dimension(NCM)       :: dH_n, dLnPhi_P = 0, dLnPhi_T, lnPhiL, lnPhiV
	real(8), dimension(NCM,NCM,2) :: dLnPhi_n


	character                     :: calc_type
	logical                       :: error_init
	logical, dimension(NC)        :: gas, heavy, trace

	call get_critical_properties (NC, Tc, Pc, omega)
	error_init = .FALSE.
	do i = 1, NC

	  if (gas(i)) then

	    k(i) = 6.02E23

	  elseif (heavy(i)) then

	    K(i) = 0

	  else

	    K(i) = exp(-log(P/Pc(i)) + 5.373D0*(1 + omega(i))*(1 - Tc(i)/T))

	  endif

	enddo
!	
!	Direct substitution loop
	do iteration = 1, 4
	  
	  call molar_fractions (NC, trace, gas, heavy, beta, K(:NC), feed(:NC), w(:NC, :2), denom)
	  SumV = dlog(sum(w(:NC,V)))
	  SumL = -dlog(sum(w(:NC,L)))
	  if (calc_type == "T") then

	    P = dexp( dlog(P) + (SUMV + SumL)/2 )

	  elseif (calc_type == "P") then

	    step = min(dabs(dlog(T)/40), dabs((SUMV + SumL)/100)) * (-SUMV - SumL)/dabs(SUMV + SumL)
	    T = dexp( dlog(T) + step )

	  endif

	  if (T < maxval(Tc(:NC))) then
	    
	    w_loc = 0.0d0
	    w_loc(1:NC,:) = w
	    call THERMO(model, NC, L-1, calc, guess, T, P, w_loc(:NCM,L), Z(L), lnPhiL, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, outCondL)
	    call THERMO(model, NC, V-2, calc, guess, T, P, w_loc(:NCM,V), Z(V), lnPhiV, DLnPhi_N, DLnPhi_T, DLnPhi_P, H, DH_N, Cp, B, outCondV)
	    w = w_loc(1:NC,:)
	    
	  endif
	  if (outCondL < 0 .AND. outCondV < 0) then

	    error_init = .true.
	    return

	  elseif (outCondL < 0 .AND. outCondV > 0) then

	    lnPhiL(:NC) = 0

	  elseif (outCondL > 0 .AND. outCondV < 0) then

	    lnPhiV(:NC) = 0

	  endif
	  K(:NC) = exp(lnPhiL(:NC) - lnPhiV(:NC))

	enddo
	return
	endsubroutine initialization
! -------------------------------------------------------------------------------------------------------
!
!
  subroutine guessPressure(NC,T,P,feed)
  
    implicit none
    
    integer, intent(in)       :: NC
    real(kind=8), intent(in)  :: T
    real(kind=8), intent(in)  :: feed(:)
    real(kind=8), intent(out) :: P
    
    integer :: ii
    real(kind=8) :: fact
    real(kind=8), dimension(NC) :: Tc, Pc, omega
    
    call get_critical_properties(NC, Tc, Pc, omega)
    
    fact=7.0d0/(3.0d0*T)
    P=0.0d0
    do ii=1,NC
      P=P+feed(ii)*Pc(ii)*dexp(fact*(omega(ii)+1.0d0)*(T-Tc(ii)))
    end do
      
  end subroutine guessPressure
  
  subroutine molar_fractions (NC, trace, gas, heavy, beta, K, z, w, denom)

    implicit none

    integer, parameter            :: L = 2, V = 1
    integer                       :: i, NC

    real(8)                       :: beta
    real(8), dimension(NC)        :: denom, K, z
    real(8), dimension(NC,2)      :: w

    logical, dimension(NC)        :: gas, heavy, trace


    do i = 1, NC

      if (.NOT. trace(i)) then
        if (gas(i)) then

          w(i, V) = z(i)/beta
          w(i, L) = 0
          denom(i) = 6.02E23

        elseif (heavy(i)) then

          w(i, V) = 0
          denom(i) = 1.D0 - beta
          w(i, L) = z(i)/denom(i)

        else

          denom(i) = 1.D0 - beta + beta*K(i)
          w(i, L) = z(i)/denom(i)
          w(i, V) = K(i)*w(i, L)

        endif
      else

        denom(i) = 1
        w(i, L) = z(i)
        w(i, V) = z(i)*K(i)

      endif

    enddo

    return
  endsubroutine molar_fractions
!---------------------------------------------------------------------------
!
!
!
!
! 	subroutine sensitivity_vector (NDim, NC, jacobian, pivot, dF_dS, z, x, y, L_dlnPhi_dl, V_dlnPhi_dv, &
! 	                               gas, heavy, dXi_dS)
! 
! 	implicit none
! 
! 	integer                                   :: i, j
! 
! 	integer, intent(in)                       :: NDim, NC
! 	integer, dimension(NDim), intent(in)      :: pivot
! 	real(8), dimension(NC), intent(in)        :: x, y, z
! 	real(8), dimension(NC,NC), intent(in)     :: L_dlnPhi_dl, V_dlnPhi_dv
! 	real(8), dimension(NDim), intent(in)      :: dF_dS
! 	real(8), dimension(NDim,NDim), intent(in) :: jacobian
! 
! 	logical, dimension(NC), intent(in)        :: gas, heavy
! 
! 	real(8), dimension(NC+1,NC+1), intent(out):: dXi_dS
! 	
! 	dXi_dS=0.0d0
! 
! 	do j = 1, NC
! 
! 	  do i = 1, NC
! 	    if ((.NOT. gas(i)) .AND. (.NOT. heavy(i))) then
! 
! 	      dxi_dS(i, j) = -(V_dlnPhi_dv(i, j)*y(j) - L_dlnPhi_dl(i, j)*x(j)) !/z(j)
! 
! 	    else
! 
! 	      dxi_dS(i, j) = 0.0D0
! 
! 	    endif
! 	  enddo
! 	  dxi_dS(NDim, j) = -(y(j) - x(j)) !/z(j)
! 	  call LUBksb (jacobian, NDim, NDim, pivot, dXi_dS(:NDim, j))
! 	  if (NDim < NC+1) then
! 	    dxi_dS(NC+1, j) = dxi_dS(NDim, j)
! 	    dxi_dS(NDim, j) = 0D0
! 	  endif
! 
! 	enddo
! 	dXi_dS(:NDim, NC+1) = -dF_dS
! 	call LUBksb (jacobian, NDim, NDim, pivot, dXi_dS(:NDim, NC+1))
! 	if (NDim < NC+1) then
! 	  dxi_dS(NC+1, NC+1) = dxi_dS(NDim, NC+1)
! 	  dxi_dS(NDim, NC+1) = 0D0
! 	endif
! 	return
! 	endsubroutine sensitivity_vector
!---------------------------------------------------------------------------
!
!
!
!
	subroutine sensitivity_vector2(NDim, NC, jacobian, pivot, dF_dS, beta, K, x, y, L_dlnPhi_dl, V_dlnPhi_dv, &
	                               gas, heavy, var_spec, dXi_dz)

	implicit none

	integer                                   :: i, j

	integer, intent(in)                       :: NDim, NC
	integer, dimension(NDim), intent(in)      :: pivot
	real(8), intent(in)                       :: beta
	real(8), dimension(NC), intent(in)        :: K, x, y
	real(8), dimension(NC,NC), intent(in)     :: L_dlnPhi_dl, V_dlnPhi_dv
	real(8), dimension(NDim), intent(in)      :: dF_dS
	real(8), dimension(NDim,NDim), intent(in) :: jacobian

	logical, dimension(NC), intent(in)        :: gas, heavy
	logical, dimension(NC+2), intent(in)      :: var_spec

	real(8), dimension(NC+1,NC+2), intent(out):: dXi_dz
	
	real(8), dimension(NC)                    :: denom, dx_dbeta, dy_dbeta
	
	dXi_dz = 0.D0
	
	do i = 1, NC
	  
	  denom(i) = 1.D0/(1.D0 - beta + beta*K(i))
	  if (var_spec(NC+2)) then
	    
	    dx_dbeta(i) = - (y(i) - x(i))*denom(i)
	    dy_dbeta(i) = K(i)*dx_dbeta(i)
	    
	  endif	
	  
	end do
    
	do j = 1, NC

	  do i = 1, NC
	    if ((.NOT. gas(i)) .AND. (.NOT. heavy(i))) then

	      if (var_spec(i)) then
	  
	        dxi_dz(i, j) = -(V_dlnPhi_dv(i, j)*K(j) - L_dlnPhi_dl(i, j)) *denom(j)
	  
	      endif
            if (var_spec(NC+2)) then		
	
! 	        ntot*d[lnPhi(i)]/dn(j) is a symmetrical matrix:
	        dxi_dz(j, NC+2) = dxi_dz(j, NC+2) - V_dlnPhi_dv(i, j)*dy_dbeta(i) + L_dlnPhi_dl(i, j)*dx_dbeta(i)
	  
            endif
            
	    else
		
	      dxi_dz(i, j) = 0.0D0

	    endif
	  enddo
	  if (var_spec(j)) then
	  
	    dxi_dz(NDim, j) = -(K(j) - 1.0d0)*denom(j)
	    call LUBksb (jacobian, NDim, NDim, pivot, dXi_dz(:NDim, j))
	    if (NDim < NC+1) then
	
	      dXi_dz(NC+1, j) = dXi_dz(NDim, j)
	      dXi_dz(NDim, j) = 0D0
	
	    endif
	  
	  endif

	enddo
	if (var_spec(NC+1)) then
	  
	  dXi_dz(:NDim, NC+1) = -dF_dS
	  call LUBksb (jacobian, NDim, NDim, pivot, dXi_dz(:NDim, NC+1))
	  if (NDim < NC+1) then
	    
	    dXi_dz(NC+1, NC+1) = dXi_dz(NDim, NC+1)
	    dXi_dz(NDim, NC+1) = 0D0
	    
	  endif
	  
	endif
	if (var_spec(NC+2)) then

!       Evaluate sensitivity with respect to beta
	  dxi_dz(NDim, NC+2) = -sum(dy_dbeta(:NC) - dx_dbeta(:NC))
	  call LUBksb (jacobian, NDim, NDim, pivot, dXi_dz(:NDim, NC+2))
	  if (NDim < NC+1) then
	    
	    dxi_dz(NC+1, NC+2) = dxi_dz(NDim, NC+2)
	    dxi_dz(NDim, NC+2) = 0D0
	    
	  endif

	endif	
	return
	
	endsubroutine sensitivity_vector2                             
                                 
  ! -------------------------------------------------------------------------------------------------------
  !
  !
  !
  subroutine gcaWriteResuls (model, NC, calc_type, outputFile, error_thermo, T, P, feed, trace, gas, heavy, beta, K, w, Z)

    implicit none

    integer, parameter            :: inputFile = 5, thermoOut = 2, L = 2, NCM = 30, V = 1

    integer                       :: i, input_status = 0, j, model, NC, outputFile

    real(8)                       :: beta, Ntot, P, T

    real(8), dimension(2)         :: Z
    real(8), dimension(NCM)       :: feed, K, n
    real(8), dimension(NCM,2)     :: w

    character                     :: calc_type
    logical                       :: error_init, error_thermo
    logical, dimension(NC)        :: gas, heavy, trace


    if (error_thermo) then

      write (outputFile, '(2X, F7.2, 2X, G10.4, 2X, F5.3, " No VLE found at this conditions", ///)') T, P, beta
      error_thermo = .false. !continue to next point

    else

      do i = 1, NC
        if (trace(i)) then

  !  Compounds in small amounts are set to zero, since they don't really do anything to the mixture. The
  !  magnitude, however, is arbitrary and may be changed.
          w(i,:2) = 0

        endif
        if (gas(i)) then

          K(i) = 1D100

        elseif (heavy(i)) then

          K(i) = 1D-100

        endif
      enddo
      write (outputFile, '(2X, F7.2, 2X, G10.4, 2X, G10.4, "  y = [ ", <NC>(G14.8, X), " ]", F10.7, /, &
  &                                               3(" -"), 28X, " x = [ ", <NC>(G14.8, X), " ]", F10.7, /, &
  &                                           3(" -"), 24X, " ln(K) = [ ", <NC>(G14.8, X), " ]", / )') &
  &                                                      T, P, beta, (w(:NC, i), Z(i), i = 1, 2), log(K(:NC))
      write (151, '(2X, F7.2, 2X, G12.6, 2X, F5.3, "  y = [ ", <NC>(G10.4, X), " ]", F10.7,  &
  &                                                  X, " x = [ ", <NC>(G10.4, X), " ]", F10.7,  &
  &                                              X, " ln(K) = [ ", <NC>(G10.4, X), " ]" )') &
                                                        T, P, beta, (w(:NC, i), Z(i), i = 1, 2), log(K(:NC))
    endif

    return

  end subroutine gcaWriteResuls
  !---------------------------------------------------------------------------
  !
  !
  subroutine writeLog (iteration, NDim, NC, alpha_iteration, alpha, T, P, K, F, error)
    implicit none
    integer                       :: alpha_iteration, iteration, NC, NDim

    real(8)                       :: alpha, error, P, T
    real(8), dimension(NC)        :: K
    real(8), dimension(NC+1)      :: F_tot
    real(8), dimension(NDim)      :: F

    logical, dimension(NC)        :: gas, heavy

    F_tot = 0
    F_tot(:NDim-1) = F(:NDim-1)
    F_tot(NC+1) = F(NDim)
    write (logFile, '(X, 2(I2, X), G11.4, F7.2, X, <2*NC+3>(G11.4, X) )') iteration, alpha_iteration, alpha, T, P, dlog(K(:NC)), F_tot(:NC+1), error

    return
  endsubroutine writeLog
  !---------------------------------------------------------------------------
  !
  !
  !  Norm of a vector v. n denotes dimension and typ (char*1) types: 1 to 9 and infinite (i).
  !
  double precision function GCAnorm (v, n, typ)

    integer n
    real(8) v(n), p
    character(1) typ
    integer ip

    if(typ == 'i')then
  !
  !  Norm infinity
      GCAnorm = maxval(v)

    else
  !
  !  Any other norm, with power between 1 and 9. Of course, dabs(v(i))**ip = |v(i)| if ip = 1.
      read(typ,'(I1)',ERR=100,END=100) ip

      GCAnorm = sum(dabs(v(:N))**ip)
      p = 1.D0/dfloat(ip)
      GCAnorm = GCAnorm**p

    endif

    return
  !----------------------------------------------------------
  100 write(*,'(//,"*** ERROR ***",//,"Function NORM requires an integer from 1 to 9 or i in typ variable to work",//)')
    stop

  endfunction GCAnorm

end module vaporFraction_flash
  
!  
!program fractionFlash
!
!  use m_gca
!
!  implicit none
!
!  integer                              :: i, input_status = 0, j!, strlen
!
!  real(8)                              :: beta, Ntot, P, T
!
!  real(8), dimension(2)                :: Z
!  real(8), dimension(NCM)              :: feed, K, n
!  real(8), dimension(NCM,2)            :: lnPhi, w
!  real(8), dimension(:,:), allocatable :: dXi_dS
!
!  character                            :: calc_type
!  character(len=256), dimension(:), allocatable :: flashInput
!  
!  character(10), dimension(NCM)        :: CompoundName
!  common /NAME/                           CompoundName
!
!  logical                              :: error_init, error_thermo = .false.
!  logical, dimension(:), allocatable   :: trace  
!  
!  
!  !initialize GCA
!  call gcaInit('input.txt',flashInput)
!  
!  allocate (trace(NC), dXi_dS(NC+1,NC+1))
!
!
!! Succesive flashes will be performed, at user specified conditions until an error occure:
!! * Error while reading (iostat <> 0)
!! * T,P = 0
!! * Sum (n(i)) = 0)
! 
!!DIR$ IF(verbose>=2)
!  open (unit = logFile, file = 'log.txt')
!!DIR$ ENDIF
!  open (unit = outputFile, file = outputFilename)
!  write (outputFile, '(3X, "T(K)    P(atm)      V/F", 3X, 3(" -"), 6X, <NC>(A10, 5X), 5X, "Z")') CompoundName(:NC)
!
!  j=1
!  do
!    read (flashInput(j), *, iostat = input_status) (n(i), i = 1, NC), T, P, beta, calc_type
!    Ntot = sum(n(:NC))
!    if (input_status /= 0 .OR. P <= 0 .OR. T <= 0 .OR. Ntot <= 0) then
!
!      exit !exit main loop and then stop.
!
!    endif
!    j=j+1
!
!! Input mole fractions:
!
!    feed(:NC) = n(:NC)/sum(n(:NC))
!    trace(:NC) = .FALSE.
!    do i = 1, NC
!
!      !zero amounts not allowed... perhaps thermo subroutines should be able to
!      !bypass zero amount components unless the user ask for "infinite dilution values".
!      !Here, only are omited in the evaluation of their isofugacity.
!      if (feed(i) <= 0D0) then
!
!        feed(i) = sum(n(:NC))/6.02e23   !1 molecule per mole.
!        trace(i) = .TRUE.
!
!      endif
!
!    enddo
!    if (beta < 0) then
!
!      beta = 0
!
!    elseif (beta > 1) then
!
!      beta = 1
!
!    endif
!    
!!DIR$ IF(verbose>=2)
!    write (logFile, *)
!!DIR$ ENDIF
!    
!    call initialization (model, NC, calc_type, T, P, feed, trace, gas, heavy, beta, K, w, Z, error_init)
!    if(error_init) then
!      write(error_unit,'(A,I0,A)') "Error while initializing flash. Skipping input ",j,"."
!      CYCLE
!    end if
!      
!    call gcaFlash (model, NC, calc_type, T, P, feed, trace, gas, heavy, beta, error_thermo, K, w, Z, dXi_dS, lnPhi)
!    call gcaWriteResuls (model, NC, calc_type, outputFile, error_thermo, T, P, feed, trace, gas, heavy, beta, K, w, Z)
!    do i = 1, NC+1
!      write (outputFile, '(X, <NC+1>(G11.4, X))') dxi_dS(i,:NC+1)
!    enddo
!    write (outputFile, *)
!    do i = 1, 2
!      write (outputFile, '(" - - - ln(�) = [", 7X, <NC+1>(G11.4, X))') lnPhi(:NC, i)
!    enddo
!    write (outputFile, *)
!  enddo
!  
!!DIR$ IF(verbose>=2)
!  close (unit = logFile)
!!DIR$ ENDIF
!  close (unit = outputFile)
!  stop
!
!end program

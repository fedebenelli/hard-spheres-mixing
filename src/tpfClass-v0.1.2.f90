module tpfClass
  
	use, intrinsic :: iso_fortran_env, only : ERROR_UNIT, INPUT_UNIT, OUTPUT_UNIT
	use vaporFraction_flash
	
	implicit none
	
	real(kind = 8), parameter :: zeroTol = 1.0d-16
	
	type t_flashData
	
	  logical, private                                   :: initialized = .false.
	  logical, private                                   :: calcStatus = .false.
	  integer                                            :: NC, model
	  character                                          :: calc_type
	  real(kind = 8)                                     :: T, P, beta, Z(2)
	  real(kind = 8), allocatable, dimension(:)          :: feed, K
	  real(kind = 8), allocatable, dimension(:), private :: Spec, xi !Spec and xi are used for extrapolation and only updated by a full flash calculation
	  real(kind = 8), allocatable, dimension(:,:)        :: w, lnPhi, dXi_dz
	  logical, allocatable, dimension(:)                 :: trace, var_spec
	  
	  contains
	    procedure :: init => flashData_init
	    procedure :: calc => flashData_calc
	    procedure :: extrapolate => flashData_extOpt
! 	    procedure :: extrapolateExpLog => flashData_extExpLog
! 	    procedure :: extrapolateExpLin => flashData_extExpLin
! 	    procedure :: extrapolateLinLin => flashData_extLinLin
	    procedure :: print => flashData_print
	    procedure :: trace_definition => trace_compounds
	    final     :: flashData_delete
	    
	end type t_flashData
  
  !constructor
!   interface t_flashData
!     procedure :: flashData_create
    !declaration of multiple functions is possible here, if using Fortran 2008 standar (IFORT > 12)
!   end interface
  
!============================================================================== 
  contains
!==============================================================================
  
	function flashData_create(NC) result(this)
	  type(t_flashData) :: this
	  integer :: NC
	  
	  allocate(this%feed(NC), this%K(NC), this%Spec(NC+2), this%xi(NC+1), this%w(NC,2),&
	           this%lnPhi(NC,2), this%dXi_dz(NC+1,NC+2),                                &
	           this%trace(NC), this%var_spec(NC+2))
	end function
	
	subroutine flashData_delete(flashData)
	  
	  type(t_flashData) :: flashData
	  
	  if (allocated(flashData%feed)) deallocate(flashData%feed)
	  if (allocated(flashData%K)) deallocate(flashData%K)
	  if (allocated(flashData%Spec)) deallocate(flashData%Spec)
	  if (allocated(flashData%xi)) deallocate(flashData%xi)
	  if (allocated(flashData%w)) deallocate(flashData%w)
	  if (allocated(flashData%trace)) deallocate(flashData%trace)
	  if (allocated(flashData%var_spec)) deallocate(flashData%var_spec)
	  
	  flashData%initialized = .false.   
	  
	end subroutine flashData_delete
	  
	subroutine flashData_init(flashData, model, NC, calc_type, T, P, feed, trace, beta, error_init)
    
    implicit none

    class(t_flashData)                  :: flashData
    integer, intent(in)                 :: model, NC
    real(8), intent(in)                 :: beta
    real(8), intent(inout)              :: P, T
    real(8), dimension(:), intent(in)   :: feed    
    character, intent(in)               :: calc_type    
    logical, intent(out)                :: error_init
    logical, dimension(NC), intent(in)  :: trace  
    
    real(kind = 8) :: w(NC,2)
    
    if (.not.GCAinitialized) then
      write(ERROR_UNIT, *) "GCA Class: GCA must be initialized prior to the initialization of flash objects"
      stop
    end if
        
    flashData%model = model
    flashData%NC = NC
    flashData%calc_type = calc_type
    flashData%K = 1.0d0
    
    call initialization (model, NC, calc_type, T, P, feed, trace, gas, heavy, beta,&
                        flashData%K, flashData%w, flashData%Z, error_init)
    
    if (error_init) then
	
      write(ERROR_UNIT,*) "GCA Class: Initialization of flash data returned error"
	
    else
	
      flashData%T = T
      flashData%P = P
      flashData%feed = feed(1:NC)
      flashData%beta = beta
      flashData%trace = trace
      flashData%initialized = .TRUE. 
      
    end if
        
  end subroutine flashData_init
!---------------------------------------------------------------------------------
!
!
!
	subroutine flashData_calc (flashData, T,P, feed, trace, beta, var_spec, error_thermo)
    
	implicit none
    
	class(t_flashData), intent(inout)              :: flashData
	real(kind = 8), intent(inout)                  :: T, P
	real(kind = 8), dimension(:), intent(in)       :: feed
	logical, dimension(:), intent(in)              :: trace, var_spec
	real(kind = 8), intent(in)                     :: beta
	
	logical, intent(out)                           :: error_thermo
	
	real(kind = 8), dimension(NC,2)                :: w
	real(kind = 8), dimension(2)                   :: Z
	
	if (flashData%initialized) then
	  
	  if (flashData%calc_type == 'T') then
	    P = flashData%P
	  else
	    T = flashData%T
	  end if
	  flashData%var_spec = var_spec      
	  call gcaFlash (flashData%model, flashData%NC, flashData%calc_type, &
	                 T, P, feed, trace, gas, heavy, beta,                &
	                 error_thermo, flashData%K, w, flashData%Z,          &
	                 flashData%var_spec, flashData%dXi_dz, flashData%lnPhi)
	  
	  flashData%calcStatus = error_thermo
	  
	  if (.NOT. error_thermo) then
	    !update result data in flashData object
	    flashData%T = T
	    flashData%P = P
	    flashData%beta = beta
	    flashData%Spec(1:flashData%NC) = flashData%feed
	    flashData%xi(1:flashData%NC) = flashData%K
	    if (flashData%calc_type == 'T') then
		
	      flashData%Spec(flashData%NC+1) = T
	      flashData%xi(flashData%NC+1) = P
		
	    else
		
	      flashData%Spec(flashData%NC+1) = P
	      flashData%xi(flashData%NC+1) = T
		
	    end if
	    flashData%Spec(flashData%NC+2) = beta
	    flashData%w = w
	    flashData%feed = feed
	    flashData%trace = trace
	    
	  end if
	else
	  
	  write (ERROR_UNIT, *) "GCA Class: flashData_calc: Data is not initialized"
	  
	end if    
  
	end subroutine flashData_calc
!------------------------------------------------------------------------------------------------------------------- 
!
!
!
!
	subroutine flashData_extOpt (flashData, T, P, beta, feed)
  
	implicit none
	
	class(t_flashData)                      :: flashData  
	real(kind=8), intent(inout)             :: P, T
	real(kind=8), intent(in)                :: beta, feed(:)
	          
	integer                                 :: ii, kk, NC
	real(kind=8)                            :: totalFeed
	real(kind=8), dimension(flashData%NC+2) :: deltaS
	real(kind=8), dimension(flashData%NC+1) :: xiUpdate, kLoc
	
	if (flashData%calcStatus) then !if true => error_thermo or no initial full calculation
	  
	  write(ERROR_UNIT,*) "GCA class: flashData_extrapolate: Successful flash calculation must be performed befor extrapolation."
	  return
	  
	end if
	  
	NC = flashData%NC
	xiUpdate = 0.0d0	
	do kk = 1, NC ! composition
	  
	  if (flashData%var_spec(kk)) then
	    deltaS(kk) = feed(kk)-flashData%Spec(kk)
	  else
	    deltaS(kk) = 0D0
	  endif
	  
	end do	
	if (flashData%var_spec(NC+1)) then
	  if (flashData%calc_type == 'T') then
	    
	    deltaS(NC+1) = dlog(T/flashData%Spec(NC+1))
	    
	  elseif (flashData%calc_type == 'P') then
	    
	    deltaS(NC+1) = dlog(P/flashData%Spec(NC+1))
	    
	  end if
	else

	  deltaS(NC+1) = 0

	endif
	if (flashData%var_spec(NC+2)) then
	  
	  deltaS(NC+2) = beta - flashData%Spec(NC+2)       
	  
	else
	  
	  deltaS(NC+2) = 0
	  
	endif
      
	!intel MKL: matrix-vector multiplication
! 	call dgemv('N', n+1, n+1, 1.0d0, flashData%dXi_dS(1:n+1,1:n+1), n+1, deltaX, &
! 		    1, 0.0d0, xiUpdate, 1)
! 	intrinsic Fortran multiplication:
	xiUpdate = matmul(flashData%dxi_dz(:NC+1,:NC+2), deltaS(:NC+2)) 
	
	!intel MKL: exponential-function subroutine for vectors
! 	call vdexp(n+1,xiUpdate,kLoc)
	
	!intel MKL: element by element multiplication of two vectors
! 	call vdmul(n+1,flashData%K0(1:n+1),kLoc,kLoc)

! 	intrinsic Fortran exponentiation and multiplication
	flashData%K(:NC) = dexp(xiUpdate(:NC))*flashData%K(:NC)
  
  
!	flashData%k(1:n) = kLoc(1:n)  
	if (flashData%calc_type == 'T') then
	  
!	  P = kLoc(n+1)
	  P = dexp(xiUpdate(NC+1))*P
!	  flashData%p = kLoc(n+1)
	  
	elseif (flashData%calc_type == 'P') then
	  
!	  t = kLoc(n+1)
!	  flashData%t = kLoc(n+1)
	  T = dexp(xiUpdate(NC+1))*T
	  
	end if
	flashData%P = P
	flashData%T = T
	totalFeed = sum(feed(1:NC))
	do ii = 1, NC
	  
! 	  flashData%w(ii,L) = feed(ii)/(totalFeed*(1.0d0 + (kLoc(ii) - 1.0d0)*flashData%beta))
! 	  flashData%w(ii,V) = kLoc(ii)*flashData%w(ii,L)
	  flashData%w(ii,L) = feed(ii)/(totalFeed*(1.0d0 + (flashData%K(ii) - 1.0d0)*flashData%beta))
	  flashData%w(ii,V) = flashData%K(ii)*flashData%w(ii,L)
	  
	end do  
	
	end subroutine
!--------------------------------------------------------------------------------------------------------------
!
!
!
!
!
	subroutine flashData_print (flashData, outputFile)
  
	implicit none
  
	class(t_flashData), intent(in) :: flashData
	integer, intent(in)           :: outputFile
  
	call gcaWriteResuls(flashData%model, flashData%NC, flashData%calc_type,        &
	                    outputFile, flashData%calcStatus, flashData%T, flashData%P,&
	                    flashData%feed, flashData%trace, gas, heavy,               &
	                    flashData%beta, flashData%K, flashData%w, flashData%Z)

	end subroutine flashData_print
!--------------------------------------------------------------------------------------------------------------
!
!
!
!
!	
	subroutine trace_compounds (flashData, NC, tolerance)

	implicit none
	
	class(t_flashData), intent(inout) :: flashData
	real(8), intent(in), optional     :: tolerance
	integer, intent(in)               :: NC
	real(8)                           :: loc_tolerance
	integer                           :: i
	
	if (present(tolerance)) then
	  
	  loc_tolerance = tolerance
	  
	else
	  
	  loc_tolerance = epsilon(1.D0)
	  
	endif
	do i = 1, NC

	  !deactivation of trace compounds:
	  if (flashData%feed(i) <= loc_tolerance) then

	    flashData%feed(i) = loc_tolerance
	    flashData%trace(i) = .TRUE.	    

	  else
	    
	    flashData%trace(i) = .FALSE.
	    
	  endif

	enddo	
	
	return
	endsubroutine trace_compounds
!--------------------------------------------------------------------------------------------------------------
!
!
!
!
!
!   subroutine flashData_extExpLog(flashData,t,p,feed)
!   
!   implicit none
!   
!   class(t_flashData) :: flashData  
!   real(kind = 8), intent(inout) :: t
!   real(kind = 8), intent(inout) :: p
!   real(kind = 8), intent(in) :: feed(:)
!   
!   integer :: ii, kk, n
!   real(kind = 8) :: fSum
!   real(kind = 8), dimension(flashData%NC+1) ::  deltaS, xiUpdate, kLoc 
!   
!   if (flashData%calcStatus) then !if true => error_thermo or no initial full calculation
!     write(ERROR_UNIT,*) "GCA class: flashData_extrapolate: Successful flash calculation must be performed befor extrapolation."
!     return
!   end if
!     
!   n = flashData%NC
!   xiUpdate = 0.0d0
!   
!   do kk = 1,n ! composition
!     if (flashData%Spec(kk).lt.zeroTol) then
!       deltaS(kk) = feed(kk)
!     elseif (feed(kk).lt.zeroTol) then      
!       deltaS(kk) = -flashData%Spec(kk)
!     else
!       deltaS(kk) = dlog(feed(kk)/flashData%Spec(kk))*flashData%Spec(kk)
!     end if
!   end do
!   if (flashData%calc_type == 'T') then
!     deltaS(n+1) = dlog(t/flashData%Spec(n+1))
!   elseif (flashData%calc_type == 'P') then
!     deltaS(n+1) = dlog(p/flashData%Spec(n+1))
!   end if
!   
!   !intel MKL: matrix-vector multiplication
!   call dgemv('N',n+1,n+1,1.0d0,flashData%dXi_dz(1:n+1,1:n+1),n+1,deltaS,1,0.0d0,xiUpdate,1)
!   
!   !intel MKL: exponential-function subroutine for vectors
!   call vdexp(n+1,xiUpdate,kLoc)
!   
!   !intel MKL: element by element multiplication of two vectors
!   call vdmul(n+1,flashData%xi(1:n+1),kLoc,kLoc)
!   
!   flashData%k(1:n) = kLoc(1:n)  
!   if (flashData%calc_type == 'T') then
!     p = kLoc(n+1)
!     flashData%p = kLoc(n+1)
!     flashData%t = t
!   elseif (flashData%calc_type == 'P') then
!     t = kLoc(n+1)
!     flashData%t = kLoc(n+1)
!     flashData%p = p
!   end if
!   
!   fSum = sum(feed(1:n))
!   do ii = 1,n
!     flashData%w(ii,L) = feed(ii)/(fSum*(1.0d0+(kLoc(ii)-1.0d0)*flashData%beta))
!     flashData%w(ii,V) = kLoc(ii)*flashData%w(ii,L)
!   end do  
!   
!   end subroutine
!   
!   subroutine flashData_extExpLin(flashData,t,p,feed)
!   
!   implicit none
!   
!   class(t_flashData) :: flashData  
!   real(kind = 8), intent(inout) :: t
!   real(kind = 8), intent(inout) :: p
!   real(kind = 8), intent(in) :: feed(:)
!   
!   integer :: ii, kk, n
!   real(kind = 8) :: fSum
!   real(kind = 8), dimension(flashData%NC+1) :: deltaX, deltaS, xiUpdate, kLoc
!   
!   if (flashData%calcStatus) then !if true => error_thermo or no initial full calculation
!     write(ERROR_UNIT,*) "GCA class: flashData_extrapolate: Successful flash calculation must be performed befor extrapolation."
!     return
!   end if
!     
!   n = flashData%NC
!   xiUpdate = 0.0d0
!   
!   
!   do kk = 1,n ! composition
!       deltaX(kk) = feed(kk)-flashData%Spec(kk)
!   end do
!   
!   !use also linearized version for temperature and pressure (by dividing by t0 or p0 respectively)
!   if (flashData%calc_type == 'T') then
!     deltaX(n+1) = t/flashData%Spec(n+1)-1.0d0
!   elseif (flashData%calc_type == 'P') then
!     deltaX(n+1) = p/flashData%Spec(n+1)-1.0d0
!   end if
!       
!   !intel MKL: matrix-vector multiplication
!   call dgemv('N',n+1,n+1,1.0d0,flashData%dXi_dz(1:n+1,1:n+1),n+1,deltaX,1,0.0d0,xiUpdate,1)
!   
!   !intel MKL: exponential-function subroutine for vectors
!   call vdexp(n+1,xiUpdate,kLoc)
!   
!   !intel MKL: element by element multiplication of two vectors
!   call vdmul(n+1,flashData%xi(1:n+1),kLoc,kLoc)
!   
!   flashData%k(1:n) = kLoc(1:n)  
!   if (flashData%calc_type == 'T') then
!     p = kLoc(n+1)
!     flashData%p = kLoc(n+1)
!     flashData%t = t
!   elseif (flashData%calc_type == 'P') then
!     t = kLoc(n+1)
!     flashData%t = kLoc(n+1)
!     flashData%p = p
!   end if
!   
!   fSum = sum(feed(1:n))
!   do ii = 1,n
!     flashData%w(ii,L) = feed(ii)/(fSum*(1.0d0+(kLoc(ii)-1.0d0)*flashData%beta))
!     flashData%w(ii,V) = kLoc(ii)*flashData%w(ii,L)
!   end do  
!   
!   end subroutine
!   
!   subroutine flashData_extLinLin(flashData,t,p,feed)
!   
!   implicit none
!   
!   class(t_flashData) :: flashData  
!   real(kind = 8), intent(inout) :: t
!   real(kind = 8), intent(inout) :: p
!   real(kind = 8), intent(in) :: feed(:)
!   
!   integer :: ii, kk, n
!   real(kind = 8) :: fSum
!   real(kind = 8), dimension(flashData%NC+1) :: deltaX, xiUpdate, kLoc
!   
!   if (flashData%calcStatus) then !if true => error_thermo or no initial full calculation
!     write(ERROR_UNIT,*) "GCA class: flashData_extrapolate: Successful flash calculation must be performed befor extrapolation."
!     return
!   end if
!     
!   n = flashData%NC
!   xiUpdate = 0.0d0
!     
!   do kk = 1,n ! composition
!       deltaX(kk) = feed(kk)-flashData%Spec(kk)
!   end do
!   
!   if (flashData%calc_type == 'T') then
!     deltaX(n+1) = t/flashData%Spec(n+1)-1.0d0
!   elseif (flashData%calc_type == 'P') then
!     deltaX(n+1) = p/flashData%Spec(n+1)-1.0d0
!   end if
!       
!   !intel MKL: matrix-vector multiplication
!   call dgemv('N',n+1,n+1,1.0d0,flashData%dXi_dz(1:n+1,1:n+1),n+1,deltaX,1,0.0d0,xiUpdate,1)
!   
!   xiUpdate = xiUpdate+1.0d0
!     
!   !intel MKL: element by element multiplication of two vectors
!   call vdmul(n+1,flashData%xi(1:n+1),xiUpdate,kLoc)
!   
!   flashData%k(1:n) = kLoc(1:n)  
!   if (flashData%calc_type == 'T') then
!     p = kLoc(n+1)
!     flashData%p = kLoc(n+1)
!     flashData%t = t
!   elseif (flashData%calc_type == 'P') then
!     t = kLoc(n+1)
!     flashData%t = kLoc(n+1)
!     flashData%p = p
!   end if
!   
!   fSum = sum(feed(1:n))
!   do ii = 1,n
!     flashData%w(ii,L) = feed(ii)/(fSum*(1.0d0+(kLoc(ii)-1.0d0)*flashData%beta))
!     flashData%w(ii,V) = kLoc(ii)*flashData%w(ii,L)
!   end do  
!   
!   end subroutine
 
 
end module tpfClass
  
 
 
 
!program fractionFlash
!
!  use tpfClass
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
!  real(8), dimension(:,:), allocatable :: dXi_dz
!
!  character                            :: calc_type
!  character(len = strLen), dimension(:), allocatable :: flashInput
!  
!  character(10), dimension(NCM)        :: CompoundName
!  common /NAME/                           CompoundName
!
!  logical                              :: error_init, error_thermo = .false.
!  logical, dimension(:), allocatable   :: trace  
!
!  type(t_flashData)                    :: flashData
!  
!  !initialize GCA
!  call gcaInit('input.txt',flashInput)
!  
!  allocate (trace(NC), dXi_dz(NC+1,NC+1))
!  flashData = t_flashData(NC)
!
!
!! Succesive flashes will be performed, at user specified conditions until an error occure:
!! * Error while reading (iostat <> 0)
!! * T,P = 0
!! * Sum (n(i)) = 0)
! 
!!DIR$ if (verbose> = 2)
!  open (unit = logFile, file = 'log.txt')
!!DIR$ ENDIF
!  open (unit = outputFile, file = outputFilename)
!  write (outputFile, '(3X, "T(K)    P(atm)      V/F", 3X, 3(" -"), 6X, <NC>(A10, 5X), 5X, "Z")') CompoundName(:NC)
!
!  j = 1
!  do
!    read (flashInput(j), *, iostat = input_status) (n(i), i = 1, NC), T, P, beta, calc_type
!    Ntot = sum(n(:NC))
!    if (input_status /= 0 .OR. P <= 0 .OR. T <= 0 .OR. Ntot <= 0) then
!
!      exit !exit main loop and then stop.
!
!    endif
!    j = j+1
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
!!DIR$ if (verbose> = 2)
!    write (logFile, *)
!!DIR$ ENDIF
!    
!    call flashData%init(model, NC, calc_type, T, P, feed, trace, beta, error_init)
!    if (error_init) then
!      write(error_unit,'(A,I0,A)') "Error while initializing flash. Skipping input ",j,"."
!      CYCLE
!    end if
!      
!    call flashData%calc(T,P,feed,trace,beta,error_thermo)
!    
!    call flashData%print(outputFile)
!    
!    do i = 1, NC+1
!      write (outputFile, '(X, <NC+1>(G11.4, X))') flashData%dxi_dz(i,:NC+1)
!    enddo
!    write (outputFile, *)
!    do i = 1, 2
!      write (outputFile, '(" - - - ln(Ø) = [", 7X, <NC+1>(G11.4, X))') flashData%lnPhi(:NC, i)
!    enddo
!    write (outputFile, *)
!  enddo
!  
!!DIR$ if (verbose> = 2)
!  close (unit = logFile)
!!DIR$ ENDIF
!  close (unit = outputFile)
!  stop
!
!end program

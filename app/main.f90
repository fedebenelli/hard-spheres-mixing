program main
   use iso_fortran_env, only: pr => real64
   use hard_spheres_mix, only: fugacity, alloc
   implicit none

   integer :: nc
   integer :: i, n_deriv, t_deriv, p_deriv, phase_type

   integer :: funit_input

   call read_file("test.txt")
   call run
contains
   subroutine read_file(filename)
      use hard_spheres_mix, only: read_eos_parameters, alloc
      integer :: funit_output
      character(len=*), intent(in) :: filename

      open (newunit=funit_input, file=trim(filename))
      read (funit_input, *) nc
      call alloc(nc)
      call read_eos_parameters(funit_input, funit_output, nc)
   end subroutine

   subroutine run
      ! Input variables
      real(pr) :: p, t
      real(pr) :: n(nc)

      ! Output variables
      real(pr) :: v
      real(pr) :: lnphi(nc), dlnphi_dn(nc, nc), dlnphi_dt(nc), dlnphi_dp(nc)

      integer :: i, point

      read (funit_input, *) (n(i), i=1, nc), T, P
      phase_type = 1
      n_deriv = 2
      t_deriv = 1
      p_deriv = 1

      do i=1,100000


      call fugacity(nc, phase_type, &
                    n_deriv, t_deriv, p_deriv, &
                    T, P, n, &
                    v, lnPhi, dlnPhi_dn, dlnPhi_dT, dlnPhi_dP &
            )
      end do

      print "(G,2x,3(E14.5))", "lnPHI:     ", lnphi
      print "(G,2x,3(E14.5))", "dlnPHI_dT: ", dlnphi_dt
      print "(G,2x,3(E14.5))", "dlnPHI_dP: ", dlnphi_dp

      do i=1,nc
         print "(G,2x,I1,3(E14.5))", "dlnphi_dn:     ", i, dlnphi_dn(i, :)
      end do
   end subroutine
end program

! program example
!   use iso_fortran_env, only: pr => real64
!
!   implicit none
!   integer, parameter                   :: inputFile = 1
!
!   integer                              :: i, n_deriv, NC, p_deriv, phase_type, t_deriv
!   real(pr)                              :: P, T, V
!   real(pr), dimension(:), allocatable   :: ac, b, dlnPhi_dP, dlnPhi_dT, m, n, lnPhi, Tc
!   real(pr), dimension(:,:), allocatable :: dlnPhi_dn, k
!
!   open (unit = inputFile, file = "test.txt")
!   read (inputFile, *) NC
!   allocate (n(NC), lnPhi(NC), ac(NC), b(NC), m(NC), Tc(NC), k(NC,NC), dlnPhi_dP(NC), dlnPhi_dT(NC), dlnPhi_dn(NC,NC))
!   call read_eos_parameters (inputFile, NC)
!   read (inputFile, *) (n(i), i = 1, NC), T, P
! !   call evaluate_mixture_parameters (NC, T, n, Tc, ac, b, m, k)
!   phase_type = 1
!   n_deriv = 2
!   t_deriv = 1
!   p_deriv = 1
!   do i = 1, NC
!   if (n(i) <= 0._pr) n(i) = 1D-20
!   enddo
!   call fugacity (NC, phase_type, n_deriv, t_deriv, p_deriv, T, P, n, v, lnPhi, dlnPhi_dn, dlnPhi_dT, dlnPhi_dP)
!
!   stop
! endprogram example

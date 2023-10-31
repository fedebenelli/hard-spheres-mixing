module hard_spheres_mix
   use iso_fortran_env, only: pr => real64
   use hyperdual_mod
   implicit none

   real(pr), parameter :: R = 0.08314472

   type(hyperdual) :: amix, bmix, lambda(0:3)
   real(pr), allocatable :: ac(:), b(:), k(:), Tc(:), kij(:, :)

   type :: prop(n)
      integer, len :: n
      real(pr) :: val
      real(pr) :: dn(n)
      real(pr) :: dn2(n, n)

      real(pr) :: dt
      real(pr) :: dt2

      real(pr) :: dv
      real(pr) :: dv2

      real(pr) :: dtv

      real(pr) :: dnv(n)
      real(pr) :: dnt(n)
   end type

   contains
   subroutine alloc(n)
      integer, intent(in) :: n
      allocate (ac(n), b(n), k(n), Tc(n), kij(n, n))
   end subroutine

   type(hyperdual) function Helmholtz_TVn(T, V, n) result(Ar_RT)
      real(pr), parameter    :: xi = 4.0_pr
      type(hyperdual), intent(in) :: t, v, n(:)
      type(hyperdual)        :: Ar_RT_att, Ar_RT_rep
      type(hyperdual)        :: eta, l1l2_l3l0, l23_l0l32, logContribution, n_tot

      n_tot = sum(n(:))
      call evaluate_mixture_parameters(T, V, n)

      eta = bmix/v/xi

      l1l2_l3l0 = lambda(1)*lambda(2)/lambda(3)/lambda(0)
      l23_l0l32 = lambda(2)**3/lambda(0)/lambda(3)**2
      logContribution = log(1._pr - xi*eta)/xi

      Ar_RT_rep = -n_tot*((1._pr + 3._pr*l1l2_l3l0)*logContribution &
                          + 3._pr/xi*(l23_l0l32 - 1._pr/2._pr - l1l2_l3l0/2._pr) &
                          *(eta + logContribution))

      Ar_RT_att = -amix*log((v + bmix)/(v))/bmix/R/T
      Ar_RT = Ar_RT_rep + Ar_RT_att
   end function Helmholtz_TVn

   function ar_derivatives(t, v, n, comp_derivs) result(residual_helmholtz)
      real(pr), intent(in) :: t, v, n(:)
      type(prop(size(n))) :: residual_helmholtz

      type(hyperdual) :: t_hd, v_hd, n_hd(size(n)), ar_hd
      integer :: i, j, nc
      logical :: comp_derivs

      nc = size(n)

      ! First derivatives wrt T and V, crossed deriv and ar value
      call reset_vars
      t_hd%f1 = 1
      v_hd%f2 = 1

      ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)

      residual_helmholtz%val = ar_hd%f0
      residual_helmholtz%dt = ar_hd%f1
      residual_helmholtz%dv = ar_hd%f2
      residual_helmholtz%dtv = ar_hd%f12

      ! Second derivative wrt T
      call reset_vars
      t_hd%f1 = 1
      t_hd%f2 = 1
      ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)
      residual_helmholtz%dt2 = ar_hd%f12

      ! Second derivative wrt V
      call reset_vars
      v_hd%f1 = 1
      v_hd%f2 = 1
      ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)
      residual_helmholtz%dv2 = ar_hd%f12

      if (comp_derivs) then
         do i = 1, nc
            !dn and dn(i,i)2
            call reset_vars
            n_hd(i)%f1 = 1
            n_hd(i)%f2 = 1
            ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)
            residual_helmholtz%dn(i) = ar_hd%f1
            residual_helmholtz%dn2(i, i) = ar_hd%f12

            ! dnT
            call reset_vars
            n_hd(i)%f1 = 1
            t_hd%f2 = 1
            ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)
            residual_helmholtz%dnt(i) = ar_hd%f12

            ! dnV
            call reset_vars
            n_hd(i)%f1 = 1
            v_hd%f2 = 1
            ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)
            residual_helmholtz%dnv(i) = ar_hd%f12

            do j = i, nc
               call reset_vars
               n_hd(i)%f1 = 1
               n_hd(j)%f2 = 1
               ar_hd = Helmholtz_TVn(t_hd, v_hd, n_hd)
               residual_helmholtz%dn2(i, j) = ar_hd%f12
               residual_helmholtz%dn2(j, i) = ar_hd%f12
            end do
         end do
      end if
   contains
      subroutine reset_vars
         t_hd = t
         v_hd = v
         n_hd = n
      end subroutine
   end function

   function Pressure_TVn(NC, T, V, n, v_deriv) result(P)
      integer                        :: NC, v_deriv
      real(pr)                       :: delta, n_tot, T, V
      real(pr), dimension(2)         :: P
      real(pr), dimension(NC)        :: n

      type(prop(size(n))) :: residual_helmholtz
      n_tot = sum(n)

      residual_helmholtz = ar_derivatives(t, v, n, .false.)

      P(1) = R*T*(n_tot/V - residual_helmholtz%dv)
      if (v_deriv > 0) P(2) = -R*T*(n_tot/V/V + residual_helmholtz%dv2) !dP/dV
   end function Pressure_TVn

   subroutine evaluate_mixture_parameters(T, V, n)
      type(hyperdual), intent(in) :: t, v, n(:)

      type(hyperdual), dimension(size(n)) :: a, ai, z2
      type(hyperdual) :: nij

      integer :: i, j
      integer :: nc

      nc = size(n)
      amix = 0.0_pr
      bmix = 0.0_pr
      lambda = 0.0_pr

      a = ac*(1.0_pr + k*(1.0_pr - sqrt(t/tc)))**2
      ! ai = ac*a*a
      ! a = sqrt(ai)
      ! z2 = n*n

      do i = 1, nc
         do j = 1, nc
            nij = n(i)*n(j)
            amix = amix + nij * sqrt(a(i)*a(j)) * (1 - kij(i, j))
            ! bmix = bmix + nij*(b(i) + b(j)) *(1 - lij(i, j))
         end do
         lambda(1) = lambda(1) + n(i)*b(i)**(1.0_pr/3.0_pr)
         lambda(2) = lambda(2) + n(i)*b(i)**(2.0_pr/3.0_pr)
      end do
      ! amix = 2*amix + sum(z2*ai)
      bmix = sum(n*b)
      ! bmix = (bmix + sum(z2*b))/sum(n)

      lambda(0) = sum(n)
      lambda(3) = bmix

      ! lambda(1) = dot_product(n, b**(1.0_pr/3))
      ! lambda(2) = dot_product(n, b**(2.0_pr/3))
      ! lambda(3) = bMix
   end subroutine evaluate_mixture_parameters

   subroutine get_covolume(b)
      !--------------------------------------------------------------------------------------------
      !
      !  This subroutine gives back the minimum possible volume for the EoS. For most common cubic
      !  EoS this is the clasical covolume. For other EoS, it may be 4 times the van der Waals co-
      !  volume.
      !
      real(pr)                :: b
      b = bMix%f0/1000
      return
   end subroutine get_covolume

   subroutine read_eos_parameters(inputFile, outputFile, NC)
      implicit none
      integer, parameter          :: NCM = 15
      integer                     :: compound, inputFile, NC, outputFile
      real(pr)                    :: aux
      real(pr), dimension(NC)     :: omega, Pc
      real(8),  dimension(2)      :: Bgpec, TCgpec, PCgpec, DCgpec

      COMMON/CRIT/TCgpec, PCgpec, DCgpec
      COMMON/COVOL/Bgpec

      do compound = 1, NC
         read (inputFile, *) Tc(compound), Pc(compound), omega(compound)
         aux = (2**(1._pr/3) - 1._pr)/3
         ac(compound) = (R*Tc(compound))**2/Pc(compound)/aux/27
         b(compound) = aux*R*Tc(compound)/Pc(compound)
         k(compound) = 0.48_pr + omega(compound)*(1.574_pr - .176_pr*omega(compound))
         if (compound <= 2) then
            TCgpec(compound) = Tc(compound)
            PCgpec(compound) = Pc(compound)*1.01325D0
            DCgpec(compound) = PCgpec(compound)*3/0.08314D0/TCgpec(compound)
            Bgpec(compound) = b(compound)*1.01325D-3
         end if
      end do

      kij = 0
      do compound = 1, NC - 1
         read (inputFile, *) kij(compound + 1:NC, compound)
         kij(compound, compound + 1:NC) = kij(compound + 1:NC, compound)
      end do
   end subroutine

   subroutine fugacity(NC, phase_type, n_deriv, t_deriv, p_deriv, T, P, n, v, lnPhi, dlnPhi_dn, dlnPhi_dT, dlnPhi_dP)
      implicit none
      integer :: i, j, n_deriv, NC, phase_type, p_deriv, t_deriv
      real(pr), intent(in) :: t, p, n(:)

      real(pr)                      :: dP_dT, dP_dV, n_tot, RT, V, Z
      real(pr), dimension(size(n))       :: dlnPhi_dP, dlnPhi_dT, dP_dn, lnPhi
      real(pr), dimension(size(n), size(n))   :: dlnPhi_dn

      type(prop(nc)) :: residual_helmholtz

      ! Function from which all thermodynamic properties are calculated.
      call volume_calculation(V)
      n_tot = sum(n(:NC))
      Z = P*V/n_tot/R/T

      residual_helmholtz = ar_derivatives(t, v, n, .true.)

      ! fugacity
      lnPhi = residual_helmholtz%dn - log(Z)

      ! Fugacity derivatives
      RT = R*T
      dP_dn =  RT*(1._pr/V - residual_helmholtz%dnv)
      dP_dV = -RT*(n_tot/V/V + residual_helmholtz%dv2)
      dP_dT = P/T - RT*residual_helmholtz%dtv
      
      dlnPhi_dT = residual_helmholtz%dnt + 1._pr/T + dP_dT*dP_dn/dP_dV/RT
      dlnPhi_dP = -dP_dn/RT/dP_dV - 1._pr/P
      do i = 1, NC
         do j = i, NC
            dlnPhi_dn(j, i) = residual_helmholtz%dn2(j, i) &
                              + 1._pr/n_tot + dP_dn(j)*dP_dn(i)/RT/dP_dV
            dlnPhi_dn(i, j) = dlnPhi_dn(j, i)
         end do
      end do
   end subroutine fugacity

   subroutine volume_calculation(v)
      implicit none
      real(pr) :: V
      ! v = 0.0097248811017998961_pr*82.05_pr*334.57975847415042_pr
      V=0.11655858674454586088325610656019
   end subroutine volume_calculation

   ! subroutine Helmholtz_CPSAFT(n_deriv, t_deriv, n, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)
   !    implicit none
   !    integer, parameter             :: NC = 2
   !    real(8), parameter             :: R = 0.08314D0
   !    real(pr), parameter            :: quad_delta = 1.0e-8_pr

   !    integer                        :: i, n_deriv, nVar, t_deriv

   !    real(8)                        :: Ar, ArT, ArTV, ArV, ArV2, RT, T, V
   !    real(8), dimension(NC)         :: Arn, ArTn, ArVn, n
   !    real(8), dimension(NC, NC)      :: Arn2

   !    real(pr)                       :: Helmholtz_TVn
   !    real(pr), dimension(NC + 2)      :: quad_df_dx, quad_x
   !    real(pr), dimension(NC + 2, NC + 2) :: quad_d2f_dx2

   !    ! external                       :: Helmholtz_TVn

   !    RT = R*T

   !    nVar = NC + 2
   !    quad_df_dx(2:NC + 2) = 1     !d(Ar/RT)/dV + d(Ar/RT)/dni
   !    quad_d2f_dx2(2, 2:NC + 2) = 1 !d2(Ar/RT)/dV2 + d2(Ar/RT)/dVdni
   !    quad_d2f_dx2(2:NC + 2, 2) = 1 !d2(Ar/RT)/dV2 + d2(Ar/RT)/dVdni

   !    if (n_deriv > 1) quad_d2f_dx2(3:NC + 2, 3:NC + 2) = 1
   !    if (t_deriv > 0) then
   !       quad_df_dx(1) = 1
   !       quad_d2f_dx2(1, :NC + 2) = 1
   !       quad_d2f_dx2(:NC + 2, 1) = 1
   !    end if
   !    quad_x(1) = (T)
   !    quad_x(2) = (V)*1q3
   !    do i = 1, NC

   !       quad_x(i + 2) = max(n(i), 1Q - 10)

   !    end do
   !    call central_finite_difference(quad_delta, nVar, Helmholtz_TVn, quad_x, quad_df_dx, quad_d2f_dx2)

   !    Ar = RT*dble(Helmholtz_TVn(nVar, quad_x(1), quad_x(2), quad_x(3:NC)))
   !    ArV = RT*dble(quad_df_dx(2))*1D3
   !    ArV2 = RT*dble(quad_d2f_dx2(2, 2))*1D6
   !    Arn(:NC) = RT*dble(quad_df_dx(3:NC + 2))
   !    ArVn(:NC) = RT*dble(quad_d2f_dx2(2, 3:NC + 2))*1D3
   !    if (n_deriv > 1) Arn2(:NC, :NC) = RT*quad_d2f_dx2(3:NC + 2, 3:NC + 2)
   !    if (t_deriv > 0) then

   !       ArT = Ar/T + RT*dble(quad_df_dx(1))
   !       ArTV = ArV/T + RT*dble(quad_d2F_dx2(1, 2))*1D3
   !       ArTn(:NC) = Arn(:NC)/T + RT*dble(quad_d2F_dx2(1, 3:NC + 2))

   !    end if
   ! end subroutine Helmholtz_CPSAFT
end module

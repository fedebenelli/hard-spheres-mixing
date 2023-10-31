module teststuff
   ! real(16), parameter :: R = 82.05Q0
   real(16), parameter :: R = 0.08314472q0

contains
   real(16) function Helmholtz_TVn(nVar, T, V, n) result(Ar_RT)
      implicit none

      real(16), parameter         :: xi = 4

      integer                :: NC, nVar
      real(16)                :: Ar_RT_att, Ar_RT_rep, eta, l1l2_l3l0, l23_l0l32, logContribution, n_tot, T, V
      real(16), dimension(nVar - 2) :: n

      real(16)                :: a, b, lambda(0:3)
      common/EoS_parameter/a, b, lambda

      NC = nVar - 2
      n_tot = sum(n(:NC))
      lambda(0) = n_tot
      call evaluate_mixture_parameters(NC, T, V, n)
      eta = b/v/xi

      l1l2_l3l0 = lambda(1)*lambda(2)/lambda(3)/lambda(0)
      l23_l0l32 = lambda(2)**3/lambda(0)/lambda(3)**2
      logContribution = qlog(1.0Q0 - xi*eta)/xi
      Ar_RT_rep = -n_tot*((1q0 + 3*l1l2_l3l0)*logContribution + 3q0/xi*(l23_l0l32 - 1q0/2 - l1l2_l3l0/2)*(eta + logContribution))

      Ar_RT_att = -a*qlog((v + b)/(v))/b/R/T

      Ar_RT = Ar_RT_rep + Ar_RT_att

      return
   end function Helmholtz_TVn

   function Pressure_TVn(NC, T, V, n, v_deriv) result(P)
      implicit none

      integer                       :: NC, v_deriv

      real(16)                       :: delta, n_tot, T, V
      real(16), dimension(2)         :: P
      real(16), dimension(NC)        :: n
      real(16), dimension(NC + 2)      :: df_dx, x
      real(16), dimension(NC + 2, NC + 2) :: d2F_dx2

      n_tot = sum(n)
      x(1) = T
      x(2) = V
      x(3:NC + 2) = n(:NC)
      df_dx = 0
      df_dx(2) = 1
      d2F_dx2 = 0
      if (v_deriv > 0) d2f_dx2(2, 2) = 1
      delta = 1D-5

      call central_finite_difference(delta, NC, Helmholtz_TVn, x, df_dx, d2f_dx2)

      P(1) = R*T*(n_tot/V - df_dx(2))
      if (v_deriv > 0) P(2) = -R*T*(n_tot/V/V + d2F_dx2(2, 2)) !dP/dV
      return

   end function Pressure_TVn

   subroutine evaluate_mixture_parameters(NC, T, V, n)

      implicit none
      integer, parameter           :: NCM = 15
      integer                      :: i, j, NC

      real(16)                     :: V, T, Told
      real(16), dimension(NC)      :: alpha, n
      real(16), dimension(NCM, NCM) :: a

      DATA Told/0.Q0/
      save                         :: a

      real(16)                     :: aMix, bMix, lambda(0:3)
      real(16), dimension(NCM)     :: ac, b, m, Tc
      real(16), dimension(NCM, NCM) :: k
      common/EoS_parameter_pure/ac, b, m, k, Tc
      common/EoS_parameter/aMix, bMix, lambda

      !   pure
      IF (T /= Told) then

         Told = T
         do i = 1, NC

            alpha(i) = (1.Q0 + m(i)*(1.Q0 - qsqrt(T/Tc(i))))**2
            a(i, i) = ac(i)*alpha(i)

         end do

      end if

      !   mixture
      do i = 1, NC

         do j = i + 1, NC

            a(i, j) = (1Q0 - k(i, j))*qsqrt(a(i, i)*a(j, j))
            a(j, i) = a(i, j)

         end do

      end do
      aMix = dot_product(n(:NC), matmul(a(:NC, :NC), n(:NC)))
      bMix = dot_product(n(:NC), b(:NC))
      lambda(1) = dot_product(n(:NC), b(:NC)**(1q0/3))
      lambda(2) = dot_product(n(:NC), b(:NC)**(2q0/3))
      lambda(3) = bMix

      return

   end subroutine evaluate_mixture_parameters
   !--------------------------------------------------------------------------------------------
   !
   !  This subroutine gives back the minimum possible volume for the EoS. For most common cubic
   !  EoS this is the clasical covolume. For other EoS, it may be 4 times the van der Waals co-
   !  volume.
   !
   subroutine get_covolume(b)

      implicit none

      real(16)                :: b
      !Variables in COMMON verctor:
      real(16)                :: aMix, bMix
      common/EoS_parameter/aMix, bMix

      b = bMix/1000

      return

   end subroutine get_covolume
   ! subroutine mixing_rule ()
   !
   !   implicit none
   !
   !   return
   !
   ! endsubroutine mixing_rule

   subroutine read_eos_parameters(inputFile, outputFile, NC)

      implicit none
      integer, parameter          :: NCM = 15
      integer                     :: compound, inputFile, NC, outputFile
      real(16)                    :: aux
      real(16), dimension(NC)     :: omega, Pc

      real(16), dimension(NCM)    :: ac, b, m, Tc
      real(16), dimension(NCM, NCM):: k
      common/EoS_parameter_pure/ac, b, m, k, Tc

      ! real(8), dimension(2)       :: Bgpec, TCgpec, PCgpec, DCgpec
      ! COMMON/CRIT/TCgpec, PCgpec, DCgpec
      ! COMMON/COVOL/Bgpec

      do compound = 1, NC

         read (inputFile, *) Tc(compound), Pc(compound), omega(compound)
         aux = (2**(1Q0/3) - 1Q0)/3
         ac(compound) = (R*Tc(compound))**2/Pc(compound)/aux/27
         b(compound) = aux*R*Tc(compound)/Pc(compound)
         m(compound) = 0.48Q0 + omega(compound)*(1.574Q0 - .176Q0*omega(compound))
         ! TCgpec(compound) = Tc(compound)
         ! PCgpec(compound) = Pc(compound)*1.01325D0
         ! DCgpec(compound) = PCgpec(compound)*3/0.08314D0/TCgpec(compound)
         ! Bgpec(compound) = b(compound)*1.01325D-3

      end do

      k = 0
      do compound = 1, NC - 1
         read (inputFile, *) k(compound + 1:NC, compound)
         k(compound, compound + 1:NC) = k(compound + 1:NC, compound)
      end do

      return

   end subroutine
   !--------------------------------------------------------------------------------------------
   !
   !
   !
   subroutine central_finite_difference(delta, nVar, func, x, df_dx, d2f_dx2)

      implicit none

      integer                       :: i, j, nVar
      real(16)                       :: delta, f, func, f_back_backward, f_back_forward, f_backward, f_for_backward, &
                                        f_for_forward, f_forward
      real(16), dimension(nVar)      :: df_dx, x, x_pert
      real(16), dimension(nVar, nVar) :: d2f_dx2
      interface
         real(16) function fun(nVar, T, V, n) result(Ar_RT)
            integer :: nvar
            real(16) :: t, v, n(nvar-2)
         end function
      end interface

      do i = 1, nVar

         if (df_dx(i) /= 0 .OR. d2F_dx2(i, i) /= 0) then

            x_pert = x
            x_pert(i) = (1Q0 + delta)*x(i)
            f_forward = func(nVar, x_pert(1), x_pert(2), x_pert(3:nVar))

            x_pert = x
            x_pert(i) = (1Q0 - delta)*x(i)
            f_backward = func(nVar, x_pert(1), x_pert(2), x_pert(3:nVar))

            if (df_dx(i) /= 0) df_dx(i) = (f_forward - f_backward)/2/delta/x(i)
            if (d2f_dx2(i, i) /= 0) then

               f = func(nVar, x(1), x(2), x(3:nVar))
               d2F_dx2(i, i) = (f_forward - 2*f + f_backward)/delta/delta/x(i)/x(i)

            end if

         end if
         do j = i + 1, nVar

            if (d2f_dx2(i, j) /= 0) then

               x_pert = x
               x_pert(i) = (1Q0 - delta)*x(i)
               x_pert(j) = (1Q0 + delta)*x(j)
               f_back_forward = func(nVar, x_pert(1), x_pert(2), x_pert(3:nVar))

               x_pert(j) = (1Q0 - delta)*x(j)
               f_back_backward = func(nVar, x_pert(1), x_pert(2), x_pert(3:nVar))

               x_pert(i) = (1Q0 + delta)*x(i)
               f_for_backward = func(nVar, x_pert(1), x_pert(2), x_pert(3:nVar))

               x_pert(j) = (1Q0 + delta)*x(j)
               f_for_forward = func(nVar, x_pert(1), x_pert(2), x_pert(3:nVar))

               d2f_dx2(i, j) = (f_for_forward - f_for_backward - f_back_forward + f_back_backward)/4/delta/delta/x(i)/x(j)
               d2f_dx2(j, i) = d2f_dx2(i, j)

            end if

         end do

      end do

      return

   end subroutine central_finite_difference

   subroutine fugacity(NC, phase_type, n_deriv, t_deriv, p_deriv, T, P, n, v, lnPhi, dlnPhi_dn, dlnPhi_dT, dlnPhi_dP)

      implicit none

      real(16), parameter            :: delta = 1.Q-8

      integer                       :: i, j, n_deriv, NC, phase_type, p_deriv, t_deriv

      real(16)                       :: dP_dT, dP_dV, n_tot, P, RT, T, V, Z
      real(16), dimension(NC)        :: dlnPhi_dP, dlnPhi_dT, dP_dn, lnPhi, n
      real(16), dimension(NC + 2)      :: dF_dx, x
      real(16), dimension(NC, NC)     :: dlnPhi_dn
      real(16), dimension(NC + 2, NC + 2) :: d2F_dx2
      real(16) :: amix, bmix, lambda(0:3)
      common/EoS_parameter/aMix, bMix, lambda

      call volume_calculation(V)
      n_tot = sum(n(:NC))
      Z = P*V/n_tot/R/T

      dF_dx(1:2) = 0
      dF_dx(3:NC + 2) = 1 !calculation of number of moles derivatives)
      x(1) = T
      x(2) = V
      x(3:NC + 2) = n
      d2F_dx2 = 0
      if (t_deriv > 0 .OR. p_deriv > 0 .OR. n_deriv > 1) then

         d2F_dx2(2, 2) = 1
         d2f_dx2(2, 3:NC + 2) = 1; d2f_dx2(3:NC + 2, 2) = 1
         if (t_deriv > 0) d2F_dx2(1, 2) = 1; d2F_dx2(2, 1) = 1
         if (n_deriv > 1) d2F_dx2(2:NC + 2, 3:NC + 2) = 1

      end if
      call central_finite_difference(delta, NC + 2, Helmholtz_TVn, x, dF_dx, d2F_dx2)

      ! fugacity
      lnPhi(:NC) = dF_dx(3:NC + 2) - qlog(Z)

      ! Fugacity derivatives
      RT = R*T
      if (t_deriv > 0 .OR. p_deriv > 0 .OR. n_deriv > 1) then

         dP_dn(:NC) = RT*(1Q0/V - d2F_dx2(2, 3:NC + 2))
         dP_dV = -RT*(n_tot/V/V + d2F_dx2(2, 2))

         if (t_deriv > 0) then

            dP_dT = P/T - RT*d2F_dx2(1, 2)
            dlnPhi_dT(:NC) = d2f_dx2(1, 3:NC + 2) + 1Q0/T + dP_dT*dP_dn(:NC)/dP_dV/RT

         end if
         if (p_deriv > 0) dlnPhi_dP(:NC) = -dP_dn(:NC)/RT/dP_dV - 1Q0/P
         if (n_deriv > 1) then
            do i = 1, NC
               do j = i, NC

                  dlnPhi_dn(j, i) = d2F_dx2(2 + j, 2 + i) + 1Q0/n_tot + dP_dn(j)*dP_dn(i)/RT/dP_dV
                  dlnPhi_dn(i, j) = dlnPhi_dn(j, i)

               end do
            end do
         end if

      end if

      return

   end subroutine fugacity

   ! subroutine enthalpy ()
   !
   !   implicit none
   !
   !   external _Helmholtz
   !
   !   call central_finite_difference (delta, _Helmholtz,
   !
   ! endsubroutine enthalpy

   subroutine volume_calculation(v)

      implicit none
      real(16) :: V
      ! v = 0.0097248811017998961Q0*82.05Q0*334.57975847415042Q0
      v = 0.11655858674454586088325610656019q0

      return

   end subroutine volume_calculation

   subroutine Helmholtz_CPSAFT(n_deriv, t_deriv, n, V, T, Ar, ArV, ArTV, ArV2, Arn, ArVn, ArTn, Arn2)

      implicit none
      integer, parameter             :: NC = 2
      real(16), parameter            :: quad_delta = 1.Q-8

      integer                        :: i, n_deriv, nVar, t_deriv

      real(8)                        :: Ar, ArT, ArTV, ArV, ArV2, RT, T, V
      real(8), dimension(NC)         :: Arn, ArTn, ArVn, n
      real(8), dimension(NC, NC)      :: Arn2

      real(16), dimension(NC + 2)      :: quad_df_dx, quad_x
      real(16), dimension(NC + 2, NC + 2) :: quad_d2f_dx2

      RT = R*T

      nVar = NC + 2
      quad_df_dx(2:NC + 2) = 1     !d(Ar/RT)/dV + d(Ar/RT)/dni
      quad_d2f_dx2(2, 2:NC + 2) = 1 !d2(Ar/RT)/dV2 + d2(Ar/RT)/dVdni
      quad_d2f_dx2(2:NC + 2, 2) = 1 !d2(Ar/RT)/dV2 + d2(Ar/RT)/dVdni

      if (n_deriv > 1) quad_d2f_dx2(3:NC + 2, 3:NC + 2) = 1
      if (t_deriv > 0) then

         quad_df_dx(1) = 1
         quad_d2f_dx2(1, :NC + 2) = 1
         quad_d2f_dx2(:NC + 2, 1) = 1

      end if

      quad_x(1) = (T)
      quad_x(2) = (V)*1q3
      do i = 1, NC

         quad_x(i + 2) = max(n(i), 1.q-10)

      end do
      call central_finite_difference(quad_delta, nVar, Helmholtz_TVn, quad_x, quad_df_dx, quad_d2f_dx2)

      Ar = RT*dble(Helmholtz_TVn(nVar, quad_x(1), quad_x(2), quad_x(3:NC)))
      ArV = RT*dble(quad_df_dx(2))*1D3
      ArV2 = RT*dble(quad_d2f_dx2(2, 2))*1D6
      Arn(:NC) = RT*dble(quad_df_dx(3:NC + 2))
      ArVn(:NC) = RT*dble(quad_d2f_dx2(2, 3:NC + 2))*1D3
      if (n_deriv > 1) Arn2(:NC, :NC) = RT*quad_d2f_dx2(3:NC + 2, 3:NC + 2)
      if (t_deriv > 0) then

         ArT = Ar/T + RT*dble(quad_df_dx(1))
         ArTV = ArV/T + RT*dble(quad_d2F_dx2(1, 2))*1D3
         ArTn(:NC) = Arn(:NC)/T + RT*dble(quad_d2F_dx2(1, 3:NC + 2))

      end if

   end subroutine Helmholtz_CPSAFT
end module teststuff

program example
   use teststuff
   implicit none
   integer, parameter                   :: inputFile = 1

   integer                              :: i, n_deriv, NC, p_deriv, phase_type, t_deriv
   real(16)                              :: P, T, V
   real(16), dimension(:), allocatable   :: ac, b, dlnPhi_dP, dlnPhi_dT, m, n, lnPhi, Tc
   real(16), dimension(:, :), allocatable :: dlnPhi_dn, k

   open (unit=inputFile, file="test.txt")
   read (inputFile, *) NC
   allocate (n(NC), lnPhi(NC), ac(NC), b(NC), m(NC), Tc(NC), k(NC, NC), dlnPhi_dP(NC), dlnPhi_dT(NC), dlnPhi_dn(NC, NC))
   call read_eos_parameters(inputFile, 2, NC)
   read (inputFile, *) (n(i), i=1, NC), T, P
   !   call evaluate_mixture_parameters (NC, T, n, Tc, ac, b, m, k)
   phase_type = 1
   n_deriv = 2
   t_deriv = 1
   p_deriv = 1
   do i = 1, NC
      if (n(i) <= 0Q0) n(i) = 1D-20
   end do
   do i=1,100000
      call fugacity(NC, phase_type, n_deriv, t_deriv, p_deriv, T, P, n, v, lnPhi, dlnPhi_dn, dlnPhi_dT, dlnPhi_dP)
   end do
   print "(G,2x,3(E14.5))", "lnPHI:     ", lnphi
   print "(G,2x,3(E14.5))", "dlnPHI_dT: ", dlnphi_dt
   print "(G,2x,3(E14.5))", "dlnPHI_dP: ", dlnphi_dp

   do i=1,nc
      print "(G,2x,I1,3(E14.5))", "dlnphi_dn:     ", i, dlnphi_dn(i, :)
   end do

   stop

end program example

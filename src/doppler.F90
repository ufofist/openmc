module doppler

  use ace_header, only: Nuclide
  use constants, only: ZERO, ONE, TWO, FOUR, PI, K_BOLTZMANN, N_GAUSS_HERMITE
  use error, only: fatal_error
  use global, only: gh_abscissa, gh_weight

  implicit none

  real(8), parameter :: sqrt_pi_inv = ONE / sqrt(PI)

contains

!===============================================================================
! BROADEN determines whether to use the standard SIGMA1 method or a
! Gauss-Hermite quadrature for performing Doppler broadening and then calls the
! appropriate function
!===============================================================================

  function broaden(nuc, xs, i_energy, E, T) result(sigma)
    type(Nuclide), intent(in) :: nuc      ! nuclide
    real(8),       intent(in) :: xs(:)    ! unbroadened cross section
    integer,       intent(in) :: i_energy ! index on energy grid
    real(8),       intent(in) :: E        ! incoming neutron energy
    real(8),       intent(in) :: T        ! temperature (difference)
    real(8)                   :: sigma    ! broadened cross section

    real(8) :: gh_width  ! half-width of Gauss-Hermite quadrature
    real(8) :: xc        ! sqrt(A/kT)

    ! determine width of Gauss-Hermite quadrature and expand by 25%
    gh_width = 1.25_8*max(FOUR, gh_abscissa(N_GAUSS_HERMITE))

    ! Since nuc % sqrtE is just sqrt(E) for each energy on the grid, we have to
    ! multiply by sqrt(A/kT) to get to the form we want, sqrt(AE/kT)
    xc = sqrt(nuc % awr / K_BOLTZMANN * T)

    ! Check whether neutron energy is within 4 relative velocity units of 0 eV
    ! or if it is within 4 relative velocity units of 1 MeV (assumed to be the
    ! cutoff between the resolved and unresolved resonance ranges)
    if (nuc % sqrtE(i_energy - 1) * xc > gh_width .and. &
         nuc % sqrtE(i_energy + 1) * xc + gh_width < xc) then
      ! Calculate broadened cross section using a Gauss-Hermite quadrature
      sigma = gauss_hermite(nuc % sqrtE, xs, i_energy, E, xc)
    else
      ! Calculate broadened cross section using the SIGMA1 method
      sigma = sigma1(nuc % sqrtE, xs, i_energy, E, xc)
    end if
  end function broaden

!===============================================================================
! SIGMA1 takes a microscopic cross section at a temperature T_1 and Doppler
! broadens it to a higher temperature T_2 based on a method originally developed
! by Cullen and Weisbin (see "Exact Doppler Broadening of Tabulated Cross
! Sections," Nucl. Sci. Eng. 60, 199-229 (1976)). The only difference here is
! the F functions are evaluated based on complementary error functions rather
! than error functions as is done in the BROADR module of NJOY.
!===============================================================================

  function sigma1(sqrtE, xs, i_energy, E, xc) result(sigma)
    real(8), intent(in)  :: sqrtE(:)    ! sqrt(energy grid)
    real(8), intent(in)  :: xs(:)       ! unbroadened cross section
    integer, intent(in)  :: i_energy    ! index on energy grid
    real(8), intent(in)  :: E           ! incoming neutron energy
    real(8), intent(in)  :: xc
    real(8)              :: sigma       ! broadened cross section

    integer              :: k        ! loop indices
    integer              :: n        ! number of energy points
    real(8)              :: F_a(0:4) ! F(a) functions as per C&W
    real(8)              :: F_b(0:4) ! F(b) functions as per C&W
    real(8)              :: H(0:4)   ! H functions as per C&W
    real(8)              :: y        ! proportional to neutron velocity
    real(8)              :: y_sq     ! y**2
    real(8)              :: y_inv    ! 1/y
    real(8)              :: y_inv_sq ! 1/y**2
    real(8)              :: slope    ! slope of xs between adjacent points
    real(8)              :: Ak, Bk   ! coefficients at each point
    real(8)              :: a, b     ! values of x(k)-y and x(k+1)-y

    n = size(sqrtE)

    sigma    = ZERO
    y        = xc * sqrt(E)
    y_sq     = y*y
    y_inv    = ONE / y
    y_inv_sq = y_inv / y

    ! =======================================================================
    ! EVALUATE FIRST TERM FROM x(k) - y = 0 to -4

    k = i_energy
    a = ZERO
    call calculate_F(F_a, a)

    do while (a >= -FOUR .and. k > 1)
      ! Move to next point
      F_b = F_a
      k = k - 1
      a = xc*sqrtE(k) - y

      ! Calculate F and H functions
      call calculate_F(F_a, a)
      H = F_a - F_b

      ! Calculate A(k), B(k), and slope terms
      Ak = y_inv_sq*H(2) + TWO*y_inv*H(1) + H(0)
      Bk = y_inv_sq*H(4) + FOUR*y_inv*H(3) + 6.0_8*H(2) + FOUR*y*H(1) + y_sq*H(0)
      slope = (xs(k+1) - xs(k)) / (xc**2*(sqrtE(k+1)**2 - sqrtE(k)**2))

      ! Add contribution to broadened cross section
      sigma = sigma + Ak*(xs(k) - slope*xc**2*sqrtE(k)**2) + slope*Bk
    end do

    ! =======================================================================
    ! EXTEND CROSS SECTION TO 0 ASSUMING 1/V SHAPE

    if (k == 1 .and. a >= -FOUR) then
      ! Since x = 0, this implies that a = -y
      F_b = F_a
      a = -y

      ! Calculate F and H functions
      call calculate_F(F_a, a)
      H = F_a - F_b

      ! Add contribution to broadened cross section
      sigma = sigma + xs(k)*xc*sqrtE(k)*(y_inv_sq*H(1) + y_inv*H(0))
    end if

    ! =======================================================================
    ! EVALUATE FIRST TERM FROM x(k) - y = 0 to 4

    k = i_energy
    b = ZERO
    call calculate_F(F_b, b)

    do while (b <= FOUR .and. k < n)
      ! Move to next point
      F_a = F_b
      k = k + 1
      b = xc*sqrtE(k) - y

      ! Calculate F and H functions
      call calculate_F(F_b, b)
      H = F_a - F_b

      ! Calculate A(k), B(k), and slope terms
      Ak = y_inv_sq*H(2) + TWO*y_inv*H(1) + H(0)
      Bk = y_inv_sq*H(4) + FOUR*y_inv*H(3) + 6.0_8*H(2) + FOUR*y*H(1) + y_sq*H(0)
      slope = (xs(k) - xs(k-1)) / (xc**2*(sqrtE(k)**2 - sqrtE(k-1)**2))

      ! Add contribution to broadened cross section
      sigma = sigma + Ak*(xs(k) - slope*xc**2*sqrtE(k)**2) + slope*Bk
    end do

    ! =======================================================================
    ! EXTEND CROSS SECTION TO INFINITY ASSUMING CONSTANT SHAPE

    ! TODO: Need to change k == n
    if (k == n .and. b <= FOUR) then
      ! Calculate F function at last energy point
      a = xc*sqrtE(k) - y
      call calculate_F(F_a, a)

      ! Add contribution to broadened cross section
      sigma = sigma + xs(k) * (y_inv_sq*F_a(2) + TWO*y_inv*F_a(1) + F_a(0))
    end if

    ! =======================================================================
    ! EVALUATE SECOND TERM FROM x(k) + y = 0 to +4

    if (y <= FOUR) then
      ! Swap signs on y
      y = -y
      y_inv = -y_inv
      k = 1

      ! Calculate a and b based on 0 and x(1)
      a = -y
      b = xc*sqrtE(k) - y

      ! Calculate F and H functions
      call calculate_F(F_a, a)
      call calculate_F(F_b, b)
      H = F_a - F_b

      ! Add contribution to broadened cross section
      sigma = sigma - xs(k) * xc*sqrtE(k) * (y_inv_sq*H(1) + y_inv*H(0))

      ! Now progress forward doing the remainder of the second term
      do while (b <= FOUR)
        ! Move to next point
        F_a = F_b
        k = k + 1
        b = xc*sqrtE(k) - y

        ! Calculate F and H functions
        call calculate_F(F_b, b)
        H = F_a - F_b

        ! Calculate A(k), B(k), and slope terms
        Ak = y_inv_sq*H(2) + TWO*y_inv*H(1) + H(0)
        Bk = y_inv_sq*H(4) + FOUR*y_inv*H(3) + 6.0_8*H(2) + FOUR*y*H(1) &
             + y_sq*H(0)
        slope = (xs(k) - xs(k-1)) / (xc**2*(sqrtE(k)**2 - sqrtE(k-1)**2))

        ! Add contribution to broadened cross section
        sigma = sigma - Ak*(xs(k) - slope*xc**2*sqrtE(k)**2) - slope*Bk
      end do
    end if

  end function sigma1

!===============================================================================
! GAUSS_HERMITE
!===============================================================================

  function gauss_hermite(sqrtE, xs, i_energy, E, xc) result(sigma)
    real(8), intent(in)  :: sqrtE(:)    ! sqrt(energy grid)
    real(8), intent(in)  :: xs(:)       ! unbroadened cross section
    integer, intent(in)  :: i_energy    ! index on energy grid
    real(8), intent(in)  :: E           ! incoming neutron energy
    real(8), intent(in)  :: xc
    real(8)              :: sigma       ! broadened cross section

    integer :: i          ! index on energy grid
    integer :: i_gh       ! index for Gauss-Hermite quadrature
    real(8) :: y          ! beta*v^2
    real(8) :: f          ! inteprolation factor
    real(8) :: xs_interp  ! interpolated cross section

    ! Set target energy based on the index
    y = xc*sqrt(E)

    ! Find index on energy grid just below lowest abscissa
    i = i_energy
    FindEnergy: do while (xc*sqrt(E) - y >= gh_abscissa(1))
      i = i - 1
      if (i < 1) then
        call fatal_error("Gauss-Hermite quadrature used at E=0.")
      end if
    end do FindEnergy

    ! initialize sum for quadrature
    sigma = ZERO

    GaussHermiteQuadrature: do i_gh = 1, N_GAUSS_HERMITE
      ! Find index on energy grid just below current abscissa point
      do while (gh_abscissa(i_gh) >= xc*sqrtE(i+1) - y)
        i = i + 1
      end do

      ! Calculate interpolation factor and interpolated cross section
      f = (gh_abscissa(i_gh) - xc*sqrtE(i) + y)/(xc*(sqrtE(i+1) - sqrtE(i)))
      xs_interp = (ONE - f)*xs(i) + f*xs(i+1)

      ! Add contribution to integral via quadrature
      sigma = sigma + (gh_abscissa(i_gh) + y)**2 * xs_interp * gh_weight(i_gh)
    end do GaussHermiteQuadrature

    ! Multiply by 1/(y^2*sqrt(pi))
    sigma = sigma * sqrt_pi_inv/(y*y)
  end function gauss_hermite

!===============================================================================
! CALCULATE_F evaluates the function:
!
!    F(n,a) = 1/sqrt(pi)*int(z^n*exp(-z^2), z = a to infinity)
!
! The five values returned in a vector correspond to the integral for n = 0
! through 4. These functions are called over and over during the Doppler
! broadening routine.
!===============================================================================

  subroutine calculate_F(F, a)

    real(8), intent(inout) :: F(0:4)
    real(8), intent(in)    :: a

#ifndef NO_F2008
    F(0) = 0.5*erfc(a)
#endif
    F(1) = 0.5*sqrt_pi_inv*exp(-a*a)
    F(2) = 0.5*F(0) + a*F(1)
    F(3) = F(1)*(1.0 + a*a)
    F(4) = 0.75*F(0) + F(1)*a*(1.5 + a*a)

  end subroutine calculate_F

end module doppler

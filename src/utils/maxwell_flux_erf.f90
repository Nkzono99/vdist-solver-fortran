module m_maxwell_flux_erf
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    private
    public :: flux_integral_halfspace_erf, solve_density_from_flux_erf

    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    !> Evaluate the half-space first moment
    !!
    !! \[
    !!   I = \int_0^\infty v f(v)\, dv
    !! \]
    !!
    !! where the 1D shifted Maxwellian is
    !!
    !! \[
    !!   f(v) = \frac{1}{\sqrt{\pi} v_{\rm th}}
    !!          \exp\!\left[-\left(\frac{v - v_0}{v_{\rm th}}\right)^2\right].
    !! \]
    !!
    !! This quantity corresponds to the "one-sided flux moment" and
    !! is widely used for particle flux estimations to a surface.
    !!
    !! The analytical closed-form expression is
    !!
    !! \[
    !!   I = \frac{1}{2} v_0 \left[ 1 + {\rm erf}\!\left(\frac{v_0}{v_{\rm th}}\right) \right]
    !!       + \frac{v_{\rm th}}{2\sqrt{\pi}}
    !!         \exp\!\left[-\left(\frac{v_0}{v_{\rm th}}\right)^2\right].
    !! \]
    !!
    !! @param[in] vth  Thermal speed \(v_{\rm th}\)
    !! @param[in] v0   Drift speed \(v_0\)
    !! @return    I    Half-space flux moment
    pure function flux_integral_halfspace_erf(vth, v0) result(I)
        real(dp), intent(in) :: vth
        real(dp), intent(in) :: v0
        real(dp) :: I
        real(dp) :: a, expa

        a = v0/vth
        expa = exp(-a*a)

        ! Analytical formula using erf
        I = 0.5_dp*v0*(1.0_dp + erf(a)) &
            + 0.5_dp*vth/sqrt(pi)*expa
    end function flux_integral_halfspace_erf

    !> Compute number density \(n\) from the relation
    !!
    !! \[
    !!    $j = n I, \quad I = \int_0^\infty v f(v)\,dv.$
    !! \]
    !!
    !! Therefore
    !!
    !! \[
    !!    n = \frac{j}{I}.
    !! \]
    !!
    !! @param[in] j   Flux-like quantity (particle flux; if current, divide by charge)
    !! @param[in] vth Thermal speed
    !! @param[in] v0  Drift speed
    !! @return    n   Number density
    pure function solve_density_from_flux_erf(j, vth, v0) result(n)
        real(dp), intent(in) :: j
        real(dp), intent(in) :: vth, v0
        real(dp) :: n
        real(dp) :: I

        I = flux_integral_halfspace_erf(vth, v0)

        if (I > 0.0_dp) then
            n = j/I
        else
            n = 0.0_dp
        end if
    end function solve_density_from_flux_erf

end module

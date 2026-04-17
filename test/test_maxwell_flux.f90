program test_maxwell_flux
    use m_maxwell_flux_erf, only: flux_integral_halfspace_erf, &
                                  solve_density_from_flux_erf
    use m_test_helpers, only: assert_close

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)

    call test_halfspace_integral_zero_drift()
    call test_density_at_zero_drift()
    call test_density_large_drift_limit()
    call test_density_negative_drift_returns_zero()
    call test_density_for_zero_flux_is_zero()

    print *, "test_maxwell_flux: all tests passed."

contains

    subroutine test_halfspace_integral_zero_drift()
        !! With v0 = 0, analytic result is I = vth / (2*sqrt(pi)).
        double precision :: I, expected

        I = flux_integral_halfspace_erf(1.0d0, 0.0d0)
        expected = 1.0d0/(2d0*sqrt(pi))

        call assert_close("flux_integral (vth=1, v0=0) = 1/(2*sqrt(pi))", I, expected)
    end subroutine

    subroutine test_density_at_zero_drift()
        !! j = n * I, so n = j * 2 * sqrt(pi) / vth at v0 = 0.
        double precision :: n, expected

        n = solve_density_from_flux_erf(1.0d0, 1.0d0, 0.0d0)
        expected = 2d0*sqrt(pi)

        call assert_close("density for j=1, vth=1, v0=0", n, expected)
    end subroutine

    subroutine test_density_large_drift_limit()
        !! For v0 >> vth the integral approaches v0, so n -> j / v0.
        double precision :: v0, vth, j_flux, n

        vth = 1.0d0
        v0 = 100d0
        j_flux = 50d0

        n = solve_density_from_flux_erf(j_flux, vth, v0)

        call assert_close("density for strongly-drifting beam approaches j/v0", &
                          n, j_flux/v0, tolerance=1d-3)
    end subroutine

    subroutine test_density_negative_drift_returns_zero()
        !! Strongly negative drift makes the half-space integral vanish;
        !! the routine should fall back to n = 0 rather than blow up.
        double precision :: n

        n = solve_density_from_flux_erf(1.0d0, 1.0d0, -50d0)

        call assert_close("density with large negative drift is 0", n, 0d0)
    end subroutine

    subroutine test_density_for_zero_flux_is_zero()
        double precision :: n

        n = solve_density_from_flux_erf(0.0d0, 1.0d0, 0.5d0)

        call assert_close("density when flux=0 is 0", n, 0d0)
    end subroutine

end program

program test_probabilities
    use m_probabilities, only: t_ZeroProbability, new_ZeroProbability, &
                               t_MaxwellianProbability, new_MaxwellianProbability, &
                               tp_Probability, t_Probability
    use m_test_helpers, only: assert_close

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)

    call test_zero_probability_is_zero()
    call test_maxwellian_peak_value()
    call test_maxwellian_offset_symmetry()
    call test_maxwellian_coefficient_scales_output()
    call test_tp_probability_forwards_to_ref()

    print *, "test_probabilities: all tests passed."

contains

    subroutine test_zero_probability_is_zero()
        type(t_ZeroProbability) :: prob

        prob = new_ZeroProbability()

        call assert_close("ZeroProbability at origin", &
                          prob%at([0d0, 0d0, 0d0], [0d0, 0d0, 0d0]), 0d0)
        call assert_close("ZeroProbability at arbitrary point", &
                          prob%at([3d0, -2d0, 7d0], [1d0, 1d0, 1d0]), 0d0)
    end subroutine

    subroutine test_maxwellian_peak_value()
        !! Peak of the 3D Maxwellian (all velocities equal to the mean) is
        !! (1 / sqrt(2*pi))**3 when all scales are 1 and coefficient = 1.
        type(t_MaxwellianProbability) :: prob
        double precision :: expected

        prob = new_MaxwellianProbability([0d0, 0d0, 0d0], [1d0, 1d0, 1d0])

        expected = (1d0/sqrt(2d0*pi))**3

        call assert_close("Maxwellian peak = (2*pi)^(-3/2)", &
                          prob%at([0d0, 0d0, 0d0], [0d0, 0d0, 0d0]), expected)
    end subroutine

    subroutine test_maxwellian_offset_symmetry()
        !! The 1-sigma offset on a single axis should scale the peak by
        !! exp(-1/2); applied across three axes (1-sigma each) it becomes
        !! exp(-3/2). The spatial position argument is not consulted.
        type(t_MaxwellianProbability) :: prob
        double precision :: peak, offset_value, expected_ratio

        prob = new_MaxwellianProbability([0d0, 0d0, 0d0], [1d0, 1d0, 1d0])

        peak = prob%at([0d0, 0d0, 0d0], [0d0, 0d0, 0d0])
        offset_value = prob%at([0d0, 0d0, 0d0], [1d0, 1d0, 1d0])

        expected_ratio = exp(-1.5d0)

        call assert_close("Maxwellian drops by exp(-3/2) at 1-sigma corner", &
                          offset_value/peak, expected_ratio)
    end subroutine

    subroutine test_maxwellian_coefficient_scales_output()
        type(t_MaxwellianProbability) :: prob_a, prob_b
        double precision :: val_a, val_b

        prob_a = new_MaxwellianProbability([0d0, 0d0, 0d0], [1d0, 1d0, 1d0])
        prob_b = new_MaxwellianProbability([0d0, 0d0, 0d0], [1d0, 1d0, 1d0], coefficient=7d0)

        val_a = prob_a%at([0d0, 0d0, 0d0], [0.3d0, -0.1d0, 0.2d0])
        val_b = prob_b%at([0d0, 0d0, 0d0], [0.3d0, -0.1d0, 0.2d0])

        call assert_close("Maxwellian coefficient multiplies the PDF", val_b, 7d0*val_a)
    end subroutine

    subroutine test_tp_probability_forwards_to_ref()
        type(tp_Probability) :: wrapper

        allocate (wrapper%ref, source=new_ZeroProbability())

        call assert_close("tp_Probability forwards to referenced ZeroProbability", &
                          wrapper%at([1d0, 2d0, 3d0], [4d0, 5d0, 6d0]), 0d0)
    end subroutine

end program

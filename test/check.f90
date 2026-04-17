program check
    use m_dust_charge_simulator
    use m_field

    implicit none

    call test_negative_potential_electron_current()
    call test_negative_potential_ion_current()

    print *, "All tests passed."

contains

    subroutine test_negative_potential_electron_current()
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: simulator
        type(t_DustParticle) :: dust
        type(t_DustParticle) :: dust_new

        double precision :: values(3, 0:1, 0:1, 0:1)
        double precision :: expected_charge
        double precision, parameter :: pi = acos(-1.0d0)

        values(:, :, :, :) = 0d0
        values(1, :, :, :) = 1d0

        currents = new_VectorFieldGrid(3, 1, 1, 1, values)
        simulator = new_DustChargeSimulator(1, 1, 1, 1, [2d0], currents)

        dust = new_DustParticle(-4d0*pi, 1d0, 1d0, [0.5d0, 0.5d0, 0.5d0], [0d0, 0d0, 0d0], 0d0)
        dust_new = simulator%update(dust, 1d0)

        expected_charge = dust%charge + (-exp(-0.5d0))*4d0*pi*(-1d0)

        call assert_close("electron current for negative dust potential", dust_new%charge, expected_charge)
    end subroutine

    subroutine test_negative_potential_ion_current()
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: simulator
        type(t_DustParticle) :: dust
        type(t_DustParticle) :: dust_new

        double precision :: values(6, 0:1, 0:1, 0:1)
        double precision :: expected_charge
        double precision, parameter :: pi = acos(-1.0d0)

        values(:, :, :, :) = 0d0
        values(4, :, :, :) = 1d0

        currents = new_VectorFieldGrid(6, 1, 1, 1, values)
        simulator = new_DustChargeSimulator(1, 1, 1, 2, [2d0, 2d0], currents)

        dust = new_DustParticle(-4d0*pi, 1d0, 1d0, [0.5d0, 0.5d0, 0.5d0], [0d0, 0d0, 0d0], 0d0)
        dust_new = simulator%update(dust, 1d0)

        expected_charge = dust%charge + (1d0 - (-1d0)/2d0)*4d0*pi*(-1d0)

        call assert_close("ion current for negative dust potential", dust_new%charge, expected_charge)
    end subroutine

    subroutine assert_close(label, actual, expected)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: actual
        double precision, intent(in) :: expected

        double precision, parameter :: tolerance = 1d-12

        if (abs(actual - expected) > tolerance) then
            print *, "FAILED:", trim(label)
            print *, "  actual  =", actual
            print *, "  expected=", expected
            error stop 1
        end if
    end subroutine

end program check

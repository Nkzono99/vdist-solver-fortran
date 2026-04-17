program test_dust_charging
    use m_dust_charge_simulator, only: t_DustParticle, new_DustParticle, &
                                       t_DustChargeSimulator, new_DustChargeSimulator
    use m_field, only: t_VectorFieldGrid, new_VectorFieldGrid
    use m_test_helpers, only: assert_close

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)

    call test_dust_potential_formula()
    call test_update_is_noop_when_dt_is_zero()

    call test_negative_potential_electron_current()
    call test_negative_potential_ion_current()

    call test_positive_potential_electron_current()
    call test_positive_potential_ion_current()

    call test_photoelectron_current_positive_potential()
    call test_photoelectron_current_negative_potential()

    print *, "test_dust_charging: all tests passed."

contains

    function make_currents(values) result(currents)
        !! Helper for a uniform (1,1,1)-cell VectorFieldGrid whose shape is
        !! (n_elements, 0:1, 0:1, 0:1) and all corners share the same vector.
        double precision, intent(in) :: values(:)
        type(t_VectorFieldGrid) :: currents

        double precision :: buf(size(values), 0:1, 0:1, 0:1)
        integer :: i, ix, iy, iz

        do iz = 0, 1; do iy = 0, 1; do ix = 0, 1
                    do i = 1, size(values)
                        buf(i, ix, iy, iz) = values(i)
                    end do
                end do; end do; end do

        currents = new_VectorFieldGrid(size(values), 1, 1, 1, buf)
    end function

    function make_dust(charge) result(dust)
        double precision, intent(in) :: charge
        type(t_DustParticle) :: dust

        ! radius = 1, mass = 1, at the centre of the unit cell.
        dust = new_DustParticle(charge, 1d0, 1d0, &
                                [0.5d0, 0.5d0, 0.5d0], &
                                [0d0, 0d0, 0d0], &
                                0d0)
    end function

    subroutine test_dust_potential_formula()
        type(t_DustParticle) :: dust

        dust = make_dust(4d0*pi)  ! charge = 4*pi, radius = 1 -> potential = 1

        call assert_close("potential = charge / (4 pi r)", dust%potential(), 1d0)
    end subroutine

    subroutine test_update_is_noop_when_dt_is_zero()
        !! Gravity term scales with dt too, so dt = 0 must leave charge and
        !! velocity untouched.
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust, dust_new

        currents = make_currents([1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 1, [2d0], currents)

        dust = make_dust(-4d0*pi)
        dust_new = sim%update(dust, 0d0)

        call assert_close("dt=0 keeps charge", dust_new%charge, dust%charge)
        call assert_close("dt=0 keeps vz", dust_new%particle%velocity(3), 0d0)
    end subroutine

    subroutine test_negative_potential_electron_current()
        !! phid = -1, Te = 2: electron current  = -je0 * exp(phid/Te).
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust_new

        currents = make_currents([1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 1, [2d0], currents)

        dust_new = sim%update(make_dust(-4d0*pi), 1d0)

        call assert_close("electron current for phid<0", &
                          dust_new%charge, -4d0*pi + exp(-0.5d0)*4d0*pi)
    end subroutine

    subroutine test_negative_potential_ion_current()
        !! phid = -1, Ti = 2, species 2 only: ion current = ji0 * (1 - phid/Ti).
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust_new

        currents = make_currents([0d0, 0d0, 0d0, 1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 2, [2d0, 2d0], currents)

        dust_new = sim%update(make_dust(-4d0*pi), 1d0)

        call assert_close("ion current for phid<0", &
                          dust_new%charge, -4d0*pi + (1d0 - (-1d0)/2d0)*4d0*pi*(-1d0))
    end subroutine

    subroutine test_positive_potential_electron_current()
        !! phid = +1, Te = 2: electron current = -je0 * (1 + phid/Te).
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust_new

        currents = make_currents([1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 1, [2d0], currents)

        dust_new = sim%update(make_dust(4d0*pi), 1d0)

        ! net_current = -1.5 ; delta = -net_current * 4*pi * (-(-1)) = +1.5 * 4*pi
        call assert_close("electron current for phid>0", &
                          dust_new%charge, 4d0*pi + 1.5d0*4d0*pi)
    end subroutine

    subroutine test_positive_potential_ion_current()
        !! phid = +1, Ti = 2, species 2: ion current = ji0 * exp(-phid/Ti).
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust_new
        double precision :: ion_current

        currents = make_currents([0d0, 0d0, 0d0, 1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 2, [2d0, 2d0], currents)

        dust_new = sim%update(make_dust(4d0*pi), 1d0)

        ion_current = exp(-0.5d0)  ! ji0 = 1
        ! delta = net_current * 4*pi * (-(-1)) = ion_current * 4*pi * (-1)
        call assert_close("ion current for phid>0", &
                          dust_new%charge, 4d0*pi + ion_current*4d0*pi*(-1d0))
    end subroutine

    subroutine test_photoelectron_current_positive_potential()
        !! phid = +1, Tph = 2, species 3 only, jph0 = 0.5, j0 = 1.
        !! ret = jph0 * exp(-0.5) * 1.5 - 1 * 1.5
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust_new
        double precision :: net_current

        currents = make_currents([0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 3, [2d0, 2d0, 2d0], currents, jph0=0.5d0)

        dust_new = sim%update(make_dust(4d0*pi), 1d0)

        net_current = 0.5d0*exp(-0.5d0)*1.5d0 - 1.5d0
        call assert_close("photoelectron current for phid>0", &
                          dust_new%charge, 4d0*pi + net_current*4d0*pi*(-1d0))
    end subroutine

    subroutine test_photoelectron_current_negative_potential()
        !! phid = -1, Tph = 2, jph0 = 0.5, j0 = 1.
        !! ret = jph0 - 1 * (1 + phid/Tph) = 0.5 - 0.5 = 0
        type(t_VectorFieldGrid) :: currents
        type(t_DustChargeSimulator) :: sim
        type(t_DustParticle) :: dust_new
        double precision :: net_current

        currents = make_currents([0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0])
        sim = new_DustChargeSimulator(1, 1, 1, 3, [2d0, 2d0, 2d0], currents, jph0=0.5d0)

        dust_new = sim%update(make_dust(-4d0*pi), 1d0)

        net_current = 0.5d0 - 1d0*(1d0 + (-1d0)/2d0)
        call assert_close("photoelectron current for phid<0", &
                          dust_new%charge, -4d0*pi + net_current*4d0*pi*(-1d0))
    end subroutine

end program

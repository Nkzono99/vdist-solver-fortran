program test_emission_normalization
    use, intrinsic :: iso_c_binding, only: c_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_simulator, only: t_ESSimulator
    use m_probabilities, only: t_MaxwellianProbability
    use m_test_helpers, only: assert_close, assert_true
    implicit none
    integer :: axis, direction, drift_case
    double precision :: drifts(3)

    drifts = [0d0, 0.6d0, -0.4d0]
    do axis = 1, 3
        do direction = -1, 1, 2
            do drift_case = 1, size(drifts)
                call check_flux(direction*axis, drifts(drift_case), drift_case == 2)
            end do
        end do
    end do
    print *, 'test_emission_normalization: all tests passed.'
contains
    subroutine check_flux(normal, drift, override)
        integer, intent(in) :: normal
        double precision, intent(in) :: drift
        logical, intent(in) :: override
        character(len=*), parameter :: filename = 'test-emission-normalization.inp'
        character(kind=c_char) :: cpath(len(filename))
        type(t_ESSimulator) :: simulator
        double precision :: fields(9, 3, 3, 3), velocity(3), flux, expected, patch_current
        double precision :: ds, dtan, vn, tangent_integral, normal_integral, peak, pi
        integer :: unit, i, a, b, c, tag

        expected = 1.7d0
        patch_current = -1d0
        if (override) then
            expected = 0.37d0
            patch_current = expected
        end if
        open(newunit=unit, file=filename, status='replace')
        write(unit, *) '&esorem emflag=0 /'
        write(unit, *) '&plasma wp=10*0, wc=0, phiz=0, phixy=0 /'
        write(unit, *) '&tmgrid nx=2, ny=2, nz=2, dt=0.1 /'
        write(unit, *) '&system nspec=1, npbnd=30*2 /'
        write(unit, *) '&intp qm=10*1, peth=0.7, path=1.3, spa=', drift, &
                       ', spe=0, speth=0, vdri=10*0, vdthz=10*0, vdthxy=10*0 /'
        write(unit, *) "&ptcond boundary_type='none' /"
        write(unit, *) '&emissn use_raycast=.false., nflag_emit=1, nepl=1, nemd=', normal, &
                       ', curf=1.7, curfs=', patch_current, &
                       ', xmine=1, xmaxe=1, ymine=1, ymaxe=1, zmine=1, zmaxe=1 /'
        close(unit)
        do i = 1, len(filename)
            cpath(i) = filename(i:i)
        end do
        fields = 0d0
        simulator = create_simulator(cpath, len(filename), 2, 2, 2, fields, 1, 10)
        open(newunit=unit, file=filename, status='old')
        close(unit, status='delete')
        tag = 0
        do i = 1, size(simulator%probability_functions)
            if (associated(simulator%probability_functions(i)%ref)) tag = i
        end do
        a = abs(normal)
        b = mod(a, 3) + 1
        c = mod(a + 1, 3) + 1
        select type (p => simulator%probability_functions(tag)%ref)
        type is (t_MaxwellianProbability)
            call assert_close('normal sigma is path', p%scales(a), 1.3d0)
            call assert_close('tangential sigma is peth', p%scales(b), 0.7d0)
            ! Numerically integrate each separable factor of the actual built
            ! 3D distribution. No analytical flux helper is used in this oracle.
            velocity = p%locs
            peak = p%at([1d0, 1d0, 1d0], velocity)
            ds = (abs(drift) + 12d0*1.3d0)/20000
            normal_integral = 0d0
            do i = 1, 20000
                vn = (i - 0.5d0)*ds
                velocity(a) = sign(1d0, dble(normal))*vn
                normal_integral = normal_integral + vn*p%at([1d0, 1d0, 1d0], velocity)*ds
            end do
            velocity = p%locs
            dtan = 24d0*0.7d0/2000
            tangent_integral = 0d0
            do i = 1, 2000
                velocity(b) = p%locs(b) + (i - 1000.5d0)*dtan
                tangent_integral = tangent_integral + p%at([1d0, 1d0, 1d0], velocity)*dtan/peak
            end do
            flux = normal_integral*tangent_integral**2
            call assert_true('one-sided emission flux is finite', ieee_is_finite(flux))
            call assert_close('one-sided emission flux recovers curf/curfs', flux, expected, 1d-7)
            if (drift == 0d0) then
                pi = acos(-1d0)
                call assert_close('zero drift uses Gaussian sigma normalization', &
                                  p%coefficient, expected*sqrt(2d0*pi)/1.3d0)
            end if
        class default
            error stop 'emission source must be Maxwellian'
        end select
        call destroy_simulator(simulator)
    end subroutine
end program

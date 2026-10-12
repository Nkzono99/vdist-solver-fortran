program test_emses_thermal
    use, intrinsic :: iso_c_binding, only: c_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_simulator, only: t_ESSimulator
    use m_allcom, only: vth_vector
    use m_test_helpers, only: assert_close, assert_close_vec, assert_true

    implicit none
    double precision, parameter :: pi = acos(-1d0)

    call test_isotropic_boundary_is_rotation_invariant()
    call test_oblique_anisotropic_boundary_pdf_and_moments()
    call test_photoelectron_uses_rotated_thermal_distribution()
    print *, "test_emses_thermal: all tests passed."

contains

    function build_source(theta, phi, perpendicular, parallel, photoelectron) result(simulator)
        double precision, intent(in) :: theta, phi, perpendicular, parallel
        logical, intent(in), optional :: photoelectron
        type(t_ESSimulator) :: simulator
        character(len=*), parameter :: filename = 'test-emses-thermal.inp'
        character(kind=c_char) :: cpath(len(filename))
        double precision :: fields(9, 3, 3, 3)
        integer :: unit, i, emit
        logical :: raycast

        emit = 0
        raycast = .false.
        if (present(photoelectron)) then
            if (photoelectron) then
                emit = 2
                raycast = .true.
            end if
        end if
        open(newunit=unit, file=filename, status='replace')
        write(unit, *) '&esorem emflag=0 /'
        write(unit, *) '&plasma wp=10*0, wc=0, phiz=', theta, ', phixy=', phi, ' /'
        write(unit, *) '&tmgrid nx=2, ny=2, nz=2, dt=0.1 /'
        write(unit, *) '&system nspec=1, npbnd=30*2 /'
        write(unit, *) '&intp qm=10*1, spa=10*0, spe=10*0, speth=10*0, ', &
                       'vdri=10*0, vdthz=10*0, vdthxy=10*0, peth=', perpendicular, ', path=', parallel, ' /'
        write(unit, *) "&ptcond boundary_type='none' /"
        write(unit, *) '&emissn nepl=10*0, nflag_emit=', emit, ', use_raycast=', raycast, &
                       ', ray_zenith_angle_deg=180, ray_azimuth_angle_deg=0 /'
        close(unit)
        do i = 1, len(filename)
            cpath(i) = filename(i:i)
        end do
        fields = 0d0
        simulator = create_simulator(cpath, len(filename), 2, 2, 2, fields, 1, 10)
        open(newunit=unit, file=filename, status='old')
        close(unit, status='delete')
    end function

    subroutine test_isotropic_boundary_is_rotation_invariant()
        type(t_ESSimulator) :: simulator
        double precision :: angles(4), velocity(3), expected, actual
        integer :: i

        angles = [0d0, 30d0, 45d0, 90d0]
        velocity = [1d0, -0.3d0, 0.2d0]
        expected = exp(-sum(velocity**2)/2d0)/(2d0*pi)**1.5d0
        do i = 1, size(angles)
            simulator = build_source(angles(i), 31d0, 1d0, 1d0)
            actual = simulator%probability_functions(3)%at([0d0, 1d0, 1d0], velocity)
            call assert_true("rotated isotropic density is finite", ieee_is_finite(actual))
            call assert_close("isotropic boundary density is rotation invariant", actual, expected)
            call assert_close_vec("isotropic autorange widths remain one", vth_vector(1), [1d0, 1d0, 1d0])
            call destroy_simulator(simulator)
        end do
    end subroutine

    subroutine test_oblique_anisotropic_boundary_pdf_and_moments()
        type(t_ESSimulator) :: simulator
        double precision :: b(3), velocity(3), perpendicular, parallel, expected, actual
        double precision :: norm, covariance(3, 3), expected_covariance(3, 3), mean(3), weight
        double precision, parameter :: spacing = 0.25d0
        integer :: ix, iy, iz, i, j

        perpendicular = 0.5d0
        parallel = 1.4d0
        b = [-sin(pi/4d0)*cos(31d0*pi/180d0), sin(pi/4d0)*sin(31d0*pi/180d0), cos(pi/4d0)]
        simulator = build_source(45d0, 31d0, perpendicular, parallel)
        velocity = [0.7d0, -0.2d0, 0.4d0]
        expected = exp(-0.5d0*(sum(velocity**2)/perpendicular**2 &
                   + (1d0/parallel**2 - 1d0/perpendicular**2)*dot_product(velocity, b)**2)) &
                   /((2d0*pi)**1.5d0*perpendicular**2*parallel)
        actual = simulator%probability_functions(3)%at([0d0, 1d0, 1d0], velocity)
        call assert_close("anisotropic density includes oblique covariance", actual, expected)

        norm = 0d0
        mean = 0d0
        covariance = 0d0
        do iz = 1, 64
            do iy = 1, 64
                do ix = 1, 64
                    velocity = (dble([ix, iy, iz]) - 32.5d0)*spacing
                    weight = simulator%probability_functions(3)%at([0d0, 1d0, 1d0], velocity)*spacing**3
                    norm = norm + weight
                    mean = mean + weight*velocity
                    do j = 1, 3
                        do i = 1, 3
                            covariance(i, j) = covariance(i, j) + weight*velocity(i)*velocity(j)
                        end do
                    end do
                end do
            end do
        end do
        call assert_close("rotated density integrates to one", norm, 1d0, 1d-8)
        call assert_close_vec("rotated density has zero mean", mean, [0d0, 0d0, 0d0], 1d-8)
        do j = 1, 3
            expected_covariance(:, j) = (parallel**2 - perpendicular**2)*b*b(j)
            expected_covariance(j, j) = expected_covariance(j, j) + perpendicular**2
            call assert_close_vec("rotated covariance matches analytic tensor", covariance(:, j), &
                                  expected_covariance(:, j), 1d-7)
        end do
        call assert_close_vec("autorange uses marginal standard deviations", vth_vector(1), &
            sqrt([expected_covariance(1, 1), expected_covariance(2, 2), expected_covariance(3, 3)]))
        call destroy_simulator(simulator)
    end subroutine

    subroutine test_photoelectron_uses_rotated_thermal_distribution()
        type(t_ESSimulator) :: simulator
        double precision :: actual, expected

        simulator = build_source(45d0, 31d0, 1d0, 1d0, photoelectron=.true.)
        actual = simulator%probability_functions(1)%at([1d0, 1d0, 1d0], [0d0, 0d0, 1d0])
        expected = 2d0*exp(-0.5d0)/(2d0*pi)**1.5d0
        call assert_true("rotated photoelectron density is finite", ieee_is_finite(actual))
        call assert_close("photoelectron builder preserves isotropic density", actual, expected)
        call destroy_simulator(simulator)
    end subroutine
end program

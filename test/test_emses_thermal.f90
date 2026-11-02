program test_emses_thermal
    use, intrinsic :: iso_c_binding, only: c_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_simulator, only: t_ESSimulator
    use m_allcom, only: vth_vector, thermal_basis, vdri_vector, spa, spe, speth, vdri, vdthz, vdthxy
    use m_test_helpers, only: assert_close, assert_close_vec, assert_true

    implicit none
    double precision, parameter :: pi = acos(-1d0)

    call test_magnetic_frame_matches_emses_background_field()
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

    subroutine test_magnetic_frame_matches_emses_background_field()
        !! Issue #42: the parallel axis of the native source frame must be the
        !! background field direction that EMSES and the Python wrapper use,
        !! B0hat = (sin(phiz)cos(phixy), sin(phiz)sin(phixy), cos(phiz)).
        type(t_ESSimulator) :: simulator
        double precision :: angles(4, 2), b(3), basis(3, 3), drift(3), theta, phi
        integer :: i, j, k

        angles(:, 1) = [0d0, 30d0, 45d0, 60d0]
        angles(:, 2) = [0d0, 31d0, 0d0, 120d0]
        do i = 1, size(angles, 1)
            theta = angles(i, 1)*pi/180d0
            phi = angles(i, 2)*pi/180d0
            b = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
            simulator = build_source(angles(i, 1), angles(i, 2), 1d0, 1d0)
            basis = thermal_basis()
            call assert_close_vec("parallel axis equals the EMSES background field direction", basis(:, 3), b)
            do j = 1, 3
                do k = 1, 3
                    call assert_close("thermal basis is orthonormal", dot_product(basis(:, j), basis(:, k)), &
                                      merge(1d0, 0d0, j == k))
                end do
            end do
            call assert_close("thermal basis is right-handed", &
                              dot_product(basis(:, 3), cross3(basis(:, 1), basis(:, 2))), 1d0)

            ! spa is a drift along B0; spe/speth lie in the perpendicular plane (speth in degrees).
            spa(1) = 3d0
            spe(1) = 0d0
            speth(1) = 0d0
            drift = vdri_vector(1)
            call assert_close_vec("spa drifts along the packed background field", drift, 3d0*b)
            spa(1) = 0d0
            spe(1) = 2d0
            speth(1) = 90d0
            drift = vdri_vector(1)
            call assert_close("spe drift is perpendicular to B0", dot_product(drift, b), 0d0)
            call assert_close_vec("speth=90 deg selects the second perpendicular axis", drift, 2d0*basis(:, 2))
            spe(1) = 0d0
            speth(1) = 0d0

            ! vdri uses zenith vdthz from +z and azimuth vdthxy from +x.
            vdri(1) = 1.5d0
            vdthz(1) = 90d0
            vdthxy(1) = 90d0
            call assert_close_vec("vdri with zenith 90 and azimuth 90 points along +y", vdri_vector(1), &
                                  [0d0, 1.5d0, 0d0])
            vdthz(1) = angles(i, 1)
            vdthxy(1) = angles(i, 2)
            call assert_close_vec("vdri at the field angles is parallel to B0", vdri_vector(1), 1.5d0*b)
            vdri(1) = 0d0
            vdthz(1) = 0d0
            vdthxy(1) = 0d0
            call destroy_simulator(simulator)
        end do
    end subroutine

    pure function cross3(a, b) result(c)
        double precision, intent(in) :: a(3), b(3)
        double precision :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
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
        ! EMSES background-field direction for phiz=45, phixy=31.
        b = [sin(pi/4d0)*cos(31d0*pi/180d0), sin(pi/4d0)*sin(31d0*pi/180d0), cos(pi/4d0)]
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

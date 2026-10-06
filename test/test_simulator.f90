program test_simulator
    !! Physics regressions for the Boris push and its inverse.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use finbound, only: t_BoundaryList, new_BoundaryList, t_CollisionRecord
    use m_particle, only: t_Particle, new_Particle
    use m_field, only: new_VectorFieldGrid
    use m_probabilities, only: tp_Probability
    use m_simulator, only: t_ESSimulator, new_ESSimulator
    use m_test_helpers, only: assert_close, assert_close_vec, assert_true

    implicit none

    call test_axis_aligned_gyro_sign()
    call test_oblique_magnetic_invariants()
    call test_linear_electric_inverse()
    call test_linear_magnetic_inverse()
    call test_nonuniform_mixed_fields_roundtrip()
    call test_periodic_crossing_roundtrip()
    call test_zero_step()

    print *, "test_simulator: all tests passed."

contains

    function build_simulator(values, periodic) result(simulator)
        double precision, intent(in) :: values(6, 0:4, 0:4, 0:4)
        logical, intent(in), optional :: periodic
        type(t_ESSimulator) :: simulator
        type(t_BoundaryList) :: boundaries
        type(tp_Probability) :: probabilities(0)
        integer :: boundary_conditions(3)

        boundary_conditions = [2, 2, 2]
        if (present(periodic)) then
            if (periodic) boundary_conditions = [0, 0, 0]
        end if
        boundaries = new_BoundaryList()
        simulator = new_ESSimulator(4, 4, 4, boundary_conditions, &
                                   new_VectorFieldGrid(6, 4, 4, 4, values), boundaries, probabilities)
    end function

    function advance(simulator, particle, dt) result(ret)
        type(t_ESSimulator), intent(in) :: simulator
        type(t_Particle), intent(in) :: particle
        double precision, intent(in) :: dt
        type(t_Particle) :: ret
        type(t_CollisionRecord) :: record

        ret = simulator%update(particle, dt, record)
        call assert_true("empty boundary list has no collision", .not. record%is_collided)
        call assert_true("updated velocity is finite", all(ieee_is_finite(ret%velocity)))
    end function

    subroutine test_axis_aligned_gyro_sign()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: p, ret
        double precision :: values(6, 0:4, 0:4, 0:4), u, expected(3)
        integer :: sign_q, sign_dt

        values = 0d0
        values(6, :, :, :) = 1d0
        simulator = build_simulator(values)
        do sign_q = -1, 1, 2
            p = new_Particle(dble(sign_q), [2d0, 2d0, 2d0], [1d0, 0d0, 0d0])
            do sign_dt = -1, 1, 2
                u = dble(sign_q*sign_dt)*0.1d0
                expected = [(1d0 - u*u)/(1d0 + u*u), 2d0*u/(1d0 + u*u), 0d0]
                ret = advance(simulator, p, dble(sign_dt)*0.2d0)
                call assert_close_vec("gyro sign for either charge and time direction", ret%velocity, expected)
            end do
        end do
    end subroutine

    subroutine test_oblique_magnetic_invariants()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: p, initial
        double precision :: values(6, 0:4, 0:4, 0:4), b(3), speed_error, parallel_error
        integer :: sign_q, sign_dt, istep

        b = [1d0, 2d0, 3d0]
        values = 0d0
        values(4, :, :, :) = b(1)
        values(5, :, :, :) = b(2)
        values(6, :, :, :) = b(3)
        simulator = build_simulator(values, periodic=.true.)
        do sign_q = -1, 1, 2
            do sign_dt = -1, 1, 2
                initial = new_Particle(dble(sign_q), [2d0, 2d0, 2d0], [1d0, 0d0, 0d0])
                p = initial
                speed_error = 0d0
                parallel_error = 0d0
                do istep = 1, 1000
                    p = advance(simulator, p, dble(sign_dt)*0.2d0)
                    speed_error = max(speed_error, abs(norm2(p%velocity) - norm2(initial%velocity)))
                    parallel_error = max(parallel_error, &
                                         abs(dot_product(p%velocity - initial%velocity, b)/norm2(b)))
                end do
                call assert_close("oblique B preserves speed for 1000 steps", speed_error, 0d0, tolerance=1d-11)
                call assert_close("oblique B preserves parallel velocity", parallel_error, 0d0, tolerance=1d-11)
            end do
        end do
    end subroutine

    subroutine test_linear_electric_inverse()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: p, forward, recovered
        double precision :: values(6, 0:4, 0:4, 0:4)
        integer :: i

        values = 0d0
        do i = 0, 4
            values(1, i, :, :) = dble(i)
        end do
        simulator = build_simulator(values)
        p = new_Particle(1d0, [0.5d0, 2d0, 2d0], [0.5d0, 0d0, 0d0])
        forward = advance(simulator, p, -0.1d0)
        call assert_close("Ex=x forward velocity", forward%velocity(1), 0.55d0)
        call assert_close("forward drift uses updated velocity", forward%position(1), 0.555d0)
        recovered = advance(simulator, forward, 0.1d0)
        call assert_close_vec("Ex=x inverse position", recovered%position, p%position)
        call assert_close_vec("Ex=x inverse velocity", recovered%velocity, p%velocity)
    end subroutine

    subroutine test_linear_magnetic_inverse()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: p, forward, recovered
        double precision :: values(6, 0:4, 0:4, 0:4)
        integer :: i

        values = 0d0
        do i = 0, 4
            values(6, i, :, :) = dble(i)
        end do
        simulator = build_simulator(values)
        p = new_Particle(1d0, [0.5d0, 2d0, 2d0], [0.5d0, 0d0, 0d0])
        forward = advance(simulator, p, -0.2d0)
        recovered = advance(simulator, forward, 0.2d0)
        call assert_close_vec("Bz=x inverse position", recovered%position, p%position)
        call assert_close_vec("Bz=x inverse velocity", recovered%velocity, p%velocity)
    end subroutine

    subroutine test_nonuniform_mixed_fields_roundtrip()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: initial, p
        double precision :: values(6, 0:4, 0:4, 0:4), dt
        integer :: i, j, k, sign_q, sign_dt, istep

        do k = 0, 4
            do j = 0, 4
                do i = 0, 4
                    values(:, i, j, k) = [0.3d0*i + 0.1d0*j, 0.1d0*i + 0.2d0*j, -0.1d0*k, &
                                         0.4d0 + 0.1d0*j, -0.3d0 + 0.05d0*k, 0.2d0 + 0.1d0*i]
                end do
            end do
        end do
        simulator = build_simulator(values)
        do sign_q = -1, 1, 2
            do sign_dt = -1, 1, 2
                initial = new_Particle(0.7d0*sign_q, [2d0, 2d0, 2d0], [0.15d0, -0.2d0, 0.1d0])
                p = initial
                dt = 0.1d0*sign_dt
                do istep = 1, 10
                    p = advance(simulator, p, dt)
                end do
                do istep = 1, 10
                    p = advance(simulator, p, -dt)
                end do
                call assert_close_vec("mixed nonuniform fields roundtrip position", p%position, initial%position)
                call assert_close_vec("mixed nonuniform fields roundtrip velocity", p%velocity, initial%velocity)
                call assert_close("roundtrip trace time", p%t, initial%t)
                call assert_close("roundtrip charge-to-mass ratio", p%q_m, initial%q_m)
            end do
        end do
    end subroutine

    subroutine test_periodic_crossing_roundtrip()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: p, forward, recovered
        double precision :: values(6, 0:4, 0:4, 0:4), ex(0:4), bz(0:4)
        integer :: i, direction

        ex = [0d0, 0.4d0, -0.1d0, -0.3d0, 0d0]
        bz = [0.2d0, 1d0, 0.7d0, 0.5d0, 0.2d0]
        values = 0d0
        do i = 0, 4
            values(1, i, :, :) = ex(i)
            values(6, i, :, :) = bz(i)
        end do
        simulator = build_simulator(values, periodic=.true.)
        do direction = -1, 1, 2
            p = new_Particle(1d0, [2d0 + 1.95d0*direction, 2d0, 2d0], [dble(direction), 0.2d0, 0d0])
            forward = advance(simulator, p, -0.2d0)
            call assert_true("test trajectory crossed a periodic face", &
                             abs(forward%position(1) - p%position(1)) > 3d0)
            recovered = advance(simulator, forward, 0.2d0)
            call assert_close_vec("periodic crossing inverse position", recovered%position, p%position)
            call assert_close_vec("periodic crossing inverse velocity", recovered%velocity, p%velocity)
        end do
    end subroutine

    subroutine test_zero_step()
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: p, ret
        double precision :: values(6, 0:4, 0:4, 0:4)

        values = 1d0
        simulator = build_simulator(values)
        p = new_Particle(-2d0, [2d0, 2d0, 2d0], [0.1d0, 0.2d0, 0.3d0])
        p%t = 1.23d0
        ret = advance(simulator, p, 0d0)
        call assert_close_vec("zero step position", ret%position, p%position)
        call assert_close_vec("zero step velocity", ret%velocity, p%velocity)
        call assert_close("zero step time", ret%t, p%t)
    end subroutine

end program

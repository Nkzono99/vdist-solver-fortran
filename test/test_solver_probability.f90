program test_solver_probability
    !! Integration test for `t_Solver%calculate_probability`.
    !!
    !! Builds the smallest possible electrostatic simulator that exercises
    !! the full backtrace -> collision -> probability lookup path, bypassing
    !! the namelist-driven builder: a 2x1x1 zero-field box with a single
    !! Maxwellian-tagged plane at x = 2, a particle travelling toward it,
    !! and a deterministic expected Maxwell PDF value.
    use m_particle, only: t_Particle, new_Particle
    use m_field, only: t_VectorFieldGrid, new_VectorFieldGrid
    use m_probabilities, only: t_Probability, tp_Probability, &
                               t_ZeroProbability, new_ZeroProbability, &
                               t_MaxwellianProbability, new_MaxwellianProbability
    use m_simulator, only: t_ESSimulator, new_ESSimulator
    use m_solver, only: t_Solver, new_Solver, t_BacktraceRecord, t_ProbabilityRecord
    use finbound, only: t_BoundaryList, new_BoundaryList, &
                        t_Boundary, t_PlaneXYZ, new_PlaneX, t_CollisionRecord

    use m_test_helpers, only: assert_close, assert_close_vec, assert_true

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)

    call test_maxwell_plane_returns_pdf_at_hit_velocity()
    call test_backtrace_without_collision_returns_invalid_probability()
    call test_probability_without_collision_returns_final_particle()
    call test_collision_time_uses_signed_step_duration()

    print *, "test_solver_probability: all tests passed."

contains

    subroutine test_collision_time_uses_signed_step_duration()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl, updated
        type(t_CollisionRecord) :: collision
        type(t_ProbabilityRecord) :: probability
        type(t_BacktraceRecord) :: trace
        double precision :: dt, speed, expected, endpoint(3), actual_dt
        integer :: direction, adaptive, step_case

        simulator = build_simulator_with_maxwell_plane_at_x2(sigma=1d0)
        solver = new_Solver(simulator)
        do direction = -1, 1, 2
            do adaptive = 0, 1
                do step_case = 1, 2
                    speed = 2d0
                    dt = dble(direction*step_case)*0.1d0
                    pcl = new_Particle(1d0, [1.95d0, 0.5d0, 0.5d0], &
                                       [-dble(direction)*speed, 0d0, 0d0])
                    pcl%t = 3d0
                    expected = pcl%t + dble(direction)*0.05d0/speed
                    actual_dt = dt
                    if (adaptive == 1) actual_dt = dt/speed
                    updated = simulator%update(pcl, actual_dt, collision, endpoint)
                    call assert_true("simulator reports collision", collision%is_collided)
                    call assert_close("unwrapped endpoint is collision clipped", endpoint(1), 2d0)
                    call assert_close("direct simulator collision time", updated%t, expected)
                    probability = solver%calculate_probability(pcl, dt, 3, adaptive == 1)
                    call assert_true("straight trace hits plane", probability%is_valid)
                    call assert_close("collision time equals signed distance/speed", probability%t, expected)
                    call assert_close("collision particle time agrees", probability%particle%t, expected)
                    trace = solver%backtrace(pcl, dt, 3, 1, adaptive == 1)
                    call assert_close("backtrace collision time agrees", trace%ts(trace%last_step), expected)
                    call assert_close("backtrace collision position", trace%traces(trace%last_step)%position(1), 2d0)
                end do
            end do
        end do
    end subroutine

    subroutine test_maxwell_plane_returns_pdf_at_hit_velocity()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl
        type(t_ProbabilityRecord) :: record
        double precision :: expected_pdf
        double precision :: position(3), velocity(3)
        double precision :: sigma

        simulator = build_simulator_with_maxwell_plane_at_x2(sigma=1d0)
        solver = new_Solver(simulator)

        ! Back-trace takes position forward as position - dt*velocity, so a
        ! negative x-velocity moves the particle toward the +x plane at x=2.
        position = [1d0, 0.5d0, 0.5d0]
        velocity = [-1d0, 0d0, 0d0]
        pcl = new_Particle(1d0, position, velocity)

        record = solver%calculate_probability(pcl, dt=2d0, max_step=5, use_adaptive_dt=.false.)

        call assert_true("probability record is valid after collision", record%is_valid)

        sigma = 1d0
        ! Maxwellian PDF at (vx=-1, vy=0, vz=0) with loc=0 and sigma=1:
        !   f(v) = (1/sqrt(2*pi))^3 * exp(-|v|^2 / 2)
        expected_pdf = (1d0/sqrt(2d0*pi))**3*exp(-0.5d0)

        call assert_close("probability equals Maxwell PDF at collision velocity", &
                          record%probability, expected_pdf, tolerance=1d-10)
    end subroutine

    subroutine test_backtrace_without_collision_returns_invalid_probability()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl
        type(t_BacktraceRecord) :: record
        double precision :: position(3), velocity(3)

        simulator = build_simulator_without_boundaries()
        solver = new_Solver(simulator)

        position = [1d0, 0.5d0, 0.5d0]
        velocity = [0d0, 1d0, 0d0]
        pcl = new_Particle(1d0, position, velocity)

        record = solver%backtrace(pcl, dt=0.25d0, max_step=3, output_interval=1, use_adaptive_dt=.false.)

        call assert_close("backtrace no-collision probability is invalid sentinel", &
                          record%probability, -1d0)
    end subroutine

    subroutine test_probability_without_collision_returns_final_particle()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl
        type(t_ProbabilityRecord) :: record
        double precision :: position(3), velocity(3)

        simulator = build_simulator_without_boundaries()
        solver = new_Solver(simulator)

        position = [1d0, 0.5d0, 0.5d0]
        velocity = [0d0, 1d0, 0d0]
        pcl = new_Particle(1d0, position, velocity)

        record = solver%calculate_probability(pcl, dt=0.25d0, max_step=3, use_adaptive_dt=.false.)

        call assert_true("probability record is invalid without collision", .not. record%is_valid)
        call assert_close("probability no-collision probability is invalid sentinel", &
                          record%probability, -1d0)
        call assert_close_vec("probability no-collision returns final propagated position", &
                              record%particle%position, [1d0, -0.25d0, 0.5d0])
        call assert_close_vec("probability no-collision returns final propagated velocity", &
                              record%particle%velocity, [0d0, 1d0, 0d0])
    end subroutine

    function build_simulator_with_maxwell_plane_at_x2(sigma) result(simulator)
        double precision, intent(in) :: sigma
        type(t_ESSimulator) :: simulator

        type(t_VectorFieldGrid) :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probs(:)

        double precision :: eb_values(6, 0:2, 0:1, 0:1)
        double precision :: locs(3), scales(3)
        integer :: boundary_conditions(3)

        class(t_Boundary), pointer :: pbound
        type(t_PlaneXYZ), pointer :: pplane

        eb_values = 0d0
        eb = new_VectorFieldGrid(6, 2, 1, 1, eb_values)

        allocate (probs(2))
        allocate (probs(1)%ref, source=new_ZeroProbability())
        locs = [0d0, 0d0, 0d0]
        scales = [sigma, sigma, sigma]
        allocate (probs(2)%ref, source=new_MaxwellianProbability(locs, scales))

        boundaries = new_BoundaryList()

        allocate (pplane)
        pplane = new_PlaneX(2d0)
        pbound => pplane
        pbound%material%tag = 2  ! index into probs()
        call boundaries%add_boundary(pbound)

        boundary_conditions = [2, 0, 0]
        simulator = new_ESSimulator(2, 1, 1, boundary_conditions, eb, boundaries, probs)
    end function

    function build_simulator_without_boundaries() result(simulator)
        type(t_ESSimulator) :: simulator

        type(t_VectorFieldGrid) :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probs(:)

        double precision :: eb_values(6, 0:2, 0:1, 0:1)
        integer :: boundary_conditions(3)

        eb_values = 0d0
        eb = new_VectorFieldGrid(6, 2, 1, 1, eb_values)

        allocate (probs(1))
        allocate (probs(1)%ref, source=new_ZeroProbability())

        boundaries = new_BoundaryList()

        boundary_conditions = [2, 2, 2]
        simulator = new_ESSimulator(2, 1, 1, boundary_conditions, eb, boundaries, probs)
    end function

end program

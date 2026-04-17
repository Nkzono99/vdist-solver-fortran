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
    use m_solver, only: t_Solver, new_Solver, t_ProbabilityRecord
    use finbound, only: t_BoundaryList, new_BoundaryList, &
                        t_Boundary, t_PlaneXYZ, new_PlaneX

    use m_test_helpers, only: assert_close, assert_true

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)

    call test_maxwell_plane_returns_pdf_at_hit_velocity()

    print *, "test_solver_probability: all tests passed."

contains

    subroutine test_maxwell_plane_returns_pdf_at_hit_velocity()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl
        type(t_ProbabilityRecord) :: record
        double precision :: expected_pdf
        double precision :: sigma

        simulator = build_simulator_with_maxwell_plane_at_x2(sigma=1d0)
        solver = new_Solver(simulator)

        ! Back-trace takes position forward as position - dt*velocity, so a
        ! negative x-velocity moves the particle toward the +x plane at x=2.
        pcl = new_Particle(1d0, [1d0, 0.5d0, 0.5d0], [-1d0, 0d0, 0d0])

        record = solver%calculate_probability(pcl, dt=2d0, max_step=5, use_adaptive_dt=.false.)

        call assert_true("probability record is valid after collision", record%is_valid)

        sigma = 1d0
        ! Maxwellian PDF at (vx=-1, vy=0, vz=0) with loc=0 and sigma=1:
        !   f(v) = (1/sqrt(2*pi))^3 * exp(-|v|^2 / 2)
        expected_pdf = (1d0/sqrt(2d0*pi))**3*exp(-0.5d0)

        call assert_close("probability equals Maxwell PDF at collision velocity", &
                          record%probability, expected_pdf, tolerance=1d-10)
    end subroutine

    function build_simulator_with_maxwell_plane_at_x2(sigma) result(simulator)
        double precision, intent(in) :: sigma
        type(t_ESSimulator) :: simulator

        type(t_VectorFieldGrid) :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probs(:)

        double precision :: eb_values(6, 0:2, 0:1, 0:1)

        class(t_Boundary), pointer :: pbound
        type(t_PlaneXYZ), pointer :: pplane

        eb_values = 0d0
        eb = new_VectorFieldGrid(6, 2, 1, 1, eb_values)

        allocate (probs(2))
        allocate (probs(1)%ref, source=new_ZeroProbability())
        allocate (probs(2)%ref, source=new_MaxwellianProbability([0d0, 0d0, 0d0], &
                                                                 [sigma, sigma, sigma]))

        boundaries = new_BoundaryList()

        allocate (pplane)
        pplane = new_PlaneX(2d0)
        pbound => pplane
        pbound%material%tag = 2  ! index into probs()
        call boundaries%add_boundary(pbound)

        simulator = new_ESSimulator(2, 1, 1, [2, 0, 0], eb, boundaries, probs)
    end function

end program

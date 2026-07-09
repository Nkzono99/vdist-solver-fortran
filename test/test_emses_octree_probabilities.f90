program test_emses_octree_probabilities
    use m_allcom, only: qm
    use finbound, only: t_BoundaryList, new_BoundaryList, &
                        t_Boundary, t_PlaneXYZ, new_PlaneX
    use m_emses_octree_probabilities, only: STATUS_OK, STATUS_SAMPLE_LIMIT_REACHED, &
                                            STATUS_NO_SIGNAL, &
                                            STATUS_OUTPUT_CAPACITY_EXCEEDED, &
                                            STATUS_INVALID_BOUNDS, &
                                            get_probabilities_octree_impl
    use m_field, only: t_VectorFieldGrid, new_VectorFieldGrid
    use m_particle, only: new_Particle
    use m_probabilities, only: new_MaxwellianProbability, tp_Probability
    use m_simulator, only: t_ESSimulator, new_ESSimulator
    use m_solver, only: t_ProbabilityRecord, t_Solver, new_Solver
    use m_test_helpers, only: assert_equal_int, assert_true

    implicit none

    call test_octree_samples_two_velocity_lobes()
    call test_octree_reports_invalid_bounds()
    call test_octree_reports_no_signal()
    call test_octree_reports_sample_limit()
    call test_octree_reports_leaf_capacity()
    call test_octree_parallel_matches_serial_counts()

    print *, "test_emses_octree_probabilities: all tests passed."

contains

    subroutine test_octree_samples_two_velocity_lobes()
        integer, parameter :: nspatial = 1
        integer, parameter :: max_samples_per_cell = 4096
        integer, parameter :: max_leaves_per_cell = 128
        integer, parameter :: sample_capacity = nspatial*max_samples_per_cell
        integer, parameter :: leaf_capacity = nspatial*max_leaves_per_cell

        type(t_ESSimulator) :: simulator
        real(8) :: spatial_points(3, nspatial)
        real(8) :: velocity_bounds(6, nspatial)
        integer :: sample_spatial_index(sample_capacity)
        real(8) :: velocities(3, sample_capacity)
        real(8) :: probabilities(sample_capacity)
        integer :: leaf_spatial_index(leaf_capacity)
        real(8) :: leaf_bounds(6, leaf_capacity)
        real(8) :: leaf_value_min(leaf_capacity)
        real(8) :: leaf_value_max(leaf_capacity)
        integer :: leaf_depth(leaf_capacity)
        integer :: leaf_sample_start(leaf_capacity)
        integer :: leaf_sample_count(leaf_capacity)
        integer :: status(nspatial)
        integer :: sample_count(nspatial)
        integer :: leaf_count(nspatial)
        integer :: actual_sample_count
        integer :: actual_leaf_count
        integer :: i
        logical :: has_left_lobe
        logical :: has_right_lobe
        logical :: has_invalid_sample
        type(t_Solver) :: solver
        type(t_ProbabilityRecord) :: record
        real(8) :: trace_dt
        real(8) :: velocity(3)
        integer :: scout_bins(3)

        qm(1) = 1d0
        simulator = build_two_plane_simulator()
        spatial_points(:, 1) = [1d0, 0.5d0, 0.5d0]
        velocity_bounds(:, 1) = [-2d0, 2d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0]
        solver = new_Solver(simulator)
        trace_dt = 0.3d0
        scout_bins = [5, 3, 3]

        velocity = [-1d0, 0d0, 0d0]
        record = solver%calculate_probability( &
            new_Particle(1d0, spatial_points(:, 1), velocity), &
            dt=trace_dt, max_step=16, use_adaptive_dt=.false.)
        call assert_true("direct solver samples negative-vx boundary", &
                         record%is_valid .and. record%probability > 0d0)

        velocity = [1d0, 0d0, 0d0]
        record = solver%calculate_probability( &
            new_Particle(1d0, spatial_points(:, 1), velocity), &
            dt=trace_dt, max_step=16, use_adaptive_dt=.false.)
        call assert_true("direct solver samples positive-vx boundary", &
                         record%is_valid .and. record%probability > 0d0)

        call get_probabilities_octree_impl( &
            simulator=simulator, &
            ispec=1, &
            nspatial=nspatial, &
            spatial_points=spatial_points, &
            velocity_bounds=velocity_bounds, &
            dt=trace_dt, &
            max_step=16, &
            use_adaptive_dt=0, &
            scout_bins=scout_bins, &
            max_depth=2, &
            max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, &
            refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, &
            expand_factor=1.25d0, &
            max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, &
            return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, &
            return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, &
            return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, &
            return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, &
            return_status=status, &
            return_sample_count=sample_count, &
            return_leaf_count=leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, &
            n_threads=1)

        call assert_equal_int("octree spatial point status is OK", status(1), STATUS_OK)
        call assert_true("octree returns samples", sample_count(1) > 0)
        call assert_true("octree returns leaves", leaf_count(1) > 0)

        has_left_lobe = .false.
        has_right_lobe = .false.
        has_invalid_sample = .false.
        do i = 1, sample_count(1)
            if (probabilities(i) == -1d0) has_invalid_sample = .true.
            if (probabilities(i) <= 0d0) cycle
            if (velocities(1, i) < -0.5d0) has_left_lobe = .true.
            if (velocities(1, i) > 0.5d0) has_right_lobe = .true.
        end do

        call assert_true("octree samples negative-vx lobe", has_left_lobe)
        call assert_true("octree samples positive-vx lobe", has_right_lobe)
        call assert_true("octree preserves invalid probability sentinel", has_invalid_sample)
    end subroutine

    subroutine test_octree_reports_invalid_bounds()
        integer, parameter :: nspatial = 1
        integer, parameter :: max_samples_per_cell = 8
        integer, parameter :: max_leaves_per_cell = 4
        integer, parameter :: sample_capacity = nspatial*max_samples_per_cell
        integer, parameter :: leaf_capacity = nspatial*max_leaves_per_cell

        type(t_ESSimulator) :: simulator
        real(8) :: spatial_points(3, nspatial)
        real(8) :: velocity_bounds(6, nspatial)
        integer :: sample_spatial_index(sample_capacity)
        real(8) :: velocities(3, sample_capacity)
        real(8) :: probabilities(sample_capacity)
        integer :: leaf_spatial_index(leaf_capacity)
        real(8) :: leaf_bounds(6, leaf_capacity)
        real(8) :: leaf_value_min(leaf_capacity)
        real(8) :: leaf_value_max(leaf_capacity)
        integer :: leaf_depth(leaf_capacity)
        integer :: leaf_sample_start(leaf_capacity)
        integer :: leaf_sample_count(leaf_capacity)
        integer :: status(nspatial)
        integer :: sample_count(nspatial)
        integer :: leaf_count(nspatial)
        integer :: actual_sample_count
        integer :: actual_leaf_count
        integer :: scout_bins(3)

        qm(1) = 1d0
        simulator = build_two_plane_simulator()
        spatial_points(:, 1) = [1d0, 0.5d0, 0.5d0]
        velocity_bounds(:, 1) = [2d0, -2d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0]
        scout_bins = [3, 3, 3]

        call get_probabilities_octree_impl( &
            simulator=simulator, ispec=1, nspatial=nspatial, &
            spatial_points=spatial_points, velocity_bounds=velocity_bounds, &
            dt=0.3d0, max_step=16, use_adaptive_dt=0, scout_bins=scout_bins, &
            max_depth=1, max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, expand_factor=1.25d0, max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, return_status=status, &
            return_sample_count=sample_count, return_leaf_count=leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, n_threads=1)

        call assert_equal_int("octree reports invalid bounds", status(1), STATUS_INVALID_BOUNDS)
        call assert_equal_int("invalid bounds produce no samples", sample_count(1), 0)
        call assert_equal_int("invalid bounds produce no leaves", leaf_count(1), 0)
    end subroutine

    subroutine test_octree_reports_no_signal()
        integer, parameter :: nspatial = 1
        integer, parameter :: max_samples_per_cell = 64
        integer, parameter :: max_leaves_per_cell = 8
        integer, parameter :: sample_capacity = nspatial*max_samples_per_cell
        integer, parameter :: leaf_capacity = nspatial*max_leaves_per_cell

        type(t_ESSimulator) :: simulator
        real(8) :: spatial_points(3, nspatial)
        real(8) :: velocity_bounds(6, nspatial)
        integer :: sample_spatial_index(sample_capacity)
        real(8) :: velocities(3, sample_capacity)
        real(8) :: probabilities(sample_capacity)
        integer :: leaf_spatial_index(leaf_capacity)
        real(8) :: leaf_bounds(6, leaf_capacity)
        real(8) :: leaf_value_min(leaf_capacity)
        real(8) :: leaf_value_max(leaf_capacity)
        integer :: leaf_depth(leaf_capacity)
        integer :: leaf_sample_start(leaf_capacity)
        integer :: leaf_sample_count(leaf_capacity)
        integer :: status(nspatial)
        integer :: sample_count(nspatial)
        integer :: leaf_count(nspatial)
        integer :: actual_sample_count
        integer :: actual_leaf_count
        integer :: scout_bins(3)

        qm(1) = 1d0
        simulator = build_two_plane_simulator()
        spatial_points(:, 1) = [1d0, 0.5d0, 0.5d0]
        velocity_bounds(:, 1) = [-0.1d0, 0.1d0, 1d0, 2d0, -0.1d0, 0.1d0]
        scout_bins = [3, 3, 3]

        call get_probabilities_octree_impl( &
            simulator=simulator, &
            ispec=1, &
            nspatial=nspatial, &
            spatial_points=spatial_points, &
            velocity_bounds=velocity_bounds, &
            dt=0.3d0, &
            max_step=16, &
            use_adaptive_dt=0, &
            scout_bins=scout_bins, &
            max_depth=0, &
            max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, &
            refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, &
            expand_factor=1.25d0, &
            max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, &
            return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, &
            return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, &
            return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, &
            return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, &
            return_status=status, &
            return_sample_count=sample_count, &
            return_leaf_count=leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, &
            n_threads=1)

        call assert_equal_int("octree reports no signal", status(1), STATUS_NO_SIGNAL)
        call assert_true("octree still returns scout samples for no-signal boxes", sample_count(1) > 0)
        call assert_true("octree no-signal samples keep invalid sentinel", &
                         all(probabilities(1:sample_count(1)) == -1d0))
    end subroutine

    subroutine test_octree_reports_sample_limit()
        integer, parameter :: nspatial = 1
        integer, parameter :: max_samples_per_cell = 16
        integer, parameter :: max_leaves_per_cell = 128
        integer, parameter :: sample_capacity = nspatial*max_samples_per_cell
        integer, parameter :: leaf_capacity = nspatial*max_leaves_per_cell

        type(t_ESSimulator) :: simulator
        real(8) :: spatial_points(3, nspatial)
        real(8) :: velocity_bounds(6, nspatial)
        integer :: sample_spatial_index(sample_capacity)
        real(8) :: velocities(3, sample_capacity)
        real(8) :: probabilities(sample_capacity)
        integer :: leaf_spatial_index(leaf_capacity)
        real(8) :: leaf_bounds(6, leaf_capacity)
        real(8) :: leaf_value_min(leaf_capacity)
        real(8) :: leaf_value_max(leaf_capacity)
        integer :: leaf_depth(leaf_capacity)
        integer :: leaf_sample_start(leaf_capacity)
        integer :: leaf_sample_count(leaf_capacity)
        integer :: status(nspatial)
        integer :: sample_count(nspatial)
        integer :: leaf_count(nspatial)
        integer :: actual_sample_count
        integer :: actual_leaf_count
        integer :: scout_bins(3)

        qm(1) = 1d0
        simulator = build_two_plane_simulator()
        spatial_points(:, 1) = [1d0, 0.5d0, 0.5d0]
        velocity_bounds(:, 1) = [-2d0, 2d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0]
        scout_bins = [5, 3, 3]

        call get_probabilities_octree_impl( &
            simulator=simulator, &
            ispec=1, &
            nspatial=nspatial, &
            spatial_points=spatial_points, &
            velocity_bounds=velocity_bounds, &
            dt=0.3d0, &
            max_step=16, &
            use_adaptive_dt=0, &
            scout_bins=scout_bins, &
            max_depth=2, &
            max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, &
            refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, &
            expand_factor=1.25d0, &
            max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, &
            return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, &
            return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, &
            return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, &
            return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, &
            return_status=status, &
            return_sample_count=sample_count, &
            return_leaf_count=leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, &
            n_threads=1)

        call assert_equal_int("octree reports sample limit", status(1), STATUS_SAMPLE_LIMIT_REACHED)
        call assert_equal_int("octree caps samples at configured stride", &
                              sample_count(1), max_samples_per_cell)
    end subroutine

    subroutine test_octree_reports_leaf_capacity()
        integer, parameter :: nspatial = 1
        integer, parameter :: max_samples_per_cell = 128
        integer, parameter :: max_leaves_per_cell = 1
        integer, parameter :: sample_capacity = nspatial*max_samples_per_cell
        integer, parameter :: leaf_capacity = nspatial*max_leaves_per_cell

        type(t_ESSimulator) :: simulator
        real(8) :: spatial_points(3, nspatial)
        real(8) :: velocity_bounds(6, nspatial)
        integer :: sample_spatial_index(sample_capacity)
        real(8) :: velocities(3, sample_capacity)
        real(8) :: probabilities(sample_capacity)
        integer :: leaf_spatial_index(leaf_capacity)
        real(8) :: leaf_bounds(6, leaf_capacity)
        real(8) :: leaf_value_min(leaf_capacity)
        real(8) :: leaf_value_max(leaf_capacity)
        integer :: leaf_depth(leaf_capacity)
        integer :: leaf_sample_start(leaf_capacity)
        integer :: leaf_sample_count(leaf_capacity)
        integer :: status(nspatial)
        integer :: sample_count(nspatial)
        integer :: leaf_count(nspatial)
        integer :: actual_sample_count
        integer :: actual_leaf_count
        integer :: scout_bins(3)

        qm(1) = 1d0
        simulator = build_two_plane_simulator()
        spatial_points(:, 1) = [1d0, 0.5d0, 0.5d0]
        velocity_bounds(:, 1) = [-2d0, 2d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0]
        scout_bins = [5, 3, 3]

        call get_probabilities_octree_impl( &
            simulator=simulator, ispec=1, nspatial=nspatial, &
            spatial_points=spatial_points, velocity_bounds=velocity_bounds, &
            dt=0.3d0, max_step=16, use_adaptive_dt=0, scout_bins=scout_bins, &
            max_depth=2, max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, expand_factor=1.25d0, max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, return_status=status, &
            return_sample_count=sample_count, return_leaf_count=leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, n_threads=1)

        call assert_equal_int("octree reports leaf capacity", status(1), STATUS_OUTPUT_CAPACITY_EXCEEDED)
        call assert_equal_int("octree caps leaves at configured stride", &
                              leaf_count(1), max_leaves_per_cell)
    end subroutine

    subroutine test_octree_parallel_matches_serial_counts()
        integer, parameter :: nspatial = 2
        integer, parameter :: max_samples_per_cell = 1024
        integer, parameter :: max_leaves_per_cell = 128
        integer, parameter :: sample_capacity = nspatial*max_samples_per_cell
        integer, parameter :: leaf_capacity = nspatial*max_leaves_per_cell

        type(t_ESSimulator) :: simulator
        real(8) :: spatial_points(3, nspatial)
        real(8) :: velocity_bounds(6, nspatial)
        integer :: sample_spatial_index(sample_capacity)
        real(8) :: velocities(3, sample_capacity)
        real(8) :: probabilities(sample_capacity)
        integer :: leaf_spatial_index(leaf_capacity)
        real(8) :: leaf_bounds(6, leaf_capacity)
        real(8) :: leaf_value_min(leaf_capacity)
        real(8) :: leaf_value_max(leaf_capacity)
        integer :: leaf_depth(leaf_capacity)
        integer :: leaf_sample_start(leaf_capacity)
        integer :: leaf_sample_count(leaf_capacity)
        integer :: serial_status(nspatial)
        integer :: serial_sample_count(nspatial)
        integer :: serial_leaf_count(nspatial)
        integer :: parallel_status(nspatial)
        integer :: parallel_sample_count(nspatial)
        integer :: parallel_leaf_count(nspatial)
        integer :: actual_sample_count
        integer :: actual_leaf_count
        integer :: scout_bins(3)

        qm(1) = 1d0
        simulator = build_two_plane_simulator()
        spatial_points(:, 1) = [1d0, 0.5d0, 0.5d0]
        spatial_points(:, 2) = [1d0, 0.25d0, 0.5d0]
        velocity_bounds(:, 1) = [-2d0, 2d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0]
        velocity_bounds(:, 2) = velocity_bounds(:, 1)
        scout_bins = [5, 3, 3]

        call get_probabilities_octree_impl( &
            simulator=simulator, ispec=1, nspatial=nspatial, &
            spatial_points=spatial_points, velocity_bounds=velocity_bounds, &
            dt=0.3d0, max_step=16, use_adaptive_dt=0, scout_bins=scout_bins, &
            max_depth=1, max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, expand_factor=1.25d0, max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, return_status=serial_status, &
            return_sample_count=serial_sample_count, return_leaf_count=serial_leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, n_threads=1)

        call get_probabilities_octree_impl( &
            simulator=simulator, ispec=1, nspatial=nspatial, &
            spatial_points=spatial_points, velocity_bounds=velocity_bounds, &
            dt=0.3d0, max_step=16, use_adaptive_dt=0, scout_bins=scout_bins, &
            max_depth=1, max_samples_per_cell=max_samples_per_cell, &
            max_leaves_per_cell=max_leaves_per_cell, refine_threshold_rel=1d-4, &
            edge_threshold_rel=1d-4, expand_factor=1.25d0, max_expansions=0, &
            return_sample_spatial_index=sample_spatial_index, &
            return_velocities=velocities, return_probabilities=probabilities, &
            return_leaf_spatial_index=leaf_spatial_index, return_leaf_bounds=leaf_bounds, &
            return_leaf_value_min=leaf_value_min, return_leaf_value_max=leaf_value_max, &
            return_leaf_depth=leaf_depth, return_leaf_sample_start=leaf_sample_start, &
            return_leaf_sample_count=leaf_sample_count, return_status=parallel_status, &
            return_sample_count=parallel_sample_count, return_leaf_count=parallel_leaf_count, &
            return_actual_sample_count=actual_sample_count, &
            return_actual_leaf_count=actual_leaf_count, n_threads=2)

        call assert_true("parallel octree status matches serial", all(parallel_status == serial_status))
        call assert_true("parallel octree sample counts match serial", &
                         all(parallel_sample_count == serial_sample_count))
        call assert_true("parallel octree leaf counts match serial", &
                         all(parallel_leaf_count == serial_leaf_count))
    end subroutine

    function build_two_plane_simulator() result(simulator)
        type(t_ESSimulator) :: simulator

        type(t_VectorFieldGrid) :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probs(:)
        real(8) :: eb_values(6, 0:2, 0:1, 0:1)
        integer :: boundary_conditions(3)
        real(8) :: locs(3)
        real(8) :: scales(3)

        class(t_Boundary), pointer :: pbound
        type(t_PlaneXYZ), pointer :: left_plane
        type(t_PlaneXYZ), pointer :: right_plane

        eb_values = 0d0
        eb = new_VectorFieldGrid(6, 2, 1, 1, eb_values)

        allocate (probs(3))
        locs = [0d0, 0d0, 0d0]
        scales = [1d0, 1d0, 1d0]
        allocate (probs(1)%ref, source=new_MaxwellianProbability(locs, scales))
        locs = [-1d0, 0d0, 0d0]
        scales = [0.25d0, 0.25d0, 0.25d0]
        allocate (probs(2)%ref, source=new_MaxwellianProbability(locs, scales))
        locs = [1d0, 0d0, 0d0]
        allocate (probs(3)%ref, source=new_MaxwellianProbability(locs, scales))

        boundaries = new_BoundaryList()

        allocate (left_plane)
        left_plane = new_PlaneX(0d0)
        pbound => left_plane
        pbound%material%tag = 3
        call boundaries%add_boundary(pbound)

        allocate (right_plane)
        right_plane = new_PlaneX(2d0)
        pbound => right_plane
        pbound%material%tag = 2
        call boundaries%add_boundary(pbound)

        boundary_conditions = [2, 0, 0]
        simulator = new_ESSimulator(2, 1, 1, boundary_conditions, eb, boundaries, probs)
    end function

end program

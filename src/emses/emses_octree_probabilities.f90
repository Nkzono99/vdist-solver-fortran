module m_emses_octree_probabilities
    !! Adaptive velocity-octree probability evaluation for EMSES solvers.

    use, intrinsic :: iso_c_binding, only: c_double, c_int

!$  use omp_lib, only: omp_set_num_threads

    use m_allcom, only: qm
    use m_particle, only: t_Particle, new_Particle
    use m_simulator, only: t_ESSimulator
    use m_solver, only: t_ProbabilityRecord, t_Solver, new_Solver

    implicit none

    private
    public get_probabilities_octree_impl
    public STATUS_OK
    public STATUS_NO_SIGNAL
    public STATUS_SAMPLE_LIMIT_REACHED
    public STATUS_ROOT_EXPANDED
    public STATUS_EXPANSION_LIMIT_REACHED
    public STATUS_OUTPUT_CAPACITY_EXCEEDED
    public STATUS_INVALID_BOUNDS

    integer(c_int), parameter :: STATUS_OK = 0
    integer(c_int), parameter :: STATUS_NO_SIGNAL = 1
    integer(c_int), parameter :: STATUS_SAMPLE_LIMIT_REACHED = 2
    integer(c_int), parameter :: STATUS_ROOT_EXPANDED = 3
    integer(c_int), parameter :: STATUS_EXPANSION_LIMIT_REACHED = STATUS_ROOT_EXPANDED
    integer(c_int), parameter :: STATUS_OUTPUT_CAPACITY_EXCEEDED = 4
    integer(c_int), parameter :: STATUS_INVALID_BOUNDS = 5

contains

    subroutine get_probabilities_octree_impl( &
        simulator, &
        ispec, &
        nspatial, &
        spatial_points, &
        velocity_bounds, &
        dt, &
        max_step, &
        use_adaptive_dt, &
        scout_bins, &
        max_depth, &
        max_samples_per_cell, &
        max_leaves_per_cell, &
        refine_threshold_rel, &
        edge_threshold_rel, &
        expand_factor, &
        max_expansions, &
        return_sample_spatial_index, &
        return_velocities, &
        return_probabilities, &
        return_leaf_spatial_index, &
        return_leaf_bounds, &
        return_leaf_value_min, &
        return_leaf_value_max, &
        return_leaf_depth, &
        return_leaf_sample_start, &
        return_leaf_sample_count, &
        return_status, &
        return_sample_count, &
        return_leaf_count, &
        return_actual_sample_count, &
        return_actual_leaf_count, &
        n_threads)
        type(t_ESSimulator), intent(in) :: simulator
        integer(c_int), value, intent(in) :: ispec
        integer(c_int), value, intent(in) :: nspatial
        real(c_double), intent(in) :: spatial_points(3, nspatial)
        real(c_double), intent(in) :: velocity_bounds(6, nspatial)
        real(c_double), value, intent(in) :: dt
        integer(c_int), value, intent(in) :: max_step
        integer(c_int), value, intent(in) :: use_adaptive_dt
        integer(c_int), intent(in) :: scout_bins(3)
        integer(c_int), value, intent(in) :: max_depth
        integer(c_int), value, intent(in) :: max_samples_per_cell
        integer(c_int), value, intent(in) :: max_leaves_per_cell
        real(c_double), value, intent(in) :: refine_threshold_rel
        real(c_double), value, intent(in) :: edge_threshold_rel
        real(c_double), value, intent(in) :: expand_factor
        integer(c_int), value, intent(in) :: max_expansions
        integer(c_int), intent(out) :: return_sample_spatial_index(nspatial*max_samples_per_cell)
        real(c_double), intent(out) :: return_velocities(3, nspatial*max_samples_per_cell)
        real(c_double), intent(out) :: return_probabilities(nspatial*max_samples_per_cell)
        integer(c_int), intent(out) :: return_leaf_spatial_index(nspatial*max_leaves_per_cell)
        real(c_double), intent(out) :: return_leaf_bounds(6, nspatial*max_leaves_per_cell)
        real(c_double), intent(out) :: return_leaf_value_min(nspatial*max_leaves_per_cell)
        real(c_double), intent(out) :: return_leaf_value_max(nspatial*max_leaves_per_cell)
        integer(c_int), intent(out) :: return_leaf_depth(nspatial*max_leaves_per_cell)
        integer(c_int), intent(out) :: return_leaf_sample_start(nspatial*max_leaves_per_cell)
        integer(c_int), intent(out) :: return_leaf_sample_count(nspatial*max_leaves_per_cell)
        integer(c_int), intent(out) :: return_status(nspatial)
        integer(c_int), intent(out) :: return_sample_count(nspatial)
        integer(c_int), intent(out) :: return_leaf_count(nspatial)
        integer(c_int), intent(out) :: return_actual_sample_count
        integer(c_int), intent(out) :: return_actual_leaf_count
        integer(c_int), value, optional, intent(in) :: n_threads

        integer :: ispatial

        return_sample_spatial_index = -1
        return_velocities = 0d0
        return_probabilities = -1d0
        return_leaf_spatial_index = -1
        return_leaf_bounds = 0d0
        return_leaf_value_min = 0d0
        return_leaf_value_max = 0d0
        return_leaf_depth = 0
        return_leaf_sample_start = 0
        return_leaf_sample_count = 0
        return_status = STATUS_NO_SIGNAL
        return_sample_count = 0
        return_leaf_count = 0

!$      if (present(n_threads)) then
!$          call omp_set_num_threads(max(1, int(n_threads)))
!$      end if

        !$omp parallel do schedule(dynamic, 1)
        do ispatial = 1, nspatial
            call evaluate_spatial_octree(ispatial)
        end do
        !$omp end parallel do

        return_actual_sample_count = sum(return_sample_count)
        return_actual_leaf_count = sum(return_leaf_count)

    contains

        subroutine evaluate_spatial_octree(ispatial)
            integer, intent(in) :: ispatial

            type(t_Solver) :: solver
            real(c_double) :: root_bounds(6)
            integer :: status
            logical :: has_signal

            solver = new_Solver(simulator)
            root_bounds = velocity_bounds(:, ispatial)
            status = STATUS_OK

            if (.not. valid_bounds(root_bounds)) then
                return_status(ispatial) = STATUS_INVALID_BOUNDS
                return
            end if

            call expand_root_if_needed(solver, ispatial, root_bounds, status)
            call traverse_octree(solver, ispatial, root_bounds, status)
            has_signal = has_positive_signal(ispatial)

            if (return_sample_count(ispatial) <= 0) then
                return_status(ispatial) = STATUS_NO_SIGNAL
            else if (status /= STATUS_OK) then
                return_status(ispatial) = status
            else if (.not. has_signal) then
                return_status(ispatial) = STATUS_NO_SIGNAL
            else
                return_status(ispatial) = STATUS_OK
            end if
        end subroutine

        logical function valid_bounds(bounds) result(ret)
            real(c_double), intent(in) :: bounds(6)

            ret = bounds(1) < bounds(2) .and. bounds(3) < bounds(4) .and. bounds(5) < bounds(6)
        end function

        logical function has_positive_signal(ispatial) result(ret)
            integer, intent(in) :: ispatial

            integer :: leaf_start
            integer :: leaf_stop

            if (return_leaf_count(ispatial) <= 0) then
                ret = .false.
                return
            end if

            leaf_start = (ispatial - 1)*max_leaves_per_cell + 1
            leaf_stop = leaf_start + return_leaf_count(ispatial) - 1
            ret = any(return_leaf_value_max(leaf_start:leaf_stop) > 0d0)
        end function

        subroutine expand_root_if_needed(solver, ispatial, bounds, status)
            type(t_Solver), intent(in) :: solver
            integer, intent(in) :: ispatial
            real(c_double), intent(inout) :: bounds(6)
            integer, intent(inout) :: status

            integer :: iexpand
            real(c_double) :: pmax, edge_pmax

            if (max_expansions <= 0 .or. expand_factor <= 1d0) then
                return
            end if

            do iexpand = 1, max_expansions
                call score_box(solver, spatial_points(:, ispatial), bounds, scout_bins, pmax, edge_pmax)
                if (pmax <= 0d0) return
                if (edge_pmax/max(pmax, tiny(1d0)) <= edge_threshold_rel) return
                call expand_bounds(bounds, expand_factor)
                status = STATUS_ROOT_EXPANDED
            end do
        end subroutine

        subroutine traverse_octree(solver, ispatial, root_bounds, status)
            type(t_Solver), intent(in) :: solver
            integer, intent(in) :: ispatial
            real(c_double), intent(in) :: root_bounds(6)
            integer, intent(inout) :: status

            real(c_double), allocatable :: queue_bounds(:, :)
            integer(c_int), allocatable :: queue_depth(:)
            integer :: head, tail
            real(c_double) :: bounds(6)
            integer :: depth
            real(c_double) :: pmin, pmax
            integer :: sample_start
            integer :: sample_added
            logical :: should_split

            allocate (queue_bounds(6, max_leaves_per_cell))
            allocate (queue_depth(max_leaves_per_cell))

            head = 1
            tail = 1
            queue_bounds(:, 1) = root_bounds
            queue_depth(1) = 0

            do while (head <= tail)
                bounds = queue_bounds(:, head)
                depth = queue_depth(head)
                head = head + 1

                if (return_leaf_count(ispatial) >= max_leaves_per_cell) then
                    status = STATUS_OUTPUT_CAPACITY_EXCEEDED
                    exit
                end if

                sample_start = return_sample_count(ispatial)
                call sample_leaf(solver, ispatial, bounds, depth, pmin, pmax, sample_added, status)
                call append_leaf(ispatial, bounds, depth, pmin, pmax, sample_start, sample_added)

                if (status == STATUS_SAMPLE_LIMIT_REACHED) exit
                should_split = depth < max_depth .and. pmax > 0d0 &
                               .and. (pmax - pmin) >= refine_threshold_rel*max(pmax, tiny(1d0))
                if (should_split) then
                    call append_children(bounds, depth + 1, queue_bounds, queue_depth, tail, status)
                    if (status == STATUS_OUTPUT_CAPACITY_EXCEEDED) exit
                end if
            end do
        end subroutine

        subroutine sample_leaf(solver, ispatial, bounds, depth, pmin, pmax, sample_added, status)
            type(t_Solver), intent(in) :: solver
            integer, intent(in) :: ispatial
            real(c_double), intent(in) :: bounds(6)
            integer, intent(in) :: depth
            real(c_double), intent(out) :: pmin
            real(c_double), intent(out) :: pmax
            integer, intent(out) :: sample_added
            integer, intent(inout) :: status

            integer :: bins(3)

            if (depth == 0) then
                bins = max(2, scout_bins)
            else
                bins = [3, 3, 3]
            end if
            call sample_grid(solver, ispatial, bounds, bins, pmin, pmax, sample_added, status, append_samples=.true.)
        end subroutine

        subroutine score_box(solver, position, bounds, bins, pmax, edge_pmax)
            type(t_Solver), intent(in) :: solver
            real(c_double), intent(in) :: position(3)
            real(c_double), intent(in) :: bounds(6)
            integer, intent(in) :: bins(3)
            real(c_double), intent(out) :: pmax
            real(c_double), intent(out) :: edge_pmax

            integer :: ivx, ivy, ivz
            real(c_double) :: velocity(3)
            real(c_double) :: prob
            real(c_double) :: score
            integer :: index(3)
            integer :: safe_bins(3)

            pmax = 0d0
            edge_pmax = 0d0
            safe_bins = max(2, bins)
            do ivz = 1, safe_bins(3)
                do ivy = 1, safe_bins(2)
                    do ivx = 1, safe_bins(1)
                        index = [ivx, ivy, ivz]
                        velocity = velocity_at_grid(bounds, index, safe_bins)
                        prob = evaluate_probability(solver, position, velocity)
                        score = max(0d0, prob)
                        pmax = max(pmax, score)
                        if (ivx == 1 .or. ivx == safe_bins(1) &
                            .or. ivy == 1 .or. ivy == safe_bins(2) &
                            .or. ivz == 1 .or. ivz == safe_bins(3)) then
                            edge_pmax = max(edge_pmax, score)
                        end if
                    end do
                end do
            end do
        end subroutine

        subroutine sample_grid(solver, ispatial, bounds, bins, pmin, pmax, sample_added, status, append_samples)
            type(t_Solver), intent(in) :: solver
            integer, intent(in) :: ispatial
            real(c_double), intent(in) :: bounds(6)
            integer, intent(in) :: bins(3)
            real(c_double), intent(out) :: pmin
            real(c_double), intent(out) :: pmax
            integer, intent(out) :: sample_added
            integer, intent(inout) :: status
            logical, intent(in) :: append_samples

            integer :: ivx, ivy, ivz
            real(c_double) :: velocity(3)
            real(c_double) :: prob
            real(c_double) :: score
            integer :: index(3)

            pmin = huge(1d0)
            pmax = 0d0
            sample_added = 0
            do ivz = 1, bins(3)
                do ivy = 1, bins(2)
                    do ivx = 1, bins(1)
                        index = [ivx, ivy, ivz]
                        velocity = velocity_at_grid(bounds, index, bins)
                        prob = evaluate_probability(solver, spatial_points(:, ispatial), velocity)
                        if (append_samples) then
                            call append_sample(ispatial, velocity, prob, status)
                            if (status == STATUS_SAMPLE_LIMIT_REACHED) then
                                if (pmin == huge(1d0)) pmin = 0d0
                                return
                            end if
                            sample_added = sample_added + 1
                        end if
                        score = max(0d0, prob)
                        pmin = min(pmin, score)
                        pmax = max(pmax, score)
                    end do
                end do
            end do
            if (pmin == huge(1d0)) pmin = 0d0
        end subroutine

        real(c_double) function evaluate_probability(solver, position, velocity) result(ret)
            type(t_Solver), intent(in) :: solver
            real(c_double), intent(in) :: position(3)
            real(c_double), intent(in) :: velocity(3)

            type(t_Particle) :: particle
            type(t_ProbabilityRecord) :: record

            particle = new_Particle(qm(ispec), position, velocity)
            record = solver%calculate_probability(particle, dt, max_step, use_adaptive_dt == 1)
            if (record%is_valid) then
                ret = max(0d0, record%probability)
            else
                ret = -1d0
            end if
        end function

        subroutine append_sample(ispatial, velocity, probability, status)
            integer, intent(in) :: ispatial
            real(c_double), intent(in) :: velocity(3)
            real(c_double), intent(in) :: probability
            integer, intent(inout) :: status

            integer :: local_index
            integer :: global_index

            if (return_sample_count(ispatial) >= max_samples_per_cell) then
                status = STATUS_SAMPLE_LIMIT_REACHED
                return
            end if

            local_index = return_sample_count(ispatial) + 1
            global_index = (ispatial - 1)*max_samples_per_cell + local_index

            return_sample_spatial_index(global_index) = ispatial - 1
            return_velocities(:, global_index) = velocity
            return_probabilities(global_index) = probability
            return_sample_count(ispatial) = local_index
        end subroutine

        subroutine append_leaf(ispatial, bounds, depth, pmin, pmax, sample_start, sample_added)
            integer, intent(in) :: ispatial
            real(c_double), intent(in) :: bounds(6)
            integer, intent(in) :: depth
            real(c_double), intent(in) :: pmin
            real(c_double), intent(in) :: pmax
            integer, intent(in) :: sample_start
            integer, intent(in) :: sample_added

            integer :: local_index
            integer :: global_index

            local_index = return_leaf_count(ispatial) + 1
            global_index = (ispatial - 1)*max_leaves_per_cell + local_index

            return_leaf_spatial_index(global_index) = ispatial - 1
            return_leaf_bounds(:, global_index) = bounds
            return_leaf_value_min(global_index) = pmin
            return_leaf_value_max(global_index) = pmax
            return_leaf_depth(global_index) = depth
            return_leaf_sample_start(global_index) = sample_start
            return_leaf_sample_count(global_index) = sample_added
            return_leaf_count(ispatial) = local_index
        end subroutine

        subroutine append_children(bounds, child_depth, queue_bounds, queue_depth, tail, status)
            real(c_double), intent(in) :: bounds(6)
            integer, intent(in) :: child_depth
            real(c_double), intent(inout) :: queue_bounds(:, :)
            integer(c_int), intent(inout) :: queue_depth(:)
            integer, intent(inout) :: tail
            integer, intent(inout) :: status

            real(c_double) :: mid(3)
            integer :: sx, sy, sz
            real(c_double) :: child(6)

            mid = [0.5d0*(bounds(1) + bounds(2)), &
                   0.5d0*(bounds(3) + bounds(4)), &
                   0.5d0*(bounds(5) + bounds(6))]

            do sz = 0, 1
                do sy = 0, 1
                    do sx = 0, 1
                        if (tail >= size(queue_bounds, 2)) then
                            status = STATUS_OUTPUT_CAPACITY_EXCEEDED
                            return
                        end if
                        child = bounds
                        if (sx == 0) then
                            child(2) = mid(1)
                        else
                            child(1) = mid(1)
                        end if
                        if (sy == 0) then
                            child(4) = mid(2)
                        else
                            child(3) = mid(2)
                        end if
                        if (sz == 0) then
                            child(6) = mid(3)
                        else
                            child(5) = mid(3)
                        end if
                        tail = tail + 1
                        queue_bounds(:, tail) = child
                        queue_depth(tail) = child_depth
                    end do
                end do
            end do
        end subroutine

        function velocity_at_grid(bounds, index, bins) result(velocity)
            real(c_double), intent(in) :: bounds(6)
            integer, intent(in) :: index(3)
            integer, intent(in) :: bins(3)
            real(c_double) :: velocity(3)

            velocity(1) = linspace_value(bounds(1), bounds(2), index(1), bins(1))
            velocity(2) = linspace_value(bounds(3), bounds(4), index(2), bins(2))
            velocity(3) = linspace_value(bounds(5), bounds(6), index(3), bins(3))
        end function

        real(c_double) function linspace_value(vmin, vmax, index, count) result(ret)
            real(c_double), intent(in) :: vmin
            real(c_double), intent(in) :: vmax
            integer, intent(in) :: index
            integer, intent(in) :: count

            if (count <= 1) then
                ret = 0.5d0*(vmin + vmax)
            else
                ret = vmin + (vmax - vmin)*dble(index - 1)/dble(count - 1)
            end if
        end function

        subroutine expand_bounds(bounds, factor)
            real(c_double), intent(inout) :: bounds(6)
            real(c_double), intent(in) :: factor

            real(c_double) :: center
            real(c_double) :: half_width
            integer :: axis

            do axis = 1, 3
                center = 0.5d0*(bounds(2*axis - 1) + bounds(2*axis))
                half_width = 0.5d0*(bounds(2*axis) - bounds(2*axis - 1))*factor
                bounds(2*axis - 1) = center - half_width
                bounds(2*axis) = center + half_width
            end do
        end subroutine

    end subroutine

end module

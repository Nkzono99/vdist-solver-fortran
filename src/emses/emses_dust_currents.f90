module m_emses_dust_currents
    !! Unit-density OML flux moments in EMSES normalised units.
    use, intrinsic :: iso_c_binding, only: c_double, c_int
!$  use omp_lib, only: omp_set_num_threads
    use m_allcom, only: qm
    use m_simulator, only: t_ESSimulator
    use m_emses_octree_probabilities, only: get_probabilities_octree_impl
    implicit none
    private
    public evaluate_dust_currents_impl
contains
    subroutine evaluate_dust_currents_impl(simulator, ispec, nspatial, spatial_points, velocity_bounds, &
        npsi, psi_grid, dt, max_step, use_adaptive_dt, scout_bins, max_depth, max_samples_per_cell, &
        max_leaves_per_cell, refine_threshold_rel, edge_threshold_rel, expand_factor, max_expansions, &
        flux, flux_error, status, sample_count, terminal_box_count, n_threads)
        type(t_ESSimulator), intent(in) :: simulator
        integer(c_int), intent(in), value :: ispec, nspatial, npsi, max_step, use_adaptive_dt
        integer(c_int), intent(in), value :: max_depth, max_samples_per_cell, max_leaves_per_cell, max_expansions
        real(c_double), intent(in) :: spatial_points(3, nspatial), velocity_bounds(6, nspatial), psi_grid(npsi)
        real(c_double), intent(in), value :: dt, refine_threshold_rel, edge_threshold_rel, expand_factor
        integer(c_int), intent(in) :: scout_bins(3)
        real(c_double), intent(out) :: flux(npsi, nspatial), flux_error(npsi, nspatial)
        integer(c_int), intent(out) :: status(nspatial), sample_count(nspatial), terminal_box_count(nspatial)
        integer(c_int), intent(in), optional :: n_threads
        integer :: i
        real(c_double) :: priority_parameters(npsi + 1)

        priority_parameters(1) = qm(ispec)
        priority_parameters(2:) = psi_grid

!$      if (present(n_threads)) call omp_set_num_threads(max(1, int(n_threads)))
        !$omp parallel do schedule(dynamic, 1)
        do i = 1, nspatial
            call evaluate_point(i)
        end do
        !$omp end parallel do
    contains
        subroutine evaluate_point(i)
            integer, intent(in) :: i
            ! Scratch belongs to this spatial point only, bounded by the per-cell
            ! budgets. Never allocate nspatial * max_samples_per_cell here.
            real(c_double), allocatable :: velocities(:, :), probabilities(:), boxes(:, :), pmin(:), pmax(:)
            integer(c_int), allocatable :: si(:), bi(:), depth(:), starts(:), counts(:)
            logical, allocatable :: terminal(:)
            integer(c_int) :: st(1), ns(1), nb(1), actual_ns, actual_nb
            integer :: b, bins(3), first, last
            real(c_double) :: box_flux(npsi), box_error(npsi), priority

            allocate(velocities(3, max_samples_per_cell), probabilities(max_samples_per_cell), &
                si(max_samples_per_cell), boxes(6, max_leaves_per_cell), pmin(max_leaves_per_cell), &
                pmax(max_leaves_per_cell), bi(max_leaves_per_cell), depth(max_leaves_per_cell), &
                starts(max_leaves_per_cell), counts(max_leaves_per_cell), terminal(max_leaves_per_cell))
            call get_probabilities_octree_impl(simulator, ispec, 1_c_int, spatial_points(:, i:i), &
                velocity_bounds(:, i:i), dt, max_step, use_adaptive_dt, scout_bins, max_depth, &
                max_samples_per_cell, max_leaves_per_cell, refine_threshold_rel, edge_threshold_rel, &
                expand_factor, max_expansions, si, velocities, probabilities, bi, boxes, pmin, pmax, &
                depth, starts, counts, st, ns, nb, actual_ns, actual_nb, &
                return_leaf_is_terminal=terminal, complete_partition=.true., parallel_spatial=.false., &
                refinement_priority=moment_priority, refinement_parameters=priority_parameters)
            status(i) = st(1)
            sample_count(i) = ns(1)
            terminal_box_count(i) = count(terminal(:nb(1)))
            flux(:, i) = 0d0
            flux_error(:, i) = 0d0
            do b = 1, nb(1)
                if (.not. terminal(b)) cycle
                bins = [3, 3, 3]
                if (depth(b) == 0) bins = max(2, scout_bins)
                if (counts(b) /= product(bins)) cycle
                first = starts(b) + 1
                last = starts(b) + counts(b)
                call integrate_box(boxes(:, b), bins, velocities(:, first:last), probabilities(first:last), &
                    qm(ispec), psi_grid, box_flux, box_error, priority)
                flux(:, i) = flux(:, i) + box_flux
                flux_error(:, i) = flux_error(:, i) + box_error
            end do
        end subroutine

    end subroutine

    function moment_priority(bounds, bins, velocities, probabilities, parameters) result(priority)
        real(c_double), intent(in) :: bounds(6), velocities(:, :), probabilities(:), parameters(:)
        integer, intent(in) :: bins(3)
        real(c_double) :: priority, box_flux(size(parameters) - 1), box_error(size(parameters) - 1)
        call integrate_box(bounds, bins, velocities, probabilities, parameters(1), parameters(2:), &
            box_flux, box_error, priority)
    end function

    subroutine integrate_box(bounds, bins, velocities, probabilities, charge_mass, psi, flux, error, priority)
        real(c_double), intent(in) :: bounds(6), velocities(:, :), probabilities(:), charge_mass, psi(:)
        integer, intent(in) :: bins(3)
        real(c_double), intent(out) :: flux(size(psi)), error(size(psi)), priority
        real(c_double) :: volume, weight, speed, integrand, low, high, coarse
        integer :: k, j, ix, iy, iz, index(3), axis

        volume = product(bounds(2:6:2) - bounds(1:5:2))
        priority = 0d0
        do k = 1, size(psi)
            low = huge(1d0)
            high = 0d0
            flux(k) = 0d0
            coarse = 0d0
            j = 0
            do iz = 1, bins(3)
                do iy = 1, bins(2)
                    do ix = 1, bins(1)
                        j = j + 1
                        speed = sqrt(sum(velocities(:, j)**2))
                        integrand = 0d0
                        ! Invalid backtraces (-1) contribute zero. The isolated
                        ! v=0 node of the integrable attractive 1/|v| singularity
                        ! is assigned zero; check convergence near the origin.
                        if (speed > 0d0 .and. probabilities(j) > 0d0) then
                            integrand = probabilities(j)*max(0d0, speed - 2d0*charge_mass*psi(k)/speed)
                        end if
                        ! sample_grid includes both faces, vx fastest. Integrate
                        ! each terminal box's OWN lattice by tensor trapezoids.
                        index = [ix, iy, iz]
                        weight = volume/dble(product(bins - 1))
                        do axis = 1, 3
                            if (index(axis) == 1 .or. index(axis) == bins(axis)) weight = weight*0.5d0
                        end do
                        flux(k) = flux(k) + weight*integrand
                        if (all(index == 1 .or. index == bins)) coarse = coarse + volume*integrand/8d0
                        low = min(low, integrand)
                        high = max(high, integrand)
                    end do
                end do
            end do
            error(k) = (high - low)*volume
            ! Coarse (corner-only) vs own-lattice trapezoids estimate quadrature
            ! error for scheduling only. The output remains the own-lattice rule;
            ! flux_error retains the requested sampled range * volume definition.
            priority = max(priority, abs(flux(k) - coarse))
        end do
    end subroutine
end module

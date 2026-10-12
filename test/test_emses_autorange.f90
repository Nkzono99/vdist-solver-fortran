program test_emses_autorange
    use, intrinsic :: iso_c_binding, only: c_loc, c_null_ptr

    use finbound, only: t_BoundaryList, new_BoundaryList

    use m_allcom, only: max_nepl, max_nspec, nemd, nepl, nflag_emit, &
                        path, peth, phixy, phiz, qm, spa, spe, speth, &
                        xmine, xmaxe, ymine, ymaxe, zmine, zmaxe
    use m_emses_autorange, only: estimate_velocity_range_map_impl
    use m_test_helpers, only: assert_close, assert_equal_int, assert_true
    use m_vdsolverf_core

    implicit none

    call test_emission_surface_forward_trace_deposits_velocity_ranges()
    call test_threaded_range_map_matches_serial()
    call test_periodic_seams_do_not_deposit_interior_cells()
    call test_multiple_periodic_crossings_keep_all_deposits()

    print *, "test_emses_autorange: all tests passed."

contains

    subroutine test_multiple_periodic_crossings_keep_all_deposits()
        integer, parameter :: lx = 8, ly = 4, lz = 1, nvel = 17, nseg = 17
        type(t_ESSimulator) :: simulator
        type(t_BoundaryList) :: boundaries
        type(tp_Probability) :: probabilities(0)
        real(8) :: values(6, 0:lx, 0:ly, 0:lz)
        real(8) :: vmin(3, lx, ly, lz), vmax(3, lx, ly, lz)
        real(8) :: weights(lx, ly, lz), confidence(lx, ly, lz), x0, x
        integer :: counts(lx, ly, lz), status(lx, ly, lz), expected(lx)
        integer :: direction, i, ix

        do direction = -1, 1, 2
            call reset_emission_globals()
            path(1) = 0d0
            peth(1) = 0d0
            spa(1) = 0.1d0
            spe(1) = -dble(direction)
            speth(1) = acos(-1d0)/2d0
            nemd(1) = 2
            ymine(1) = 0.4d0
            ymaxe(1) = 0.4d0
            if (direction > 0) then
                xmine(1) = 7.8d0
                xmaxe(1) = 8d0
            else
                xmine(1) = 0d0
                xmaxe(1) = 0.2d0
            end if
            x0 = 0.5d0*(xmine(1) + xmaxe(1))
            values = 0d0
            boundaries = new_BoundaryList()
            simulator = new_ESSimulator(lx, ly, lz, [0, 2, 2], &
                new_VectorFieldGrid(6, lx, ly, lz, values), boundaries, probabilities)
            call estimate_velocity_range_map_impl(simulator, 1, lx, ly, lz, &
                trace_dt=16.2d0, coverage_sigma=1d0, safety_factor=1d0, max_step=1, &
                use_adaptive_dt=0, source_samples_per_cell=1, velocity_sample_mode=0, &
                minimum_count=1, collect_moments=0, show_progress=0, accumulator_cache_size=2, &
                return_vx_min=vmin(1, :, :, :), return_vx_max=vmax(1, :, :, :), &
                return_vy_min=vmin(2, :, :, :), return_vy_max=vmax(2, :, :, :), &
                return_vz_min=vmin(3, :, :, :), return_vz_max=vmax(3, :, :, :), &
                return_count=counts, return_weight_sum=weights, &
                return_mean_v_ptr=c_null_ptr, return_cov_v_ptr=c_null_ptr, &
                return_status=status, return_confidence=confidence)
            expected = 0
            expected(int(x0) + 1) = nvel
            do i = 1, nseg
                x = modulo(x0 + dble(direction)*16.2d0*dble(i)/dble(nseg), 8d0)
                expected(int(x) + 1) = expected(int(x) + 1) + nvel
            end do
            do ix = 1, lx
                call assert_equal_int("multiple wraps retain all physical samples", sum(counts(ix, :, :)), expected(ix))
                call assert_close("multiple wraps retain weight", sum(weights(ix, :, :)), dble(expected(ix)))
            end do
        end do
    end subroutine

    subroutine test_periodic_seams_do_not_deposit_interior_cells()
        integer, parameter :: lx = 8, ly = 1, lz = 8
        type(t_ESSimulator) :: simulator
        type(t_BoundaryList) :: boundaries
        type(tp_Probability) :: probabilities(0)
        real(8) :: values(6, 0:lx, 0:ly, 0:lz)
        real(8) :: vmin(3, lx, ly, lz), vmax(3, lx, ly, lz)
        real(8) :: weights(lx, ly, lz), confidence(lx, ly, lz)
        real(8), target :: mean(3, lx, ly, lz), covariance(3, 3, lx, ly, lz)
        integer :: counts(lx, ly, lz), status(lx, ly, lz)
        integer :: direction, adaptive, varying_field, ix, iz, source_cell, end_cell

        do direction = -1, 1, 2
            do adaptive = 0, 1
                do varying_field = 0, 1
                    call reset_emission_globals()
                    path(1) = 0d0
                    peth(1) = 0d0
                    spa(1) = 0.1d0
                    spe(1) = dble(direction)*sqrt(2d0)
                    speth(1) = -acos(-1d0)/4d0
                    nemd(1) = 2
                    ymine(1) = 0.4d0
                    ymaxe(1) = 0.4d0
                    if (direction > 0) then
                        source_cell = 8
                        end_cell = 1
                        xmine(1) = 7.8d0
                        xmaxe(1) = 8d0
                        zmine(1) = 7.8d0
                        zmaxe(1) = 8d0
                    else
                        source_cell = 1
                        end_cell = 8
                        xmine(1) = 0d0
                        xmaxe(1) = 0.2d0
                        zmine(1) = 0d0
                        zmaxe(1) = 0.2d0
                    end if
                    values = 0d0
                    if (varying_field == 1) then
                        do ix = 0, lx
                            values(1, ix, :, :) = 0.3d0*dble(min(ix, lx - ix))
                        end do
                    end if
                    boundaries = new_BoundaryList()
                    simulator = new_ESSimulator(lx, ly, lz, [0, 2, 0], &
                        new_VectorFieldGrid(6, lx, ly, lz, values), boundaries, probabilities)
                    call estimate_velocity_range_map_impl(simulator, 1, lx, ly, lz, &
                        trace_dt=0.2d0, coverage_sigma=1d0, safety_factor=1d0, max_step=1, &
                        use_adaptive_dt=adaptive, source_samples_per_cell=1, velocity_sample_mode=0, &
                        minimum_count=1, collect_moments=1, show_progress=0, accumulator_cache_size=2, &
                        return_vx_min=vmin(1, :, :, :), return_vx_max=vmax(1, :, :, :), &
                        return_vy_min=vmin(2, :, :, :), return_vy_max=vmax(2, :, :, :), &
                        return_vz_min=vmin(3, :, :, :), return_vz_max=vmax(3, :, :, :), &
                        return_count=counts, return_weight_sum=weights, &
                        return_mean_v_ptr=c_loc(mean), return_cov_v_ptr=c_loc(covariance), &
                        return_status=status, return_confidence=confidence)
                    call assert_true("periodic source receives samples", counts(source_cell, 1, source_cell) > 0)
                    call assert_equal_int("short seam step deposits once per source sample", &
                        counts(end_cell, 1, end_cell), counts(source_cell, 1, source_cell))
                    do iz = 1, lz
                        do ix = 1, lx
                            if (ix == source_cell .and. iz == source_cell) cycle
                            if (ix == end_cell .and. iz == end_cell) cycle
                            call assert_equal_int("seam crossing leaves unvisited cells empty", counts(ix, 1, iz), 0)
                            call assert_close("unvisited cell has zero weight", weights(ix, 1, iz), 0d0)
                            call assert_close("unvisited cell has zero confidence", confidence(ix, 1, iz), 0d0)
                        end do
                    end do
                    call assert_close("source velocity remains local", mean(1, source_cell, 1, source_cell), &
                                      dble(direction))
                    call assert_true("seam velocity stays finite", &
                        abs(mean(1, end_cell, 1, end_cell) - dble(direction)) < 0.01d0)
                end do
            end do
        end do
    end subroutine

    subroutine test_emission_surface_forward_trace_deposits_velocity_ranges()
        integer, parameter :: lx = 2
        integer, parameter :: ly = 1
        integer, parameter :: lz = 1

        type(t_ESSimulator) :: simulator
        type(t_VectorFieldGrid) :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probability_functions(:)
        real(8) :: values(6, 0:lx, 0:ly, 0:lz)
        real(8) :: vx_min(lx, ly, lz), vx_max(lx, ly, lz)
        real(8) :: vy_min(lx, ly, lz), vy_max(lx, ly, lz)
        real(8) :: vz_min(lx, ly, lz), vz_max(lx, ly, lz)
        real(8) :: weight_sum(lx, ly, lz)
        real(8), target :: mean_v(3, lx, ly, lz)
        real(8), target :: cov_v(3, 3, lx, ly, lz)
        real(8) :: confidence(lx, ly, lz)
        integer :: count(lx, ly, lz)
        integer :: status(lx, ly, lz)
        integer :: boundary_conditions(3)

        call reset_emission_globals()

        values = 0d0
        eb = new_VectorFieldGrid(6, lx, ly, lz, values)
        boundaries = new_BoundaryList()
        allocate (probability_functions(1))
        allocate (probability_functions(1)%ref, source=new_ZeroProbability())

        boundary_conditions = [2, 2, 2]
        simulator = new_ESSimulator(lx, ly, lz, boundary_conditions, eb, boundaries, probability_functions)

        call estimate_velocity_range_map_impl( &
            simulator=simulator, &
            ispec=1, &
            lx=lx, ly=ly, lz=lz, &
            trace_dt=0.5d0, &
            coverage_sigma=1d0, &
            safety_factor=1d0, &
            max_step=5, &
            use_adaptive_dt=0, &
            source_samples_per_cell=1, &
            velocity_sample_mode=0, &
            minimum_count=1, &
            collect_moments=1, &
            show_progress=0, &
            accumulator_cache_size=2, &
            return_vx_min=vx_min, &
            return_vx_max=vx_max, &
            return_vy_min=vy_min, &
            return_vy_max=vy_max, &
            return_vz_min=vz_min, &
            return_vz_max=vz_max, &
            return_count=count, &
            return_weight_sum=weight_sum, &
            return_mean_v_ptr=c_loc(mean_v), &
            return_cov_v_ptr=c_loc(cov_v), &
            return_status=status, &
            return_confidence=confidence)

        call assert_true("source cell receives envelope samples", count(1, 1, 1) > 0)
        call assert_true("downstream cell receives forward-traced samples", count(2, 1, 1) > 0)
        call assert_true("x velocity range is positive in source cell", vx_max(1, 1, 1) > 0d0)
        call assert_true("x velocity min/max are ordered", vx_min(1, 1, 1) <= vx_max(1, 1, 1))
        call assert_true("mean velocity diagnostic is populated", mean_v(1, 1, 1, 1) > 0d0)
    end subroutine

    subroutine reset_emission_globals()
        qm = 0d0
        path = 0d0
        peth = 0d0
        spa = 0d0
        spe = 0d0
        speth = 0d0
        phiz = 0d0
        phixy = 0d0
        nflag_emit = 0
        nepl = 0
        nemd = 0
        xmine = 0d0
        xmaxe = 0d0
        ymine = 0d0
        ymaxe = 0d0
        zmine = 0d0
        zmaxe = 0d0

        qm(1) = 1d0
        path(1) = 0.2d0
        peth(1) = 0.1d0
        spa(1) = 1d0
        nflag_emit(1) = 1
        nepl(1) = 1
        nemd(1) = 1
        xmine(1) = 0d0
        xmaxe(1) = 0d0
        ymine(1) = 0d0
        ymaxe(1) = 1d0
        zmine(1) = 0d0
        zmaxe(1) = 1d0
    end subroutine

    subroutine test_threaded_range_map_matches_serial()
        integer, parameter :: lx = 4
        integer, parameter :: ly = 3
        integer, parameter :: lz = 2

        real(8) :: vx_min_1(lx, ly, lz), vx_max_1(lx, ly, lz)
        real(8) :: vy_min_1(lx, ly, lz), vy_max_1(lx, ly, lz)
        real(8) :: vz_min_1(lx, ly, lz), vz_max_1(lx, ly, lz)
        real(8) :: vx_min_4(lx, ly, lz), vx_max_4(lx, ly, lz)
        real(8) :: vy_min_4(lx, ly, lz), vy_max_4(lx, ly, lz)
        real(8) :: vz_min_4(lx, ly, lz), vz_max_4(lx, ly, lz)
        integer :: count_1(lx, ly, lz), count_4(lx, ly, lz)
        integer :: ix, iy, iz

        call run_autorange(lx, ly, lz, 1, vx_min_1, vx_max_1, vy_min_1, vy_max_1, &
                           vz_min_1, vz_max_1, count_1)
        call run_autorange(lx, ly, lz, 4, vx_min_4, vx_max_4, vy_min_4, vy_max_4, &
                           vz_min_4, vz_max_4, count_4)

        call assert_equal_int("threaded total count matches serial", sum(count_4), sum(count_1))

        do iz = 1, lz
            do iy = 1, ly
                do ix = 1, lx
                    call assert_equal_int("threaded cell count matches serial", &
                                          count_4(ix, iy, iz), count_1(ix, iy, iz))
                    if (count_1(ix, iy, iz) <= 0) cycle

                    call assert_close("threaded vx_min matches serial", vx_min_4(ix, iy, iz), &
                                      vx_min_1(ix, iy, iz), 1d-12)
                    call assert_close("threaded vx_max matches serial", vx_max_4(ix, iy, iz), &
                                      vx_max_1(ix, iy, iz), 1d-12)
                    call assert_close("threaded vy_min matches serial", vy_min_4(ix, iy, iz), &
                                      vy_min_1(ix, iy, iz), 1d-12)
                    call assert_close("threaded vy_max matches serial", vy_max_4(ix, iy, iz), &
                                      vy_max_1(ix, iy, iz), 1d-12)
                    call assert_close("threaded vz_min matches serial", vz_min_4(ix, iy, iz), &
                                      vz_min_1(ix, iy, iz), 1d-12)
                    call assert_close("threaded vz_max matches serial", vz_max_4(ix, iy, iz), &
                                      vz_max_1(ix, iy, iz), 1d-12)
                end do
            end do
        end do
    end subroutine

    subroutine run_autorange(lx, ly, lz, n_threads, vx_min, vx_max, vy_min, vy_max, &
                             vz_min, vz_max, count)
        integer, intent(in) :: lx, ly, lz
        integer, intent(in) :: n_threads
        real(8), intent(out) :: vx_min(lx, ly, lz), vx_max(lx, ly, lz)
        real(8), intent(out) :: vy_min(lx, ly, lz), vy_max(lx, ly, lz)
        real(8), intent(out) :: vz_min(lx, ly, lz), vz_max(lx, ly, lz)
        integer, intent(out) :: count(lx, ly, lz)

        type(t_ESSimulator) :: simulator
        type(t_VectorFieldGrid) :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probability_functions(:)
        real(8) :: values(6, 0:lx, 0:ly, 0:lz)
        real(8) :: weight_sum(lx, ly, lz)
        real(8) :: confidence(lx, ly, lz)
        integer :: status(lx, ly, lz)
        integer :: boundary_conditions(3)

        call reset_emission_globals()
        ymaxe(1) = dble(ly)
        zmaxe(1) = dble(lz)

        values = 0d0
        eb = new_VectorFieldGrid(6, lx, ly, lz, values)
        boundaries = new_BoundaryList()
        allocate (probability_functions(1))
        allocate (probability_functions(1)%ref, source=new_ZeroProbability())

        boundary_conditions = [2, 2, 2]
        simulator = new_ESSimulator(lx, ly, lz, boundary_conditions, eb, boundaries, probability_functions)

        call estimate_velocity_range_map_impl( &
            simulator=simulator, &
            ispec=1, &
            lx=lx, ly=ly, lz=lz, &
            trace_dt=0.5d0, &
            coverage_sigma=1d0, &
            safety_factor=1d0, &
            max_step=8, &
            use_adaptive_dt=0, &
            source_samples_per_cell=1, &
            velocity_sample_mode=0, &
            minimum_count=1, &
            collect_moments=0, &
            show_progress=0, &
            accumulator_cache_size=2, &
            n_threads=n_threads, &
            return_vx_min=vx_min, &
            return_vx_max=vx_max, &
            return_vy_min=vy_min, &
            return_vy_max=vy_max, &
            return_vz_min=vz_min, &
            return_vz_max=vz_max, &
            return_count=count, &
            return_weight_sum=weight_sum, &
            return_mean_v_ptr=c_null_ptr, &
            return_cov_v_ptr=c_null_ptr, &
            return_status=status, &
            return_confidence=confidence)
    end subroutine

end program

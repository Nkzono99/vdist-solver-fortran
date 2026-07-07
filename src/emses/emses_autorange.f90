module m_emses_autorange
    !! Deterministic source-envelope tracing for per-cell velocity ranges.

    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_f_pointer, c_int, c_int64_t, c_ptr
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value

!$  use omp_lib, only: omp_destroy_lock, omp_get_thread_num, omp_init_lock, omp_lock_kind, &
!$                     omp_set_lock, omp_set_num_threads, omp_unset_lock

    use forbear, only: bar_object
    use finbound, only: t_CollisionRecord

    use m_vdsolverf_core
    use m_allcom, only: emission_vdri_vector, emission_vth_vector, &
                        max_nepl, nepl, nemd, nflag_emit, npbnd, qm, &
                        vdri_vector, vth_vector, xmine, xmaxe, ymine, ymaxe, &
                        zmine, zmaxe, zssurf

    implicit none

    private
    public estimate_velocity_range_map_impl

    integer(c_int), parameter :: STATUS_OK = 0
    integer(c_int), parameter :: STATUS_LOW_COUNT = 1
    integer(c_int), parameter :: STATUS_FALLBACK = 2

contains

    subroutine estimate_velocity_range_map_impl( &
        simulator, &
        ispec, &
        lx, ly, lz, &
        trace_dt, &
        coverage_sigma, &
        safety_factor, &
        max_step, &
        use_adaptive_dt, &
        source_samples_per_cell, &
        velocity_sample_mode, &
        minimum_count, &
        collect_moments, &
        show_progress, &
        accumulator_cache_size, &
        return_vx_min, &
        return_vx_max, &
        return_vy_min, &
        return_vy_max, &
        return_vz_min, &
        return_vz_max, &
        return_count, &
        return_weight_sum, &
        return_mean_v_ptr, &
        return_cov_v_ptr, &
        return_status, &
        return_confidence, &
        n_threads)
        type(t_ESSimulator), intent(in) :: simulator
        integer(c_int), value, intent(in) :: ispec
        integer(c_int), value, intent(in) :: lx
        integer(c_int), value, intent(in) :: ly
        integer(c_int), value, intent(in) :: lz
        real(c_double), value, intent(in) :: trace_dt
        real(c_double), value, intent(in) :: coverage_sigma
        real(c_double), value, intent(in) :: safety_factor
        integer(c_int), value, intent(in) :: max_step
        integer(c_int), value, intent(in) :: use_adaptive_dt
        integer(c_int), value, intent(in) :: source_samples_per_cell
        integer(c_int), value, intent(in) :: velocity_sample_mode
        integer(c_int), value, intent(in) :: minimum_count
        integer(c_int), value, intent(in) :: collect_moments
        integer(c_int), value, optional, intent(in) :: show_progress
        integer(c_int), value, optional, intent(in) :: accumulator_cache_size
        real(c_double), intent(out) :: return_vx_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vx_max(lx, ly, lz)
        real(c_double), intent(out) :: return_vy_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vy_max(lx, ly, lz)
        real(c_double), intent(out) :: return_vz_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vz_max(lx, ly, lz)
        integer(c_int), intent(out) :: return_count(lx, ly, lz)
        real(c_double), intent(out) :: return_weight_sum(lx, ly, lz)
        type(c_ptr), value, intent(in) :: return_mean_v_ptr
        type(c_ptr), value, intent(in) :: return_cov_v_ptr
        integer(c_int), intent(out) :: return_status(lx, ly, lz)
        real(c_double), intent(out) :: return_confidence(lx, ly, lz)
        integer(c_int), value, optional, intent(in) :: n_threads

        real(c_double) :: nan_value
        integer :: active_n_threads
        integer :: cache_entry_limit
        integer :: cache_hash_capacity
        logical :: collect_moments_enabled
        logical :: show_progress_enabled
        type(bar_object) :: bar
        integer :: progress_total_tasks
        integer :: progress_completed_tasks
        integer :: progress_last_percent
        integer(c_int64_t), allocatable :: cache_cell_index(:, :)
        integer(c_int), allocatable :: cache_used(:)
        integer(c_int), allocatable :: cache_count(:, :)
        real(c_double), allocatable :: cache_vx_min(:, :)
        real(c_double), allocatable :: cache_vx_max(:, :)
        real(c_double), allocatable :: cache_vy_min(:, :)
        real(c_double), allocatable :: cache_vy_max(:, :)
        real(c_double), allocatable :: cache_vz_min(:, :)
        real(c_double), allocatable :: cache_vz_max(:, :)
        real(c_double), allocatable :: cache_weight_sum(:, :)
        real(c_double), allocatable :: cache_mean_v(:, :, :)
        real(c_double), allocatable :: cache_cov_v(:, :, :, :)
        real(c_double), pointer :: return_mean_v(:, :, :, :)
        real(c_double), pointer :: return_cov_v(:, :, :, :, :)
!$      integer(omp_lock_kind), allocatable :: flush_locks(:)
        integer, parameter :: FLUSH_LOCK_COUNT = 4096

        nan_value = ieee_value(0d0, ieee_quiet_nan)
        active_n_threads = 1
        if (present(n_threads)) then
            active_n_threads = max(1, int(n_threads))
        end if
        collect_moments_enabled = collect_moments == 1
        if (collect_moments_enabled) then
            collect_moments_enabled = c_associated(return_mean_v_ptr) .and. c_associated(return_cov_v_ptr)
        end if
        show_progress_enabled = .false.
        if (present(show_progress)) then
            show_progress_enabled = show_progress == 1
        end if
        cache_entry_limit = 20000
        if (present(accumulator_cache_size)) then
            cache_entry_limit = max(1, int(accumulator_cache_size))
        end if
        cache_hash_capacity = max(8, 2*cache_entry_limit + 1)
        progress_total_tasks = 0
        progress_completed_tasks = 0
        progress_last_percent = -1

!$      call omp_set_num_threads(active_n_threads)

        if (collect_moments_enabled) then
            call c_f_pointer(return_mean_v_ptr, return_mean_v, [3, lx, ly, lz])
            call c_f_pointer(return_cov_v_ptr, return_cov_v, [3, 3, lx, ly, lz])
        end if

        call initialize_output_accumulators
        call initialize_cache_accumulators
        call initialize_flush_locks

        call initialize_progress_bar
        call trace_external_boundary_sources
        call trace_emission_surface_sources
        call flush_all_caches
        call finish_progress_bar
        call finalize_cells
        call destroy_flush_locks

    contains

        subroutine initialize_progress_bar
            if (.not. show_progress_enabled) then
                return
            end if

            progress_total_tasks = count_external_boundary_source_tasks() &
                                   + count_emission_surface_source_tasks()
            if (progress_total_tasks <= 0) then
                show_progress_enabled = .false.
                return
            end if

            call bar%initialize(filled_char_string='+', &
                                prefix_string='progress |', &
                                suffix_string='| ', &
                                add_progress_percent=.true.)
            call bar%start
        end subroutine

        subroutine finish_progress_bar
            if (.not. show_progress_enabled) then
                return
            end if

            call bar%update(current=1d0)
            call bar%destroy
        end subroutine

        subroutine record_progress_task
            real(c_double) :: progress_fraction
            integer :: progress_percent

            if (.not. show_progress_enabled) then
                return
            end if

!$omp      critical(emses_autorange_progress)
            progress_completed_tasks = progress_completed_tasks + 1
            progress_fraction = min(0.99d0, dble(progress_completed_tasks)/dble(progress_total_tasks))
            progress_percent = int(100d0*progress_fraction)
            if (progress_percent > progress_last_percent) then
                call bar%update(current=progress_fraction)
                progress_last_percent = progress_percent
            end if
!$omp      end critical(emses_autorange_progress)
        end subroutine

        integer function count_external_boundary_source_tasks() result(ret)
            integer :: axis
            real(c_double) :: patch_min(3), patch_max(3)
            real(c_double) :: vmean(3), vthermal(3)

            ret = 0
            if (nflag_emit(ispec) /= 0) then
                return
            end if

            vmean = vdri_vector(ispec)
            vthermal = abs(vth_vector(ispec))

            do axis = 1, 3
                if (npbnd(axis, ispec) /= 2) then
                    cycle
                end if

                patch_min = [0d0, 0d0, 0d0]
                patch_max = [dble(lx), dble(ly), dble(lz)]

                patch_min(axis) = 0d0
                patch_max(axis) = 0d0
                if (.not. (axis == 3 .and. zssurf >= 0d0)) then
                    ret = ret + count_source_patch_tasks(axis, 1, patch_min, patch_max, vmean, vthermal)
                end if

                patch_min = [0d0, 0d0, 0d0]
                patch_max = [dble(lx), dble(ly), dble(lz)]
                patch_min(axis) = dble(dim_size(axis))
                patch_max(axis) = patch_min(axis)
                ret = ret + count_source_patch_tasks(axis, -1, patch_min, patch_max, vmean, vthermal)
            end do
        end function

        integer function count_emission_surface_source_tasks() result(ret)
            integer :: iepl
            integer :: iepl_start
            integer :: iepl_end
            integer :: normal_axis
            integer :: normal_sign
            real(c_double) :: patch_min(3), patch_max(3)
            real(c_double) :: vmean(3), vthermal(3)

            ret = 0
            if (nepl(ispec) == 0) then
                return
            end if

            if (ispec == 1) then
                iepl_start = 1
            else
                iepl_start = sum(nepl(1:ispec - 1)) + 1
            end if
            iepl_end = sum(nepl(1:ispec))

            do iepl = iepl_start, iepl_end
                normal_axis = abs(nemd(iepl))
                if (normal_axis < 1 .or. normal_axis > 3) then
                    cycle
                end if

                if (nemd(iepl) > 0) then
                    normal_sign = 1
                else
                    normal_sign = -1
                end if

                patch_min = [xmine(iepl), ymine(iepl), zmine(iepl)]
                patch_max = [xmaxe(iepl), ymaxe(iepl), zmaxe(iepl)]
                vmean = emission_vdri_vector(ispec, iepl)
                vthermal = abs(emission_vth_vector(ispec, iepl))

                ret = ret + count_source_patch_tasks(normal_axis, normal_sign, patch_min, patch_max, vmean, vthermal)
            end do
        end function

        integer function count_source_patch_tasks(normal_axis, normal_sign, patch_min, patch_max, vmean, vthermal) result(ret)
            integer, intent(in) :: normal_axis
            integer, intent(in) :: normal_sign
            real(c_double), intent(in) :: patch_min(3), patch_max(3)
            real(c_double), intent(in) :: vmean(3), vthermal(3)

            integer :: tangent1, tangent2
            integer :: c1_start, c1_end, c2_start, c2_end
            integer :: n1, n2
            real(c_double) :: velocities(3, 32)
            integer :: nvel

            ret = 0
            call tangent_axes(normal_axis, tangent1, tangent2)
            call build_velocity_support(vmean, vthermal, normal_axis, normal_sign, velocities, nvel)
            if (nvel <= 0) then
                return
            end if

            c1_start = first_cell(patch_min(tangent1), dim_size(tangent1))
            c1_end = last_cell(patch_max(tangent1), dim_size(tangent1))
            c2_start = first_cell(patch_min(tangent2), dim_size(tangent2))
            c2_end = last_cell(patch_max(tangent2), dim_size(tangent2))

            n1 = c1_end - c1_start + 1
            n2 = c2_end - c2_start + 1
            if (n1 <= 0 .or. n2 <= 0) then
                return
            end if

            ret = n1*n2
        end function

        subroutine initialize_output_accumulators
            return_vx_min = huge(1d0)
            return_vy_min = huge(1d0)
            return_vz_min = huge(1d0)
            return_vx_max = -huge(1d0)
            return_vy_max = -huge(1d0)
            return_vz_max = -huge(1d0)
            return_count = 0
            return_weight_sum = 0d0

            if (collect_moments_enabled) then
                return_mean_v = 0d0
                return_cov_v = 0d0
            end if
            return_status = STATUS_FALLBACK
            return_confidence = 0d0
        end subroutine

        subroutine initialize_cache_accumulators
            allocate (cache_cell_index(cache_hash_capacity, active_n_threads))
            allocate (cache_used(active_n_threads))
            allocate (cache_count(cache_hash_capacity, active_n_threads))
            allocate (cache_weight_sum(cache_hash_capacity, active_n_threads))
            allocate (cache_vx_min(cache_hash_capacity, active_n_threads))
            allocate (cache_vx_max(cache_hash_capacity, active_n_threads))
            allocate (cache_vy_min(cache_hash_capacity, active_n_threads))
            allocate (cache_vy_max(cache_hash_capacity, active_n_threads))
            allocate (cache_vz_min(cache_hash_capacity, active_n_threads))
            allocate (cache_vz_max(cache_hash_capacity, active_n_threads))

            cache_cell_index = 0_c_int64_t
            cache_used = 0
            cache_count = 0
            cache_weight_sum = 0d0

            if (collect_moments_enabled) then
                allocate (cache_mean_v(3, cache_hash_capacity, active_n_threads))
                allocate (cache_cov_v(3, 3, cache_hash_capacity, active_n_threads))
                cache_mean_v = 0d0
                cache_cov_v = 0d0
            end if
        end subroutine

        subroutine initialize_flush_locks
            integer :: ilock

!$          allocate (flush_locks(FLUSH_LOCK_COUNT))
!$          do ilock = 1, FLUSH_LOCK_COUNT
!$              call omp_init_lock(flush_locks(ilock))
!$          end do
        end subroutine

        subroutine destroy_flush_locks
            integer :: ilock

!$          if (allocated(flush_locks)) then
!$              do ilock = 1, size(flush_locks)
!$                  call omp_destroy_lock(flush_locks(ilock))
!$              end do
!$              deallocate (flush_locks)
!$          end if
        end subroutine

        subroutine flush_all_caches
            integer :: ithread

            do ithread = 1, active_n_threads
                call flush_worker_cache(ithread)
            end do
        end subroutine

        subroutine flush_worker_cache(ithread)
            integer, intent(in) :: ithread

            integer :: slot
            integer :: ix, iy, iz
            integer :: lock_id
            integer(c_int64_t) :: cell_index

            if (cache_used(ithread) <= 0) then
                return
            end if

            do slot = 1, cache_hash_capacity
                cell_index = cache_cell_index(slot, ithread)
                if (cell_index == 0_c_int64_t) cycle

                call decode_cell_index(cell_index, ix, iy, iz)
                lock_id = flush_lock_index(cell_index)
!$              call omp_set_lock(flush_locks(lock_id))
                call reduce_cache_slot(slot, ithread, ix, iy, iz)
!$              call omp_unset_lock(flush_locks(lock_id))
            end do

            cache_cell_index(:, ithread) = 0_c_int64_t
            cache_used(ithread) = 0
        end subroutine

        subroutine reduce_cache_slot(slot, ithread, ix, iy, iz)
            integer, intent(in) :: slot
            integer, intent(in) :: ithread
            integer, intent(in) :: ix, iy, iz

            return_vx_min(ix, iy, iz) = min(return_vx_min(ix, iy, iz), cache_vx_min(slot, ithread))
            return_vx_max(ix, iy, iz) = max(return_vx_max(ix, iy, iz), cache_vx_max(slot, ithread))
            return_vy_min(ix, iy, iz) = min(return_vy_min(ix, iy, iz), cache_vy_min(slot, ithread))
            return_vy_max(ix, iy, iz) = max(return_vy_max(ix, iy, iz), cache_vy_max(slot, ithread))
            return_vz_min(ix, iy, iz) = min(return_vz_min(ix, iy, iz), cache_vz_min(slot, ithread))
            return_vz_max(ix, iy, iz) = max(return_vz_max(ix, iy, iz), cache_vz_max(slot, ithread))
            return_count(ix, iy, iz) = return_count(ix, iy, iz) + cache_count(slot, ithread)
            return_weight_sum(ix, iy, iz) = return_weight_sum(ix, iy, iz) + cache_weight_sum(slot, ithread)

            if (collect_moments_enabled) then
                return_mean_v(:, ix, iy, iz) = return_mean_v(:, ix, iy, iz) + cache_mean_v(:, slot, ithread)
                return_cov_v(:, :, ix, iy, iz) = return_cov_v(:, :, ix, iy, iz) + cache_cov_v(:, :, slot, ithread)
            end if
        end subroutine

        subroutine trace_external_boundary_sources
            integer :: axis
            real(c_double) :: patch_min(3), patch_max(3)
            real(c_double) :: vmean(3), vthermal(3)

            if (nflag_emit(ispec) /= 0) then
                return
            end if

            vmean = vdri_vector(ispec)
            vthermal = abs(vth_vector(ispec))

            do axis = 1, 3
                if (npbnd(axis, ispec) /= 2) then
                    cycle
                end if

                patch_min = [0d0, 0d0, 0d0]
                patch_max = [dble(lx), dble(ly), dble(lz)]

                patch_min(axis) = 0d0
                patch_max(axis) = 0d0
                if (.not. (axis == 3 .and. zssurf >= 0d0)) then
                    call trace_source_patch(axis, 1, patch_min, patch_max, vmean, vthermal)
                end if

                patch_min = [0d0, 0d0, 0d0]
                patch_max = [dble(lx), dble(ly), dble(lz)]
                patch_min(axis) = dble(dim_size(axis))
                patch_max(axis) = patch_min(axis)
                call trace_source_patch(axis, -1, patch_min, patch_max, vmean, vthermal)
            end do
        end subroutine

        subroutine trace_emission_surface_sources
            integer :: iepl
            integer :: iepl_start
            integer :: iepl_end
            integer :: normal_axis
            integer :: normal_sign
            real(c_double) :: patch_min(3), patch_max(3)
            real(c_double) :: vmean(3), vthermal(3)

            if (nepl(ispec) == 0) then
                return
            end if

            if (ispec == 1) then
                iepl_start = 1
            else
                iepl_start = sum(nepl(1:ispec - 1)) + 1
            end if
            iepl_end = sum(nepl(1:ispec))

            do iepl = iepl_start, iepl_end
                normal_axis = abs(nemd(iepl))
                if (normal_axis < 1 .or. normal_axis > 3) then
                    cycle
                end if

                if (nemd(iepl) > 0) then
                    normal_sign = 1
                else
                    normal_sign = -1
                end if

                patch_min = [xmine(iepl), ymine(iepl), zmine(iepl)]
                patch_max = [xmaxe(iepl), ymaxe(iepl), zmaxe(iepl)]
                vmean = emission_vdri_vector(ispec, iepl)
                vthermal = abs(emission_vth_vector(ispec, iepl))

                call trace_source_patch(normal_axis, normal_sign, patch_min, patch_max, vmean, vthermal)
            end do
        end subroutine

        subroutine trace_source_patch(normal_axis, normal_sign, patch_min, patch_max, vmean, vthermal)
            integer, intent(in) :: normal_axis
            integer, intent(in) :: normal_sign
            real(c_double), intent(in) :: patch_min(3), patch_max(3)
            real(c_double), intent(in) :: vmean(3), vthermal(3)

            integer :: tangent1, tangent2
            integer :: c1, c2, s1, s2
            integer :: c1_start, c1_end, c2_start, c2_end
            integer :: n1, n2, itask
            integer :: nsub
            real(c_double) :: low1, high1, low2, high2
            real(c_double) :: position(3)
            real(c_double) :: velocities(3, 32)
            integer :: ivel, nvel

            call tangent_axes(normal_axis, tangent1, tangent2)

            nsub = max(1, source_samples_per_cell)
            c1_start = first_cell(patch_min(tangent1), dim_size(tangent1))
            c1_end = last_cell(patch_max(tangent1), dim_size(tangent1))
            c2_start = first_cell(patch_min(tangent2), dim_size(tangent2))
            c2_end = last_cell(patch_max(tangent2), dim_size(tangent2))

            call build_velocity_support(vmean, vthermal, normal_axis, normal_sign, velocities, nvel)
            if (nvel <= 0) then
                return
            end if

            n1 = c1_end - c1_start + 1
            n2 = c2_end - c2_start + 1
            if (n1 <= 0 .or. n2 <= 0) then
                return
            end if

            if (collect_moments_enabled) then
                !$omp parallel do schedule(static) &
                !$omp private(itask, c2, c1, s2, s1, low1, high1, low2, high2, position, ivel)
                do itask = 0, n1*n2 - 1
                    c2 = c2_start + itask/n1
                    c1 = c1_start + mod(itask, n1)

                    low2 = max(patch_min(tangent2), dble(c2))
                    high2 = min(patch_max(tangent2), dble(c2 + 1))
                    if (high2 <= low2) then
                        call record_progress_task
                        cycle
                    end if

                    low1 = max(patch_min(tangent1), dble(c1))
                    high1 = min(patch_max(tangent1), dble(c1 + 1))
                    if (high1 <= low1) then
                        call record_progress_task
                        cycle
                    end if

                    do s2 = 1, nsub
                        do s1 = 1, nsub
                            position = 0d0
                            position(tangent1) = low1 + (dble(s1) - 0.5d0)*(high1 - low1)/dble(nsub)
                            position(tangent2) = low2 + (dble(s2) - 0.5d0)*(high2 - low2)/dble(nsub)
                            position(normal_axis) = source_normal_position(normal_axis, normal_sign, patch_min, patch_max)

                            do ivel = 1, nvel
                                call trace_support_particle(position, velocities(:, ivel), 1d0)
                            end do
                        end do
                    end do
                    call record_progress_task
                end do
                !$omp end parallel do
            else
                !$omp parallel do schedule(static) &
                !$omp private(itask, c2, c1, s2, s1, low1, high1, low2, high2, position, ivel)
                do itask = 0, n1*n2 - 1
                    c2 = c2_start + itask/n1
                    c1 = c1_start + mod(itask, n1)

                    low2 = max(patch_min(tangent2), dble(c2))
                    high2 = min(patch_max(tangent2), dble(c2 + 1))
                    if (high2 <= low2) then
                        call record_progress_task
                        cycle
                    end if

                    low1 = max(patch_min(tangent1), dble(c1))
                    high1 = min(patch_max(tangent1), dble(c1 + 1))
                    if (high1 <= low1) then
                        call record_progress_task
                        cycle
                    end if

                    do s2 = 1, nsub
                        do s1 = 1, nsub
                            position = 0d0
                            position(tangent1) = low1 + (dble(s1) - 0.5d0)*(high1 - low1)/dble(nsub)
                            position(tangent2) = low2 + (dble(s2) - 0.5d0)*(high2 - low2)/dble(nsub)
                            position(normal_axis) = source_normal_position(normal_axis, normal_sign, patch_min, patch_max)

                            do ivel = 1, nvel
                                call trace_support_particle(position, velocities(:, ivel), 1d0)
                            end do
                        end do
                    end do
                    call record_progress_task
                end do
                !$omp end parallel do
            end if
        end subroutine

        subroutine trace_support_particle(position, velocity, weight)
            real(c_double), intent(in) :: position(3)
            real(c_double), intent(in) :: velocity(3)
            real(c_double), intent(in) :: weight

            type(t_Particle) :: pcl
            type(t_Particle) :: pcl_prev
            type(t_CollisionRecord) :: collision
            real(c_double) :: step_dt
            real(c_double) :: speed
            integer :: istep

            pcl = new_Particle(qm(ispec), position, velocity)
            call deposit_point(pcl%position, pcl%velocity, weight)

            do istep = 1, max_step
                pcl_prev = pcl
                step_dt = -abs(trace_dt)
                if (use_adaptive_dt == 1) then
                    speed = sqrt(sum(pcl%velocity*pcl%velocity))
                    if (speed > tiny(1d0)) then
                        step_dt = -abs(trace_dt)/speed
                    end if
                end if

                pcl = simulator%update(pcl, step_dt, collision)
                call deposit_segment(pcl_prev%position, pcl%position, pcl_prev%velocity, pcl%velocity, weight)

                if (collision%is_collided) then
                    exit
                end if
            end do
        end subroutine

        subroutine deposit_segment(position0, position1, velocity0, velocity1, weight)
            real(c_double), intent(in) :: position0(3)
            real(c_double), intent(in) :: position1(3)
            real(c_double), intent(in) :: velocity0(3)
            real(c_double), intent(in) :: velocity1(3)
            real(c_double), intent(in) :: weight

            integer :: nseg
            integer :: i
            real(c_double) :: ratio
            real(c_double) :: position(3), velocity(3)

            nseg = max(1, ceiling(maxval(abs(position1 - position0))))
            do i = 1, nseg
                ratio = dble(i)/dble(nseg)
                position = position0*(1d0 - ratio) + position1*ratio
                velocity = velocity0*(1d0 - ratio) + velocity1*ratio
                call deposit_point(position, velocity, weight)
            end do
        end subroutine

        subroutine deposit_point(position, velocity, weight)
            real(c_double), intent(in) :: position(3)
            real(c_double), intent(in) :: velocity(3)
            real(c_double), intent(in) :: weight

            integer :: ix, iy, iz
            integer :: ithread
            integer(c_int64_t) :: cell_index

            if (position(1) < 0d0 .or. position(1) > dble(lx)) return
            if (position(2) < 0d0 .or. position(2) > dble(ly)) return
            if (position(3) < 0d0 .or. position(3) > dble(lz)) return

            ix = min(max(int(position(1)) + 1, 1), lx)
            iy = min(max(int(position(2)) + 1, 1), ly)
            iz = min(max(int(position(3)) + 1, 1), lz)

            ithread = worker_index()
            cell_index = encode_cell_index(ix, iy, iz)

            call deposit_to_cache(ithread, cell_index, velocity, weight)
        end subroutine

        subroutine deposit_to_cache(ithread, cell_index, velocity, weight)
            integer, intent(in) :: ithread
            integer(c_int64_t), intent(in) :: cell_index
            real(c_double), intent(in) :: velocity(3)
            real(c_double), intent(in) :: weight

            integer :: slot
            logical :: found

            call find_cache_slot(ithread, cell_index, slot, found)
            if (.not. found .and. cache_used(ithread) >= cache_entry_limit) then
                call flush_worker_cache(ithread)
                call find_cache_slot(ithread, cell_index, slot, found)
            end if

            if (found) then
                call update_cache_slot(slot, ithread, velocity, weight)
            else
                call initialize_cache_slot(slot, ithread, cell_index, velocity, weight)
            end if
        end subroutine

        subroutine find_cache_slot(ithread, cell_index, slot, found)
            integer, intent(in) :: ithread
            integer(c_int64_t), intent(in) :: cell_index
            integer, intent(out) :: slot
            logical, intent(out) :: found

            integer :: probe
            integer :: start_slot

            start_slot = cache_hash_slot(cell_index)
            found = .false.
            slot = start_slot

            do probe = 0, cache_hash_capacity - 1
                slot = 1 + mod(start_slot - 1 + probe, cache_hash_capacity)
                if (cache_cell_index(slot, ithread) == cell_index) then
                    found = .true.
                    return
                end if
                if (cache_cell_index(slot, ithread) == 0_c_int64_t) then
                    return
                end if
            end do

            call flush_worker_cache(ithread)
            slot = cache_hash_slot(cell_index)
            found = .false.
        end subroutine

        subroutine initialize_cache_slot(slot, ithread, cell_index, velocity, weight)
            integer, intent(in) :: slot
            integer, intent(in) :: ithread
            integer(c_int64_t), intent(in) :: cell_index
            real(c_double), intent(in) :: velocity(3)
            real(c_double), intent(in) :: weight

            integer :: i, j

            cache_cell_index(slot, ithread) = cell_index
            cache_used(ithread) = cache_used(ithread) + 1
            cache_vx_min(slot, ithread) = velocity(1)
            cache_vx_max(slot, ithread) = velocity(1)
            cache_vy_min(slot, ithread) = velocity(2)
            cache_vy_max(slot, ithread) = velocity(2)
            cache_vz_min(slot, ithread) = velocity(3)
            cache_vz_max(slot, ithread) = velocity(3)
            cache_count(slot, ithread) = 1
            cache_weight_sum(slot, ithread) = weight

            if (collect_moments_enabled) then
                do i = 1, 3
                    cache_mean_v(i, slot, ithread) = weight*velocity(i)
                end do
                do j = 1, 3
                    do i = 1, 3
                        cache_cov_v(i, j, slot, ithread) = weight*velocity(i)*velocity(j)
                    end do
                end do
            end if
        end subroutine

        subroutine update_cache_slot(slot, ithread, velocity, weight)
            integer, intent(in) :: slot
            integer, intent(in) :: ithread
            real(c_double), intent(in) :: velocity(3)
            real(c_double), intent(in) :: weight

            integer :: i, j

            cache_vx_min(slot, ithread) = min(cache_vx_min(slot, ithread), velocity(1))
            cache_vx_max(slot, ithread) = max(cache_vx_max(slot, ithread), velocity(1))
            cache_vy_min(slot, ithread) = min(cache_vy_min(slot, ithread), velocity(2))
            cache_vy_max(slot, ithread) = max(cache_vy_max(slot, ithread), velocity(2))
            cache_vz_min(slot, ithread) = min(cache_vz_min(slot, ithread), velocity(3))
            cache_vz_max(slot, ithread) = max(cache_vz_max(slot, ithread), velocity(3))
            cache_count(slot, ithread) = cache_count(slot, ithread) + 1
            cache_weight_sum(slot, ithread) = cache_weight_sum(slot, ithread) + weight

            if (collect_moments_enabled) then
                do i = 1, 3
                    cache_mean_v(i, slot, ithread) = cache_mean_v(i, slot, ithread) + weight*velocity(i)
                end do
                do j = 1, 3
                    do i = 1, 3
                        cache_cov_v(i, j, slot, ithread) = cache_cov_v(i, j, slot, ithread) &
                                                            + weight*velocity(i)*velocity(j)
                    end do
                end do
            end if
        end subroutine

        integer(c_int64_t) function encode_cell_index(ix, iy, iz) result(ret)
            integer, intent(in) :: ix, iy, iz

            ret = int(ix, c_int64_t) &
                  + int(lx, c_int64_t)*(int(iy - 1, c_int64_t) &
                  + int(ly, c_int64_t)*int(iz - 1, c_int64_t))
        end function

        subroutine decode_cell_index(cell_index, ix, iy, iz)
            integer(c_int64_t), intent(in) :: cell_index
            integer, intent(out) :: ix, iy, iz

            integer(c_int64_t) :: linear0
            integer(c_int64_t) :: lx64, ly64

            lx64 = int(lx, c_int64_t)
            ly64 = int(ly, c_int64_t)
            linear0 = cell_index - 1_c_int64_t

            ix = int(mod(linear0, lx64)) + 1
            iy = int(mod(linear0/lx64, ly64)) + 1
            iz = int(linear0/(lx64*ly64)) + 1
        end subroutine

        integer function cache_hash_slot(cell_index) result(ret)
            integer(c_int64_t), intent(in) :: cell_index

            integer(c_int64_t) :: capacity
            integer(c_int64_t) :: hash_value

            capacity = int(cache_hash_capacity, c_int64_t)
            hash_value = mod(cell_index - 1_c_int64_t, capacity)
            hash_value = mod(hash_value*1103515245_c_int64_t + 12345_c_int64_t, capacity)
            ret = 1 + int(hash_value)
        end function

        integer function flush_lock_index(cell_index) result(ret)
            integer(c_int64_t), intent(in) :: cell_index

            ret = 1 + int(mod(cell_index - 1_c_int64_t, int(FLUSH_LOCK_COUNT, c_int64_t)))
        end function

        integer function worker_index() result(ret)
            ret = 1
!$          ret = omp_get_thread_num() + 1
        end function

        subroutine finalize_cells
            integer :: ix, iy, iz
            real(c_double) :: center, half_width
            real(c_double) :: mean(3)
            real(c_double) :: min_count
            integer :: i, j

            min_count = max(1d0, dble(minimum_count))

            do iz = 1, lz
                do iy = 1, ly
                    do ix = 1, lx
                        if (return_count(ix, iy, iz) <= 0) then
                            return_vx_min(ix, iy, iz) = nan_value
                            return_vx_max(ix, iy, iz) = nan_value
                            return_vy_min(ix, iy, iz) = nan_value
                            return_vy_max(ix, iy, iz) = nan_value
                            return_vz_min(ix, iy, iz) = nan_value
                            return_vz_max(ix, iy, iz) = nan_value
                            if (collect_moments_enabled) then
                                return_mean_v(:, ix, iy, iz) = nan_value
                                return_cov_v(:, :, ix, iy, iz) = nan_value
                            end if
                            return_status(ix, iy, iz) = STATUS_FALLBACK
                            return_confidence(ix, iy, iz) = 0d0
                            cycle
                        end if

                        if (collect_moments_enabled) then
                            mean = return_mean_v(:, ix, iy, iz)/return_weight_sum(ix, iy, iz)
                            return_mean_v(:, ix, iy, iz) = mean
                            do j = 1, 3
                                do i = 1, 3
                                    return_cov_v(i, j, ix, iy, iz) = &
                                        return_cov_v(i, j, ix, iy, iz)/return_weight_sum(ix, iy, iz) &
                                        - mean(i)*mean(j)
                                end do
                            end do
                        end if

                        call expand_range(return_vx_min(ix, iy, iz), return_vx_max(ix, iy, iz))
                        call expand_range(return_vy_min(ix, iy, iz), return_vy_max(ix, iy, iz))
                        call expand_range(return_vz_min(ix, iy, iz), return_vz_max(ix, iy, iz))

                        if (return_count(ix, iy, iz) < minimum_count) then
                            return_status(ix, iy, iz) = STATUS_LOW_COUNT
                        else
                            return_status(ix, iy, iz) = STATUS_OK
                        end if
                        return_confidence(ix, iy, iz) = min(1d0, dble(return_count(ix, iy, iz))/min_count)
                    end do
                end do
            end do
        end subroutine

        subroutine expand_range(vmin, vmax)
            real(c_double), intent(inout) :: vmin
            real(c_double), intent(inout) :: vmax

            real(c_double) :: center
            real(c_double) :: half_width

            center = 0.5d0*(vmin + vmax)
            half_width = 0.5d0*(vmax - vmin)*max(1d0, safety_factor)
            half_width = max(half_width, 10d0*tiny(1d0))
            vmin = center - half_width
            vmax = center + half_width
        end subroutine

        subroutine build_velocity_support(vmean, vthermal, normal_axis, normal_sign, velocities, nvel)
            real(c_double), intent(in) :: vmean(3)
            real(c_double), intent(in) :: vthermal(3)
            integer, intent(in) :: normal_axis
            integer, intent(in) :: normal_sign
            real(c_double), intent(out) :: velocities(3, 32)
            integer, intent(out) :: nvel

            real(c_double) :: candidate(3)
            real(c_double) :: diagonal_scale
            integer :: axis
            integer :: sx, sy, sz

            nvel = 0
            call append_velocity_if_inward(vmean, normal_axis, normal_sign, velocities, nvel)

            do axis = 1, 3
                candidate = vmean
                candidate(axis) = candidate(axis) + coverage_sigma*abs(vthermal(axis))
                call append_velocity_if_inward(candidate, normal_axis, normal_sign, velocities, nvel)

                candidate = vmean
                candidate(axis) = candidate(axis) - coverage_sigma*abs(vthermal(axis))
                call append_velocity_if_inward(candidate, normal_axis, normal_sign, velocities, nvel)
            end do

            diagonal_scale = coverage_sigma/sqrt(3d0)
            do sz = -1, 1, 2
                do sy = -1, 1, 2
                    do sx = -1, 1, 2
                        candidate = vmean + diagonal_scale*abs(vthermal)*dble([sx, sy, sz])
                        call append_velocity_if_inward(candidate, normal_axis, normal_sign, velocities, nvel)
                    end do
                end do
            end do

            candidate = vmean
            candidate(normal_axis) = dble(normal_sign)*max(10d0*tiny(1d0), &
                                                           0.05d0*abs(vthermal(normal_axis)))
            call append_velocity(candidate, velocities, nvel)

            candidate = vmean
            candidate(normal_axis) = dble(normal_sign)* &
                                     (max(0d0, dble(normal_sign)*vmean(normal_axis)) &
                                      + coverage_sigma*abs(vthermal(normal_axis)))
            call append_velocity(candidate, velocities, nvel)

        end subroutine

        subroutine append_velocity_if_inward(candidate, normal_axis, normal_sign, velocities, nvel)
            real(c_double), intent(in) :: candidate(3)
            integer, intent(in) :: normal_axis
            integer, intent(in) :: normal_sign
            real(c_double), intent(inout) :: velocities(:, :)
            integer, intent(inout) :: nvel

            if (dble(normal_sign)*candidate(normal_axis) <= 0d0) then
                return
            end if
            call append_velocity(candidate, velocities, nvel)
        end subroutine

        subroutine append_velocity(candidate, velocities, nvel)
            real(c_double), intent(in) :: candidate(3)
            real(c_double), intent(inout) :: velocities(:, :)
            integer, intent(inout) :: nvel

            if (nvel >= size(velocities, 2)) then
                return
            end if
            nvel = nvel + 1
            velocities(:, nvel) = candidate
        end subroutine

        subroutine tangent_axes(normal_axis, tangent1, tangent2)
            integer, intent(in) :: normal_axis
            integer, intent(out) :: tangent1
            integer, intent(out) :: tangent2

            if (normal_axis == 1) then
                tangent1 = 2
                tangent2 = 3
            else if (normal_axis == 2) then
                tangent1 = 3
                tangent2 = 1
            else
                tangent1 = 1
                tangent2 = 2
            end if
        end subroutine

        integer function dim_size(axis) result(ret)
            integer, intent(in) :: axis

            if (axis == 1) then
                ret = lx
            else if (axis == 2) then
                ret = ly
            else
                ret = lz
            end if
        end function

        integer function first_cell(value, n) result(ret)
            real(c_double), intent(in) :: value
            integer, intent(in) :: n

            ret = min(max(floor(value), 0), n - 1)
        end function

        integer function last_cell(value, n) result(ret)
            real(c_double), intent(in) :: value
            integer, intent(in) :: n

            ret = min(max(ceiling(value) - 1, 0), n - 1)
        end function

        real(c_double) function source_normal_position(normal_axis, normal_sign, patch_min, patch_max) result(ret)
            integer, intent(in) :: normal_axis
            integer, intent(in) :: normal_sign
            real(c_double), intent(in) :: patch_min(3)
            real(c_double), intent(in) :: patch_max(3)

            real(c_double) :: eps

            eps = 1d-9
            if (normal_sign > 0) then
                ret = patch_min(normal_axis) + eps
            else
                ret = patch_max(normal_axis) - eps
            end if

            ret = max(eps, min(dble(dim_size(normal_axis)) - eps, ret))
        end function

    end subroutine

end module

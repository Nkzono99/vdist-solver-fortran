module m_emses_autorange
    !! Deterministic source-envelope tracing for per-cell velocity ranges.

    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value

!$  use omp_lib, only: omp_get_thread_num, omp_set_num_threads

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
        return_vx_min, &
        return_vx_max, &
        return_vy_min, &
        return_vy_max, &
        return_vz_min, &
        return_vz_max, &
        return_count, &
        return_weight_sum, &
        return_mean_v, &
        return_cov_v, &
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
        real(c_double), intent(out) :: return_vx_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vx_max(lx, ly, lz)
        real(c_double), intent(out) :: return_vy_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vy_max(lx, ly, lz)
        real(c_double), intent(out) :: return_vz_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vz_max(lx, ly, lz)
        integer(c_int), intent(out) :: return_count(lx, ly, lz)
        real(c_double), intent(out) :: return_weight_sum(lx, ly, lz)
        real(c_double), intent(out) :: return_mean_v(3, lx, ly, lz)
        real(c_double), intent(out) :: return_cov_v(3, 3, lx, ly, lz)
        integer(c_int), intent(out) :: return_status(lx, ly, lz)
        real(c_double), intent(out) :: return_confidence(lx, ly, lz)
        integer(c_int), value, optional, intent(in) :: n_threads

        real(c_double) :: nan_value
        integer :: active_n_threads
        logical :: collect_moments_enabled
        logical :: show_progress_enabled
        type(bar_object) :: bar
        integer :: progress_total_tasks
        integer :: progress_completed_tasks
        integer :: progress_last_percent
        real(c_double), allocatable :: local_vx_min(:, :, :, :)
        real(c_double), allocatable :: local_vx_max(:, :, :, :)
        real(c_double), allocatable :: local_vy_min(:, :, :, :)
        real(c_double), allocatable :: local_vy_max(:, :, :, :)
        real(c_double), allocatable :: local_vz_min(:, :, :, :)
        real(c_double), allocatable :: local_vz_max(:, :, :, :)
        integer(c_int), allocatable :: local_count(:, :, :, :)
        real(c_double), allocatable :: local_weight_sum(:, :, :, :)
        real(c_double), allocatable :: local_mean_v(:, :, :, :, :)
        real(c_double), allocatable :: local_cov_v(:, :, :, :, :, :)

        nan_value = ieee_value(0d0, ieee_quiet_nan)
        active_n_threads = 1
        if (present(n_threads)) then
            active_n_threads = max(1, int(n_threads))
        end if
        collect_moments_enabled = collect_moments == 1
        show_progress_enabled = .false.
        if (present(show_progress)) then
            show_progress_enabled = show_progress == 1
        end if
        progress_total_tasks = 0
        progress_completed_tasks = 0
        progress_last_percent = -1

!$      call omp_set_num_threads(active_n_threads)

        call initialize_thread_accumulators
        return_status = STATUS_FALLBACK
        return_confidence = 0d0

        call initialize_progress_bar
        call trace_external_boundary_sources
        call trace_emission_surface_sources
        call finish_progress_bar
        call merge_thread_accumulators
        call finalize_cells

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

        subroutine initialize_thread_accumulators
            allocate (local_vx_min(lx, ly, lz, active_n_threads))
            allocate (local_vx_max(lx, ly, lz, active_n_threads))
            allocate (local_vy_min(lx, ly, lz, active_n_threads))
            allocate (local_vy_max(lx, ly, lz, active_n_threads))
            allocate (local_vz_min(lx, ly, lz, active_n_threads))
            allocate (local_vz_max(lx, ly, lz, active_n_threads))
            allocate (local_count(lx, ly, lz, active_n_threads))
            allocate (local_weight_sum(lx, ly, lz, active_n_threads))

            local_vx_min = huge(1d0)
            local_vy_min = huge(1d0)
            local_vz_min = huge(1d0)
            local_vx_max = -huge(1d0)
            local_vy_max = -huge(1d0)
            local_vz_max = -huge(1d0)
            local_count = 0
            local_weight_sum = 0d0

            if (collect_moments_enabled) then
                allocate (local_mean_v(3, lx, ly, lz, active_n_threads))
                allocate (local_cov_v(3, 3, lx, ly, lz, active_n_threads))
                local_mean_v = 0d0
                local_cov_v = 0d0
            end if
        end subroutine

        subroutine merge_thread_accumulators
            integer :: ithread

            return_vx_min = huge(1d0)
            return_vy_min = huge(1d0)
            return_vz_min = huge(1d0)
            return_vx_max = -huge(1d0)
            return_vy_max = -huge(1d0)
            return_vz_max = -huge(1d0)
            return_count = 0
            return_weight_sum = 0d0
            return_mean_v = 0d0
            return_cov_v = 0d0

            do ithread = 1, active_n_threads
                return_vx_min = min(return_vx_min, local_vx_min(:, :, :, ithread))
                return_vx_max = max(return_vx_max, local_vx_max(:, :, :, ithread))
                return_vy_min = min(return_vy_min, local_vy_min(:, :, :, ithread))
                return_vy_max = max(return_vy_max, local_vy_max(:, :, :, ithread))
                return_vz_min = min(return_vz_min, local_vz_min(:, :, :, ithread))
                return_vz_max = max(return_vz_max, local_vz_max(:, :, :, ithread))
                return_count = return_count + local_count(:, :, :, ithread)
                return_weight_sum = return_weight_sum + local_weight_sum(:, :, :, ithread)

                if (collect_moments_enabled) then
                    return_mean_v = return_mean_v + local_mean_v(:, :, :, :, ithread)
                    return_cov_v = return_cov_v + local_cov_v(:, :, :, :, :, ithread)
                end if
            end do
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
            integer :: i, j

            if (position(1) < 0d0 .or. position(1) > dble(lx)) return
            if (position(2) < 0d0 .or. position(2) > dble(ly)) return
            if (position(3) < 0d0 .or. position(3) > dble(lz)) return

            ix = min(max(int(position(1)) + 1, 1), lx)
            iy = min(max(int(position(2)) + 1, 1), ly)
            iz = min(max(int(position(3)) + 1, 1), lz)

            ithread = worker_index()

            local_vx_min(ix, iy, iz, ithread) = min(local_vx_min(ix, iy, iz, ithread), velocity(1))
            local_vx_max(ix, iy, iz, ithread) = max(local_vx_max(ix, iy, iz, ithread), velocity(1))
            local_vy_min(ix, iy, iz, ithread) = min(local_vy_min(ix, iy, iz, ithread), velocity(2))
            local_vy_max(ix, iy, iz, ithread) = max(local_vy_max(ix, iy, iz, ithread), velocity(2))
            local_vz_min(ix, iy, iz, ithread) = min(local_vz_min(ix, iy, iz, ithread), velocity(3))
            local_vz_max(ix, iy, iz, ithread) = max(local_vz_max(ix, iy, iz, ithread), velocity(3))

            local_count(ix, iy, iz, ithread) = local_count(ix, iy, iz, ithread) + 1
            local_weight_sum(ix, iy, iz, ithread) = local_weight_sum(ix, iy, iz, ithread) + weight
            if (collect_moments_enabled) then
                do i = 1, 3
                    local_mean_v(i, ix, iy, iz, ithread) = local_mean_v(i, ix, iy, iz, ithread) &
                                                            + weight*velocity(i)
                end do
                do j = 1, 3
                    do i = 1, 3
                        local_cov_v(i, j, ix, iy, iz, ithread) = local_cov_v(i, j, ix, iy, iz, ithread) &
                                                                  + weight*velocity(i)*velocity(j)
                    end do
                end do
            end if
        end subroutine

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
                            return_mean_v(:, ix, iy, iz) = nan_value
                            return_cov_v(:, :, ix, iy, iz) = nan_value
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
                        else
                            return_mean_v(:, ix, iy, iz) = nan_value
                            return_cov_v(:, :, ix, iy, iz) = nan_value
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

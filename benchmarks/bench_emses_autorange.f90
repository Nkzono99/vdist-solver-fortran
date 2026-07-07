program bench_emses_autorange
    use finbound, only: t_BoundaryList, new_BoundaryList

    use m_allcom, only: max_nepl, max_nspec, nemd, nepl, nflag_emit, path, &
                        peth, phixy, phiz, qm, spa, spe, speth, xmine, xmaxe, &
                        ymine, ymaxe, zmine, zmaxe
    use m_emses_autorange, only: estimate_velocity_range_map_impl
    use m_vdsolverf_core

    implicit none

    integer :: n_threads
    integer :: lx, ly, lz
    integer :: max_step
    integer :: reps
    integer :: irep
    integer :: rate
    integer :: start_count
    integer :: end_count
    real(8) :: elapsed

    type(t_ESSimulator) :: simulator
    type(t_VectorFieldGrid) :: eb
    type(t_BoundaryList) :: boundaries
    type(tp_Probability), allocatable :: probability_functions(:)
    real(8), allocatable :: values(:, :, :, :)
    real(8), allocatable :: vx_min(:, :, :), vx_max(:, :, :)
    real(8), allocatable :: vy_min(:, :, :), vy_max(:, :, :)
    real(8), allocatable :: vz_min(:, :, :), vz_max(:, :, :)
    real(8), allocatable :: weight_sum(:, :, :)
    real(8), allocatable :: mean_v(:, :, :, :)
    real(8), allocatable :: cov_v(:, :, :, :, :)
    real(8), allocatable :: confidence(:, :, :)
    integer, allocatable :: count(:, :, :)
    integer, allocatable :: status(:, :, :)
    integer :: boundary_conditions(3)

    n_threads = read_arg_int(1, 1)
    lx = read_arg_int(2, 32)
    ly = read_arg_int(3, 32)
    lz = read_arg_int(4, 32)
    max_step = read_arg_int(5, 96)
    reps = read_arg_int(6, 5)

    call reset_emission_globals(lx, ly, lz)

    allocate (values(6, 0:lx, 0:ly, 0:lz))
    allocate (vx_min(lx, ly, lz), vx_max(lx, ly, lz))
    allocate (vy_min(lx, ly, lz), vy_max(lx, ly, lz))
    allocate (vz_min(lx, ly, lz), vz_max(lx, ly, lz))
    allocate (weight_sum(lx, ly, lz))
    allocate (mean_v(3, lx, ly, lz))
    allocate (cov_v(3, 3, lx, ly, lz))
    allocate (confidence(lx, ly, lz))
    allocate (count(lx, ly, lz))
    allocate (status(lx, ly, lz))

    values = 0d0
    eb = new_VectorFieldGrid(6, lx, ly, lz, values)
    boundaries = new_BoundaryList()
    allocate (probability_functions(1))
    allocate (probability_functions(1)%ref, source=new_ZeroProbability())

    boundary_conditions = [2, 2, 2]
    simulator = new_ESSimulator(lx, ly, lz, boundary_conditions, eb, boundaries, probability_functions)

    call run_once

    call system_clock(start_count, rate)
    do irep = 1, reps
        call run_once
    end do
    call system_clock(end_count)

    elapsed = dble(end_count - start_count)/dble(rate)
    print '(a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0,1x,a,i0)', &
        "threads=", n_threads, "lx=", lx, "ly=", ly, "lz=", lz, "max_step=", max_step, "reps=", reps
    print '(a,es14.6,1x,a,es14.6,1x,a,i0)', &
        "elapsed_s=", elapsed, "per_run_s=", elapsed/dble(reps), "total_count=", sum(count)

contains

    subroutine run_once
        call estimate_velocity_range_map_impl( &
            simulator=simulator, &
            ispec=1, &
            lx=lx, ly=ly, lz=lz, &
            trace_dt=0.5d0, &
            coverage_sigma=1d0, &
            safety_factor=1d0, &
            max_step=max_step, &
            use_adaptive_dt=0, &
            source_samples_per_cell=1, &
            velocity_sample_mode=0, &
            minimum_count=1, &
            collect_moments=0, &
            n_threads=n_threads, &
            return_vx_min=vx_min, &
            return_vx_max=vx_max, &
            return_vy_min=vy_min, &
            return_vy_max=vy_max, &
            return_vz_min=vz_min, &
            return_vz_max=vz_max, &
            return_count=count, &
            return_weight_sum=weight_sum, &
            return_mean_v=mean_v, &
            return_cov_v=cov_v, &
            return_status=status, &
            return_confidence=confidence)
    end subroutine

    subroutine reset_emission_globals(nx, ny, nz)
        integer, intent(in) :: nx, ny, nz

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
        ymaxe(1) = dble(ny)
        zmine(1) = 0d0
        zmaxe(1) = dble(nz)
    end subroutine

    integer function read_arg_int(index, default_value) result(ret)
        integer, intent(in) :: index
        integer, intent(in) :: default_value

        character(len=64) :: buffer
        integer :: status

        ret = default_value
        if (command_argument_count() < index) return

        call get_command_argument(index, buffer, status=status)
        if (status /= 0) return

        read (buffer, *, iostat=status) ret
        if (status /= 0) ret = default_value
    end function

end program

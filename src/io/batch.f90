module m_vdsolverf_batch
    !! Serial batch engine shared by the CLI and native regressions. No MPI here.
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: iso_c_binding, only: c_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!$  use omp_lib, only: omp_get_max_threads
    use m_batch_support
    use m_npy, only: read_npy, write_npy
    use m_field_bundle, only: load_field_bundle
    use m_simulator, only: t_ESSimulator
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_emses_dust_currents, only: evaluate_dust_currents_impl
    implicit none
    private
    public :: t_batch_config, t_batch, read_batch_config, open_batch, close_batch
    public :: compute_chunk, chunk_complete, write_summary, chunk_prefix
    public :: prepare_work
    type t_batch_config
        character(32) :: mode = 'dust_currents'
        character(:), allocatable :: bundle_dir, output_dir, directory
        integer :: ispec = -1, max_step = 30000, max_probability_types = 100
        integer :: scout_bins(3) = [9, 9, 9], max_depth = 6, max_samples_per_cell = 50000
        integer :: max_leaves_per_cell = 4096, max_expansions = 3, chunk_size = 16, n_threads = 1
        real(real64) :: dt = -1d0, refine_threshold_rel = 1d-4, edge_threshold_rel = 1d-4, expand_factor = 1.5d0
        logical :: use_adaptive_dt = .false., use_electric_field = .true., use_magnetic_field = .true.
        real(real64) :: coverage_sigma = 4d0, safety_factor = 1.25d0
        integer :: source_samples_per_cell = 1, minimum_count = 4, accumulator_cache_size = 20000
    end type
    type t_batch
        type(t_batch_config) :: config
        type(t_ESSimulator) :: simulator
        real(real64), allocatable :: targets(:, :), psi(:)
        integer :: n_targets = 0, npsi = 0, n_chunks = 0
    end type
contains
    subroutine read_batch_config(path, config)
        character(*), intent(in) :: path
        type(t_batch_config), intent(out) :: config
        character(4096) :: bundle_dir, output_dir
        character(32) :: mode
        character(1024) :: message
        integer :: ispec, max_step, max_probability_types, scout_bins(3), max_depth, max_samples_per_cell
        integer :: max_leaves_per_cell, max_expansions, chunk_size, n_threads, unit, ios
        real(real64) :: dt, refine_threshold_rel, edge_threshold_rel, expand_factor
        real(real64) :: coverage_sigma, safety_factor
        integer :: source_samples_per_cell, minimum_count, accumulator_cache_size
        logical :: use_adaptive_dt, use_electric_field, use_magnetic_field
        namelist /vdsolverf_batch/ bundle_dir, ispec, dt, max_step, use_adaptive_dt, max_probability_types, &
            scout_bins, max_depth, max_samples_per_cell, max_leaves_per_cell, refine_threshold_rel, &
            edge_threshold_rel, expand_factor, max_expansions, chunk_size, output_dir, n_threads, &
            use_electric_field, use_magnetic_field, mode, coverage_sigma, safety_factor, &
            source_samples_per_cell, minimum_count, accumulator_cache_size
        mode = config%mode
        coverage_sigma = config%coverage_sigma
        safety_factor = config%safety_factor
        source_samples_per_cell = config%source_samples_per_cell
        minimum_count = config%minimum_count
        accumulator_cache_size = config%accumulator_cache_size
        bundle_dir = ''
        output_dir = 'out'
        ispec = config%ispec
        dt = config%dt
        max_step = config%max_step
        max_probability_types = config%max_probability_types
        scout_bins = config%scout_bins
        max_depth = config%max_depth
        max_samples_per_cell = config%max_samples_per_cell
        max_leaves_per_cell = config%max_leaves_per_cell
        refine_threshold_rel = config%refine_threshold_rel
        edge_threshold_rel = config%edge_threshold_rel
        expand_factor = config%expand_factor
        max_expansions = config%max_expansions
        chunk_size = config%chunk_size
        n_threads = 1
!$      n_threads = omp_get_max_threads()
        use_adaptive_dt = config%use_adaptive_dt
        use_electric_field = config%use_electric_field
        use_magnetic_field = config%use_magnetic_field
        open(newunit=unit, file=path, status='old', iostat=ios, iomsg=message)
        if (ios /= 0) call batch_fail('cannot open batch namelist '//path//': '//trim(message))
        read(unit, nml=vdsolverf_batch, iostat=ios, iomsg=message)
        close(unit)
        if (ios /= 0) call batch_fail('invalid &vdsolverf_batch: '//trim(message))
        if (mode /= 'dust_currents' .and. mode /= 'velocity_range') &
            call batch_fail('mode must be dust_currents or velocity_range')
        if (.not. all(ieee_is_finite([coverage_sigma, safety_factor]))) &
            call batch_fail('coverage_sigma/safety_factor must be finite')
        if (coverage_sigma <= 0d0 .or. safety_factor < 1d0) &
            call batch_fail('coverage_sigma must be positive and safety_factor >= 1')
        if (min(source_samples_per_cell, minimum_count, accumulator_cache_size) < 1) &
            call batch_fail('source_samples_per_cell/minimum_count/accumulator_cache_size must be positive')
        if (len_trim(bundle_dir) == 0 .or. len_trim(output_dir) == 0) call batch_fail('bundle_dir/output_dir required')
        if (ispec < 0 .or. ispec == huge(ispec)) call batch_fail('ispec must be a nonnegative 0-based species index')
        if (.not. ieee_is_finite(dt) .or. dt <= 0d0 .or. max_step < 1) call batch_fail('dt/max_step must be positive')
        if (n_threads < 1 .or. max_probability_types < 1 .or. chunk_size < 1) &
            call batch_fail('n_threads/max_probability_types/chunk_size must be positive')
        if (any(scout_bins < 2)) call batch_fail('scout_bins must contain three integers >= 2')
        if (product(real(scout_bins, real64)) > real(huge(1), real64)) call batch_fail('scout_bins product too large')
        if (max_depth < 0 .or. max_expansions < 0) call batch_fail('max_depth/max_expansions must be nonnegative')
        if (max_samples_per_cell < 1 .or. max_leaves_per_cell < 1) call batch_fail('octree budgets must be positive')
        if (.not. all(ieee_is_finite([refine_threshold_rel, edge_threshold_rel, expand_factor]))) &
            call batch_fail('octree thresholds/expand_factor must be finite')
        if (min(refine_threshold_rel, edge_threshold_rel) < 0d0 .or. expand_factor < 1d0) &
            call batch_fail('thresholds must be >= 0 and expand_factor >= 1')
        config%directory = parent_path(path)
        config%mode = mode
        config%coverage_sigma = coverage_sigma
        config%safety_factor = safety_factor
        config%source_samples_per_cell = source_samples_per_cell
        config%minimum_count = minimum_count
        config%accumulator_cache_size = accumulator_cache_size
        config%bundle_dir = path_join(config%directory, trim(bundle_dir))
        config%output_dir = path_join(config%directory, trim(output_dir))
        config%ispec = ispec
        config%dt = dt
        config%max_step = max_step
        config%use_adaptive_dt = use_adaptive_dt
        config%max_probability_types = max_probability_types
        config%scout_bins = scout_bins
        config%max_depth = max_depth
        config%max_samples_per_cell = max_samples_per_cell
        config%max_leaves_per_cell = max_leaves_per_cell
        config%refine_threshold_rel = refine_threshold_rel
        config%edge_threshold_rel = edge_threshold_rel
        config%expand_factor = expand_factor
        config%max_expansions = max_expansions
        config%chunk_size = chunk_size
        config%n_threads = n_threads
        config%use_electric_field = use_electric_field
        config%use_magnetic_field = use_magnetic_field
    end subroutine

    subroutine open_batch(path, batch)
        character(*), intent(in) :: path
        type(t_batch), intent(out) :: batch
        real(real64), allocatable :: flat(:), fields(:, :, :, :)
        integer(int64), allocatable :: dims(:)
        integer :: dimensions(3), i
        character(:), allocatable :: inppath
        character(c_char), allocatable :: cpath(:)
        call read_batch_config(path, batch%config)
        call read_npy(path_join(batch%config%directory, 'targets.npy'), flat, dims)
        if (size(dims) /= 2) call batch_fail('targets.npy must have shape (n,9)')
        if (dims(1) < 1 .or. dims(1) > huge(1) .or. dims(2) /= 9) call batch_fail('targets.npy must have shape (n,9)')
        batch%n_targets = int(dims(1))
        batch%targets = reshape(flat, [9, batch%n_targets])
        if (.not. all(ieee_is_finite(batch%targets(:3, :)))) call batch_fail('target positions must be finite')
        ! Invalid velocity boxes deliberately reach the core and return status 5.
        call read_npy(path_join(batch%config%directory, 'psi.npy'), batch%psi, dims)
        if (size(dims) /= 1) call batch_fail('psi.npy must have shape (npsi,)')
        if (dims(1) < 1 .or. dims(1) > huge(1)) call batch_fail('psi.npy must be nonempty')
        if (.not. all(ieee_is_finite(batch%psi))) call batch_fail('psi.npy must contain finite values')
        batch%npsi = int(dims(1))
        batch%n_chunks = (batch%n_targets-1)/batch%config%chunk_size+1
        call load_field_bundle(batch%config%bundle_dir, batch%config%ispec, batch%config%use_electric_field, &
            batch%config%use_magnetic_field, fields, dimensions)
        inppath = path_join(batch%config%bundle_dir, 'plasma-vdsolverf.inp')
        allocate(cpath(len(inppath)))
        do i = 1, len(inppath)
            cpath(i) = inppath(i:i)
        end do
        batch%simulator = create_simulator(cpath, len(inppath), dimensions(1), dimensions(2), dimensions(3), &
            fields, batch%config%ispec+1, batch%config%max_probability_types)
    end subroutine

    subroutine close_batch(batch)
        type(t_batch), intent(inout) :: batch
        call destroy_simulator(batch%simulator)
    end subroutine

    function chunk_prefix(batch, chunk) result(prefix)
        type(t_batch), intent(in) :: batch
        integer, intent(in) :: chunk
        character(:), allocatable :: prefix
        character(32) :: buffer
        if (chunk < 0 .or. chunk >= batch%n_chunks) call batch_fail('chunk index out of range')
        write(buffer, '(i0.9)') chunk*batch%config%chunk_size
        prefix = path_join(batch%config%output_dir, 'chunk_'//trim(buffer))
    end function

    function chunk_complete(batch, chunk) result(complete)
        type(t_batch), intent(in) :: batch
        integer, intent(in) :: chunk
        logical :: complete, exists
        character(:), allocatable :: prefix
        character(24), parameter :: suffixes(6) = [character(24) :: '_flux.npy', '_flux_error.npy', &
            '_status.npy', '_sample_count.npy', '_terminal_box_count.npy', '.done']
        integer :: i
        prefix = chunk_prefix(batch, chunk)
        complete = .true.
        do i = 1, size(suffixes)
            inquire(file=prefix//trim(suffixes(i)), exist=exists)
            complete = complete .and. exists
        end do
    end function

    subroutine prepare_work(batch, resume, work, computed, skipped)
        type(t_batch), intent(in) :: batch
        logical, intent(in) :: resume
        integer, allocatable, intent(out) :: work(:)
        integer, intent(out) :: computed, skipped
        integer :: chunk
        call make_directory(batch%config%output_dir)
        call remove_file(path_join(batch%config%output_dir, 'summary.json'))
        allocate(work(batch%n_chunks))
        computed = 0
        skipped = 0
        do chunk = 0, batch%n_chunks-1
            ! Invalidate the old run before workers start a forced recomputation.
            if (.not. resume) call remove_file(chunk_prefix(batch, chunk)//'.done')
            if (resume) then
                if (chunk_complete(batch, chunk)) then
                    skipped = skipped+1
                    cycle
                end if
            end if
            computed = computed+1
            work(computed) = chunk
        end do
    end subroutine

    subroutine compute_chunk(batch, chunk)
        type(t_batch), intent(in) :: batch
        integer, intent(in) :: chunk
        real(real64), allocatable :: flux(:, :), error(:, :)
        integer(int32), allocatable :: status(:), samples(:), terminals(:)
        character(:), allocatable :: prefix, marker
        integer :: first, last, count, unit, ios
        first = chunk*batch%config%chunk_size + 1
        last = min(int(first, int64)+batch%config%chunk_size-1, int(batch%n_targets, int64))
        count = last-first+1
        prefix = chunk_prefix(batch, chunk)
        marker = prefix//'.done'
        call remove_file(marker)
        allocate(flux(batch%npsi, count), error(batch%npsi, count), status(count), samples(count), terminals(count))
        associate(c => batch%config)
            call evaluate_dust_currents_impl(batch%simulator, c%ispec+1, count, batch%targets(:3, first:last), &
                batch%targets(4:, first:last), batch%npsi, batch%psi, c%dt, c%max_step, merge(1, 0, c%use_adaptive_dt), &
                c%scout_bins, c%max_depth, c%max_samples_per_cell, c%max_leaves_per_cell, c%refine_threshold_rel, &
                c%edge_threshold_rel, c%expand_factor, c%max_expansions, flux, error, status, samples, terminals, c%n_threads)
        end associate
        call write_npy(prefix//'_flux.npy.tmp', reshape(flux, [size(flux)]), [int(count, int64), int(batch%npsi, int64)])
        call atomic_replace(prefix//'_flux.npy.tmp', prefix//'_flux.npy')
        call write_npy(prefix//'_flux_error.npy.tmp', reshape(error, [size(error)]), &
            [int(count, int64), int(batch%npsi, int64)])
        call atomic_replace(prefix//'_flux_error.npy.tmp', prefix//'_flux_error.npy')
        call write_npy(prefix//'_status.npy.tmp', status, [int(count, int64)])
        call atomic_replace(prefix//'_status.npy.tmp', prefix//'_status.npy')
        call write_npy(prefix//'_sample_count.npy.tmp', int(samples, int64), [int(count, int64)])
        call atomic_replace(prefix//'_sample_count.npy.tmp', prefix//'_sample_count.npy')
        call write_npy(prefix//'_terminal_box_count.npy.tmp', int(terminals, int64), [int(count, int64)])
        call atomic_replace(prefix//'_terminal_box_count.npy.tmp', prefix//'_terminal_box_count.npy')
        open(newunit=unit, file=marker//'.tmp', status='replace', iostat=ios)
        if (ios /= 0) call batch_fail('cannot create completion marker '//marker)
        close(unit, iostat=ios)
        if (ios /= 0) call batch_fail('cannot close completion marker '//marker)
        call atomic_replace(marker//'.tmp', marker)
    end subroutine

    subroutine write_summary(batch, seconds, ranks, computed, skipped)
        type(t_batch), intent(in) :: batch
        real(real64), intent(in) :: seconds
        integer, intent(in) :: ranks, computed, skipped
        integer(int64) :: counts(0:5)
        integer(int32), allocatable :: status(:)
        integer(int64), allocatable :: dims(:)
        integer :: chunk, code, unit, ios, count
        character(:), allocatable :: path
        counts = 0
        do chunk = 0, batch%n_chunks-1
            if (.not. chunk_complete(batch, chunk)) call batch_fail('missing completed chunk at summary time')
            call read_npy(chunk_prefix(batch, chunk)//'_status.npy', status, dims)
            count = min(batch%config%chunk_size, batch%n_targets-chunk*batch%config%chunk_size)
            if (size(dims) /= 1) call batch_fail('invalid completed chunk status shape; use --no-resume')
            if (dims(1) /= count) call batch_fail('completed chunk size changed; use --no-resume')
            if (any(status < 0) .or. any(status > 5)) call batch_fail('invalid completed chunk status; use --no-resume')
            do code = 0, 5
                counts(code) = counts(code) + int(count_status(status, code), int64)
            end do
        end do
        path = path_join(batch%config%output_dir, 'summary.json')
        open(newunit=unit, file=path//'.tmp', status='replace', iostat=ios)
        if (ios /= 0) call batch_fail('cannot create summary '//path)
        write(unit, '(a)') '{'
        write(unit, '(a)') '  "schema": "vdsolverf.batch_summary",'
        write(unit, '(a,i0,a)') '  "n_targets": ', batch%n_targets, ','
        write(unit, '(a,i0,a)') '  "npsi": ', batch%npsi, ','
        write(unit, '(a,es24.16e3,a)') '  "wall_time_seconds": ', seconds, ','
        write(unit, '(a,i0,a)') '  "ranks": ', ranks, ','
        write(unit, '(a,i0,a)') '  "threads": ', batch%config%n_threads, ','
        write(unit, '(a,i0,a)') '  "chunks_computed": ', computed, ','
        write(unit, '(a,i0,a)') '  "chunks_skipped": ', skipped, ','
        write(unit, '(a)') '  "status_counts": {'
        do code = 0, 5
            if (code < 5) then
                write(unit, '(a,i0,a,i0,a)') '    "', code, '": ', counts(code), ','
            else
                write(unit, '(a,i0,a,i0)') '    "', code, '": ', counts(code)
            end if
        end do
        write(unit, '(a)') '  }'
        write(unit, '(a)', iostat=ios) '}'
        if (ios /= 0) call batch_fail('cannot write summary '//path)
        close(unit, iostat=ios)
        if (ios /= 0) call batch_fail('cannot finish summary '//path)
        call atomic_replace(path//'.tmp', path)
    end subroutine

    integer function count_status(status, code) result(n)
        integer(int32), intent(in) :: status(:)
        integer, intent(in) :: code
        n = count(status == code)
    end function
end module

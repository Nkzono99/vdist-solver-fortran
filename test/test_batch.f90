program test_batch
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use m_vdsolverf_batch
    use m_npy, only: read_npy, write_npy
    use m_field_bundle, only: load_field_bundle
    use m_batch_support
    use m_emses_dust_currents, only: evaluate_dust_currents_impl
    use m_test_helpers, only: assert_true, assert_equal_int
    implicit none
    character(*), parameter :: root = 'test-vdsolverf-batch'
    type(t_batch) :: batch
    real(real64) :: fields(9, 3, 2, 4), targets(9, 5), psi(3), flux(3, 5), errors(3, 5)
    real(real64) :: parallel_flux(3, 5), parallel_errors(3, 5)
    real(real64), allocatable :: flat(:), loaded(:, :, :, :)
    integer(int32) :: status(5), samples(5), terminals(5)
    integer(int32) :: parallel_status(5), parallel_samples(5), parallel_terminals(5)
    integer(int32), allocatable :: stored_status(:)
    integer(int64), allocatable :: dims(:), stored_count(:)
    integer, allocatable :: work(:)
    integer :: i, u, computed, skipped, chunk, first, last, dimensions(3), repeat
    character(:), allocatable :: prefix
    call make_directory(root//'/bundle')
    open(newunit=u, file=root//'/bundle/plasma-vdsolverf.inp', status='replace')
    write(u, '(a)') '&esorem emflag=0 /'
    write(u, '(a)') '&plasma wp=10*0, wc=0, phiz=0, phixy=0 /'
    write(u, '(a)') '&tmgrid nx=2, ny=1, nz=3, dt=0.1 /'
    write(u, '(a)') '&system nspec=1, npbnd=30*2 /'
    write(u, '(a)') '&intp qm=10*1, peth=10*1, path=10*1, vdri=10*0, vdthz=10*0, vdthxy=10*0 /'
    write(u, '(a)') "&ptcond boundary_type='none' /"
    write(u, '(a)') '&emissn nepl=10*0, nflag_emit=10*0 /'
    close(u)
    open(newunit=u, file=root//'/bundle/bundle.json', status='replace')
    write(u, '(a)') '{"schema":"vdsolverf.field_bundle", "version":"0.1.0", "nx":2, "ny":1, "nz":3,'
    write(u, '(a)') '"ispecs":[0], "use_electric_field":true, "use_magnetic_field":true}'
    close(u)
    fields = reshape([(1d-5*real(i, real64), i=1, size(fields))], shape(fields))
    call write_npy(root//'/bundle/ebvalues_ispec0.npy', reshape(fields, [size(fields)]), &
        [4_int64, 2_int64, 3_int64, 9_int64], version=2)
    call load_field_bundle(root//'/bundle', 0, .true., .true., loaded, dimensions)
    call assert_true('C field dimensions reverse exactly', all(shape(loaded) == shape(fields)))
    call assert_true('C field layout exact (noncubic, nonzero)', all(loaded == fields))
    do i = 1, 5
        targets(:, i) = [0.3d0+0.2d0*i, 0.5d0, 1.5d0, -4d0, 4d0, -4d0, 4d0, -4d0, 4d0]
    end do
    targets(5, 5) = -5d0 ! status 5 must survive the disk format and summary
    psi = [-1d0, 0d0, 1d0]
    call write_npy(root//'/targets.npy', reshape(targets, [size(targets)]), [5_int64, 9_int64])
    call write_npy(root//'/psi.npy', psi, [3_int64])
    open(newunit=u, file=root//'/batch.nml', status='replace')
    write(u, '(a)') "&vdsolverf_batch bundle_dir='bundle', ispec=0, dt=0.7, max_step=32,"
    write(u, '(a)') 'scout_bins=3,4,3, max_depth=2, max_samples_per_cell=4096, max_leaves_per_cell=128,'
    write(u, '(a)') 'chunk_size=2, n_threads=4, max_expansions=0 /'
    close(u)
    call open_batch(root//'/batch.nml', batch)
    call evaluate_dust_currents_impl(batch%simulator, 1, 5, targets(:3, :), targets(4:, :), 3, psi, &
        0.7d0, 32, 0, [3, 4, 3], 2, 4096, 128, 1d-4, 1d-4, 1.5d0, 0, flux, errors, status, samples, terminals, 1)
    do repeat = 1, 3
        call evaluate_dust_currents_impl(batch%simulator, 1, 5, targets(:3, :), targets(4:, :), 3, psi, &
            0.7d0, 32, 0, [3, 4, 3], 2, 4096, 128, 1d-4, 1d-4, 1.5d0, 0, parallel_flux, parallel_errors, &
            parallel_status, parallel_samples, parallel_terminals, 4)
        call assert_true('nonuniform field parallel flux deterministic', all(parallel_flux == flux))
        call assert_true('nonuniform field parallel error deterministic', all(parallel_errors == errors))
        call assert_true('parallel diagnostics exact', all(parallel_status == status) .and. &
            all(parallel_samples == samples) .and. all(parallel_terminals == terminals))
    end do
    call prepare_work(batch, .true., work, computed, skipped)
    call assert_equal_int('initial chunks computed', computed, 3)
    do i = 1, computed
        call compute_chunk(batch, work(i))
    end do
    call write_summary(batch, 0d0, 1, computed, skipped)
    do chunk = 0, batch%n_chunks-1
        first = chunk*2+1
        last = min(first+1, 5)
        prefix = chunk_prefix(batch, chunk)
        call read_npy(prefix//'_flux.npy', flat, dims)
        call assert_true('batch flux equals core bit for bit', all(flat == reshape(flux(:, first:last), [size(flat)])))
        call read_npy(prefix//'_flux_error.npy', flat, dims)
        call assert_true('batch error equals core bit for bit', all(flat == reshape(errors(:, first:last), [size(flat)])))
        call read_npy(prefix//'_status.npy', stored_status, dims)
        call assert_true('batch status equals core', all(stored_status == status(first:last)))
        call read_npy(prefix//'_sample_count.npy', stored_count, dims)
        call assert_true('batch sample counts equal core', all(stored_count == samples(first:last)))
        call read_npy(prefix//'_terminal_box_count.npy', stored_count, dims)
        call assert_true('batch terminal counts equal core', all(stored_count == terminals(first:last)))
    end do
    call prepare_work(batch, .true., work, computed, skipped)
    call assert_equal_int('resume computes no completed chunks', computed, 0)
    call assert_equal_int('resume skips completed chunks', skipped, 3)
    prefix = chunk_prefix(batch, 1)
    call remove_file(prefix//'.done')
    call write_npy(prefix//'_flux.npy', [-99d0], [1_int64])
    call prepare_work(batch, .true., work, computed, skipped)
    call assert_equal_int('missing marker recomputes exactly one chunk', computed, 1)
    call assert_equal_int('missing marker keeps completed chunks', skipped, 2)
    call assert_equal_int('correct incomplete chunk selected', work(1), 1)
    call compute_chunk(batch, work(1))
    call read_npy(prefix//'_flux.npy', flat, dims)
    call assert_true('partial chunk replaced by valid data', all(flat == reshape(flux(:, 3:4), [size(flat)])))
    call prepare_work(batch, .false., work, computed, skipped)
    call assert_equal_int('no-resume schedules every chunk', computed, 3)
    call assert_equal_int('no-resume skips none', skipped, 0)
    do chunk = 0, batch%n_chunks-1
        call assert_true('no-resume invalidates old markers before dispatch', .not. chunk_complete(batch, chunk))
    end do
    call close_batch(batch)
    call execute_command_line('rm -rf '//root)
    print *, 'test_batch: all tests passed.'
end program

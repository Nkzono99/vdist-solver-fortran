module m_range_batch
    !! Standalone tier-2 range estimation; uses exactly the ctypes core.
    use, intrinsic :: iso_fortran_env, only: real64, int8, int32, int64
    use, intrinsic :: iso_c_binding, only: c_char, c_null_ptr
    use m_batch_support
    use m_vdsolverf_batch, only: t_batch_config
    use m_field_bundle, only: load_field_bundle
    use m_npy, only: write_npy
    use m_simulator, only: t_ESSimulator
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_emses_autorange, only: estimate_velocity_range_map_impl
    implicit none
    private
    public :: run_range_batch
    character(16), parameter :: real_names(8) = [character(16) :: &
        'vx_min', 'vx_max', 'vy_min', 'vy_max', 'vz_min', 'vz_max', 'weight_sum', 'confidence']
    character(16), parameter :: int_names(2) = [character(16) :: 'count', 'status']
contains
    function input_fingerprint(path, config) result(fingerprint)
        !! CRC64 of all input bytes, including field payload and effective threads.
        !! Detects accidental changes even if file size and timestamp are preserved.
        character(*), intent(in) :: path
        type(t_batch_config), intent(in) :: config
        character(16) :: fingerprint
        integer(int64) :: crc, table(0:255), v
        integer :: i, bit
        do i = 0, 255
            v = int(i, int64)
            do bit = 1, 8
                if (btest(v, 0)) then
                    ! Signed representation of CRC64/XZ polynomial C96C5795D7870F42.
                    ! Decimal avoids GNU 8's BOZ-to-signed conversion range error.
                    v = ieor(ishft(v, -1), -3932672073523589310_int64)
                else
                    v = ishft(v, -1)
                end if
            end do
            table(i) = v
        end do
        crc = not(0_int64)
        call hash_file(path, crc, table)
        call hash_file(path_join(config%bundle_dir, 'bundle.json'), crc, table)
        call hash_file(path_join(config%bundle_dir, 'plasma-vdsolverf.inp'), crc, table)
        call hash_file(path_join(config%bundle_dir, 'ebvalues_ispec'//number(config%ispec)//'.npy'), crc, table)
        v = int(config%n_threads, int64)
        crc = ieor(crc, v)
        write(fingerprint, '(z16.16)') not(crc)
    end function

    subroutine hash_file(path, crc, table)
        character(*), intent(in) :: path
        integer(int64), intent(inout) :: crc
        integer(int64), intent(in) :: table(0:255)
        integer(int8) :: buffer(65536)
        integer(int64) :: length, position
        integer :: unit, ios, n, i, index_byte
        inquire(file=path, size=length, iostat=ios)
        if (ios /= 0 .or. length < 0) call batch_fail('cannot fingerprint '//path)
        open(newunit=unit, file=path, status='old', access='stream', form='unformatted', iostat=ios)
        if (ios /= 0) call batch_fail('cannot fingerprint '//path)
        position = 1
        do while (position <= length)
            n = int(min(int(size(buffer), int64), length-position+1))
            read(unit, pos=position, iostat=ios) buffer(:n)
            if (ios /= 0) call batch_fail('cannot fingerprint bytes of '//path)
            do i = 1, n
                index_byte = int(iand(ieor(crc, int(buffer(i), int64)), 255_int64))
                crc = ieor(ishft(crc, -8), table(index_byte))
            end do
            position = position+n
        end do
        close(unit)
        crc = ieor(crc, length) ! Preserve file boundaries in the combined digest.
    end subroutine

    subroutine run_range_batch(path, config, resume, ranks)
        character(*), intent(in) :: path
        type(t_batch_config), intent(in) :: config
        logical, intent(in) :: resume
        integer, intent(in) :: ranks
        type(t_ESSimulator) :: simulator
        character(:), allocatable :: output, summary, inppath, metadata, old_summary, filename
        character(c_char), allocatable :: cpath(:)
        character(16) :: fingerprint
        real(real64), allocatable :: fields(:, :, :, :), values(:, :, :, :)
        integer(int32), allocatable :: integers(:, :, :, :)
        integer :: dimensions(3), i, unit, ios, code
        integer(int64) :: clock_start, clock_end, clock_rate, dims(3)
        logical :: complete, exists
        call system_clock(clock_start, clock_rate)
        fingerprint = input_fingerprint(path, config)
        output = path_join(config%output_dir, 'velocity_range')
        summary = path_join(output, 'summary.json')
        inquire(file=summary, exist=exists)
        if (resume .and. exists) then
            old_summary = read_text(summary)
            complete = json_string(old_summary, 'input_fingerprint') == fingerprint
            do i = 1, size(real_names)
                inquire(file=path_join(output, trim(real_names(i))//'.npy'), exist=exists)
                complete = complete .and. exists
            end do
            do i = 1, size(int_names)
                inquire(file=path_join(output, trim(int_names(i))//'.npy'), exist=exists)
                complete = complete .and. exists
            end do
            if (complete) then
                print '(a)', 'vdsolverf-batch: velocity_range skipped (inputs unchanged)'
                return
            end if
        end if
        call make_directory(output)
        call remove_file(summary)
        ! Range mode inherits field flags from the bundle; no EMSES files or Python.
        metadata = read_text(path_join(config%bundle_dir, 'bundle.json'))
        call load_field_bundle(config%bundle_dir, config%ispec, &
            json_logical(metadata, 'use_electric_field'), json_logical(metadata, 'use_magnetic_field'), &
            fields, dimensions)
        inppath = path_join(config%bundle_dir, 'plasma-vdsolverf.inp')
        allocate(cpath(len(inppath)))
        do i = 1, len(inppath)
            cpath(i) = inppath(i:i)
        end do
        simulator = create_simulator(cpath, len(inppath), dimensions(1), dimensions(2), dimensions(3), &
            fields, config%ispec+1, config%max_probability_types)
        deallocate(fields)
        allocate(values(dimensions(1), dimensions(2), dimensions(3), 8), &
                 integers(dimensions(1), dimensions(2), dimensions(3), 2))
        associate(c => config)
            call estimate_velocity_range_map_impl(simulator, c%ispec+1, dimensions(1), dimensions(2), dimensions(3), &
                c%dt, c%coverage_sigma, c%safety_factor, c%max_step, merge(1, 0, c%use_adaptive_dt), &
                c%source_samples_per_cell, 0, c%minimum_count, 0, 0, c%accumulator_cache_size, &
                values(:, :, :, 1), values(:, :, :, 2), values(:, :, :, 3), values(:, :, :, 4), &
                values(:, :, :, 5), values(:, :, :, 6), integers(:, :, :, 1), values(:, :, :, 7), &
                c_null_ptr, c_null_ptr, integers(:, :, :, 2), values(:, :, :, 8), c%n_threads)
        end associate
        call destroy_simulator(simulator)
        dims = int(dimensions(3:1:-1), int64)
        do i = 1, size(real_names)
            filename = path_join(output, trim(real_names(i))//'.npy')
            call write_npy(filename//'.tmp', reshape(values(:, :, :, i), [size(values(:, :, :, i))]), dims)
            call atomic_replace(filename//'.tmp', filename)
        end do
        do i = 1, size(int_names)
            filename = path_join(output, trim(int_names(i))//'.npy')
            call write_npy(filename//'.tmp', reshape(integers(:, :, :, i), [size(integers(:, :, :, i))]), dims)
            call atomic_replace(filename//'.tmp', filename)
        end do
        call system_clock(clock_end)
        open(newunit=unit, file=summary//'.tmp', status='replace', iostat=ios)
        if (ios /= 0) call batch_fail('cannot create range summary')
        write(unit, '(a)') '{"schema": "vdsolverf.velocity_range_summary", "version": "0.1.0",'
        write(unit, '(a)') '"input_fingerprint": "'//fingerprint//'",'
        write(unit, '(a)') '"directory": '//json_member(metadata, 'directory')//','
        write(unit, '(a)') '"istep": '//json_member(metadata, 'istep')//','
        write(unit, '(a,i0,a,i0,a,i0,a)') '"shape": [', dims(1), ',', dims(2), ',', dims(3), '],'
        write(unit, '(a,i0,a)') '"ranks": ', ranks, ','
        write(unit, '(a,es24.16e3,a)') '"wall_time_seconds": ', &
            real(clock_end-clock_start, real64)/real(clock_rate, real64), ','
        write(unit, '(a)') '"settings": {'
        write(unit, '(a,i0,a)') '"ispec": ', config%ispec, ','
        write(unit, '(a,es24.16e3,a)') '"dt": ', config%dt, ','
        write(unit, '(a,i0,a)') '"max_step": ', config%max_step, ','
        write(unit, '(a)') '"use_electric_field": '//json_member(metadata, 'use_electric_field')//','
        write(unit, '(a)') '"use_magnetic_field": '//json_member(metadata, 'use_magnetic_field')//','
        write(unit, '(a)') '"use_adaptive_dt": '//merge('true ', 'false', config%use_adaptive_dt)//','
        write(unit, '(a,es24.16e3,a)') '"coverage_sigma": ', config%coverage_sigma, ','
        write(unit, '(a,es24.16e3,a)') '"safety_factor": ', config%safety_factor, ','
        write(unit, '(a,i0,a)') '"source_samples_per_cell": ', config%source_samples_per_cell, ','
        write(unit, '(a,i0,a)') '"minimum_count": ', config%minimum_count, ','
        write(unit, '(a,i0,a)') '"max_probability_types": ', config%max_probability_types, ','
        write(unit, '(a,i0,a)') '"accumulator_cache_size": ', config%accumulator_cache_size, ','
        write(unit, '(a,i0,a)') '"n_threads": ', config%n_threads, '},'
        write(unit, '(a)') '"status_counts": {'
        do code = 0, 3
            write(unit, '(a,i0,a,i0)', advance='no') '"', code, '": ', &
                count(integers(:, :, :, 2) == code, kind=int64)
            if (code < 3) write(unit, '(a)', advance='no') ','
        end do
        write(unit, '(a)', iostat=ios) '}}'
        if (ios /= 0) call batch_fail('cannot write range summary')
        close(unit, iostat=ios)
        if (ios /= 0) call batch_fail('cannot close range summary')
        call atomic_replace(summary//'.tmp', summary)
        print '(a)', 'vdsolverf-batch: velocity_range computed'
    end subroutine
end module

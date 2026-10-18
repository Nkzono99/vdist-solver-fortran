module m_npy
    !! Npy v1/v2, C-order, little-endian f8/i4/i8. Values are flat in file order.
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use m_batch_support, only: batch_fail, number
    implicit none
    private
    public :: read_npy, write_npy
    interface read_npy
        module procedure read_real, read_i32, read_i64
    end interface
    interface write_npy
        module procedure write_real, write_i32, write_i64
    end interface
contains
    function shape_size(shape) result(n)
        integer(int64), intent(in) :: shape(:)
        integer(int64) :: n
        integer :: i
        if (size(shape) < 1 .or. size(shape) > 4) call batch_fail('npy supports one to four dimensions')
        n = 1
        do i = 1, size(shape)
            if (shape(i) < 0) call batch_fail('negative npy dimension')
            if (shape(i) /= 0) then
                if (n > huge(n)/shape(i)) call batch_fail('npy shape product overflow')
            end if
            n = n*shape(i)
        end do
    end function

    function header_member(header, key) result(value)
        character(*), intent(in) :: header, key
        character(:), allocatable :: value
        integer :: p, q
        p = index(header, "'"//key//"'")
        if (p == 0) p = index(header, '"'//key//'"')
        if (p == 0) call batch_fail('missing npy header key '//key)
        p = p + len(key) + 2
        do while (p <= len(header))
            if (header(p:p) /= ' ') exit
            p = p + 1
        end do
        if (p > len(header)) call batch_fail('truncated npy header')
        if (header(p:p) /= ':') call batch_fail('invalid npy header key '//key)
        p = p + 1
        do while (p <= len(header))
            if (header(p:p) /= ' ') exit
            p = p + 1
        end do
        if (p > len(header)) call batch_fail('truncated npy header')
        q = p
        if (header(p:p) == '(') then
            q = index(header(p:), ')') + p - 1
            if (q < p) call batch_fail('unterminated npy shape')
        else if (header(p:p) == "'" .or. header(p:p) == '"') then
            q = index(header(p+1:), header(p:p)) + p
            if (q == p) call batch_fail('unterminated npy dtype')
        else
            do while (q <= len(header))
                if (index(',} ', header(q:q)) /= 0) exit
                q = q + 1
            end do
            q = q - 1
        end if
        value = header(p:q)
    end function

    subroutine open_read(path, dtype, itemsize, unit, shape, n)
        character(*), intent(in) :: path, dtype
        integer, intent(in) :: itemsize
        integer, intent(out) :: unit
        integer(int64), allocatable, intent(out) :: shape(:)
        integer(int64), intent(out) :: n
        character(8) :: magic
        character(4) :: bytes
        character(:), allocatable :: header, raw
        integer :: ios, major, lengthbytes, i, p, start, nd
        integer(int64) :: headerlength, filesize, pos, dims(4)
        open(newunit=unit, file=path, status='old', access='stream', form='unformatted', &
            convert='little_endian', iostat=ios)
        if (ios /= 0) call batch_fail('cannot open npy '//path)
        read(unit, iostat=ios) magic
        if (ios /= 0) call batch_fail('truncated npy magic '//path)
        if (magic(:6) /= achar(147)//'NUMPY') call batch_fail('invalid npy magic '//path)
        major = iachar(magic(7:7))
        if ((major /= 1 .and. major /= 2) .or. iachar(magic(8:8)) /= 0) &
            call batch_fail('only npy v1.0/v2.0 supported: '//path)
        lengthbytes = 2*major
        read(unit, iostat=ios) bytes(:lengthbytes)
        if (ios /= 0) call batch_fail('truncated npy header length '//path)
        headerlength = 0
        do i = 1, lengthbytes
            headerlength = headerlength + ishft(int(iachar(bytes(i:i)), int64), 8*(i-1))
        end do
        if (headerlength < 1 .or. headerlength > 1048576) call batch_fail('invalid npy header length '//path)
        allocate(character(int(headerlength)) :: header)
        read(unit, iostat=ios) header
        if (ios /= 0) call batch_fail('truncated npy header '//path)
        if (header(len(header):) /= achar(10)) call batch_fail('npy header must end with newline: '//path)
        raw = header_member(header, 'descr')
        if (raw /= "'"//dtype//"'" .and. raw /= '"'//dtype//'"') &
            call batch_fail('npy dtype must be '//dtype//' (little-endian): '//path)
        raw = header_member(header, 'fortran_order')
        if (raw /= 'False') call batch_fail('npy Fortran-order arrays are unsupported: '//path)
        raw = header_member(header, 'shape')
        if (raw(1:1) /= '(') call batch_fail('invalid npy shape '//path)
        p = 2
        nd = 0
        do
            do while (p < len(raw))
                if (raw(p:p) /= ' ') exit
                p = p + 1
            end do
            if (p >= len(raw)) exit
            start = p
            do while (p < len(raw))
                if (index('0123456789', raw(p:p)) == 0) exit
                p = p + 1
            end do
            if (p == start .or. nd == 4) call batch_fail('invalid npy dimension '//path)
            nd = nd + 1
            read(raw(start:p-1), *, iostat=ios) dims(nd)
            if (ios /= 0) call batch_fail('npy dimension overflow '//path)
            do while (p < len(raw))
                if (raw(p:p) /= ' ') exit
                p = p + 1
            end do
            if (p >= len(raw)) exit
            if (raw(p:p) /= ',') call batch_fail('invalid npy shape separator '//path)
            p = p + 1
        end do
        shape = dims(:nd)
        n = shape_size(shape)
        inquire(unit=unit, size=filesize, pos=pos)
        if (n > (huge(n)-pos)/itemsize) call batch_fail('npy payload size overflow '//path)
        if (filesize /= pos-1+n*itemsize) call batch_fail('npy payload length mismatch '//path)
    end subroutine

    subroutine open_write(path, dtype, unit, shape, n, version)
        character(*), intent(in) :: path, dtype
        integer, intent(out) :: unit
        integer(int64), intent(in) :: shape(:), n
        integer, intent(in), optional :: version
        integer :: i, major, lengthbytes, padding, ios
        character(4) :: bytes
        character(32) :: buffer
        character(:), allocatable :: header
        if (shape_size(shape) /= n) call batch_fail('npy write shape does not match data: '//path)
        major = 1
        if (present(version)) major = version
        if (major /= 1 .and. major /= 2) call batch_fail('npy write version must be 1 or 2')
        header = "{'descr': '"//dtype//"', 'fortran_order': False, 'shape': ("
        do i = 1, size(shape)
            write(buffer, '(i0)') shape(i)
            header = header//trim(buffer)//','
        end do
        header = header//'), }'
        lengthbytes = 2*major
        padding = modulo(-(8+lengthbytes+len(header)+1), 64)
        header = header//repeat(' ', padding)//achar(10)
        if (major == 1 .and. len(header) > 65535) call batch_fail('npy v1 header too long')
        do i = 1, lengthbytes
            bytes(i:i) = achar(iand(ishft(len(header), -8*(i-1)), 255))
        end do
        open(newunit=unit, file=path, status='replace', access='stream', form='unformatted', &
            convert='little_endian', iostat=ios)
        if (ios /= 0) call batch_fail('cannot write npy '//path)
        write(unit, iostat=ios) achar(147)//'NUMPY'//achar(major)//achar(0), bytes(:lengthbytes), header
        if (ios /= 0) call batch_fail('cannot write npy header '//path)
    end subroutine

    subroutine read_real(path, values, shape)
        character(*), intent(in) :: path
        real(real64), allocatable, intent(out) :: values(:)
        integer(int64), allocatable, intent(out) :: shape(:)
        integer :: unit, ios
        integer(int64) :: n
        call open_read(path, '<f8', 8, unit, shape, n)
        allocate(values(n), stat=ios)
        if (ios /= 0) call batch_fail('cannot allocate npy '//path)
        read(unit, iostat=ios) values
        close(unit)
        if (ios /= 0) call batch_fail('cannot read npy payload '//path)
    end subroutine

    subroutine read_i32(path, values, shape)
        character(*), intent(in) :: path
        integer(int32), allocatable, intent(out) :: values(:)
        integer(int64), allocatable, intent(out) :: shape(:)
        integer :: unit, ios
        integer(int64) :: n
        call open_read(path, '<i4', 4, unit, shape, n)
        allocate(values(n), stat=ios)
        if (ios /= 0) call batch_fail('cannot allocate npy '//path)
        read(unit, iostat=ios) values
        close(unit)
        if (ios /= 0) call batch_fail('cannot read npy payload '//path)
    end subroutine

    subroutine read_i64(path, values, shape)
        character(*), intent(in) :: path
        integer(int64), allocatable, intent(out) :: values(:)
        integer(int64), allocatable, intent(out) :: shape(:)
        integer :: unit, ios
        integer(int64) :: n
        call open_read(path, '<i8', 8, unit, shape, n)
        allocate(values(n), stat=ios)
        if (ios /= 0) call batch_fail('cannot allocate npy '//path)
        read(unit, iostat=ios) values
        close(unit)
        if (ios /= 0) call batch_fail('cannot read npy payload '//path)
    end subroutine

    subroutine write_real(path, values, shape, version)
        character(*), intent(in) :: path
        real(real64), intent(in) :: values(:)
        integer(int64), intent(in) :: shape(:)
        integer, intent(in), optional :: version
        integer :: unit, ios, closing
        call open_write(path, '<f8', unit, shape, size(values, kind=int64), version)
        write(unit, iostat=ios) values
        close(unit, iostat=closing)
        if (ios /= 0 .or. closing /= 0) call batch_fail('cannot finish npy '//path)
    end subroutine

    subroutine write_i32(path, values, shape, version)
        character(*), intent(in) :: path
        integer(int32), intent(in) :: values(:)
        integer(int64), intent(in) :: shape(:)
        integer, intent(in), optional :: version
        integer :: unit, ios, closing
        call open_write(path, '<i4', unit, shape, size(values, kind=int64), version)
        write(unit, iostat=ios) values
        close(unit, iostat=closing)
        if (ios /= 0 .or. closing /= 0) call batch_fail('cannot finish npy '//path)
    end subroutine

    subroutine write_i64(path, values, shape, version)
        character(*), intent(in) :: path
        integer(int64), intent(in) :: values(:)
        integer(int64), intent(in) :: shape(:)
        integer, intent(in), optional :: version
        integer :: unit, ios, closing
        call open_write(path, '<i8', unit, shape, size(values, kind=int64), version)
        write(unit, iostat=ios) values
        close(unit, iostat=closing)
        if (ios /= 0 .or. closing /= 0) call batch_fail('cannot finish npy '//path)
    end subroutine
end module

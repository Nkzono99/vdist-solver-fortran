module m_batch_support
    !! Small POSIX file helpers and a reader for field-bundle JSON metadata.
    use, intrinsic :: iso_fortran_env, only: error_unit, int64
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    implicit none
    private
    public :: batch_fail, read_text, json_member, json_string, json_integer, json_logical
    public :: make_directory, atomic_replace, remove_file, path_join, parent_path, number
    public :: set_batch_abort
    abstract interface
        subroutine abort_handler()
        end subroutine
    end interface
    procedure(abort_handler), pointer :: abort_callback => null()
    interface
        function c_rename(old, new) bind(c, name='rename') result(code)
            import c_char, c_int
            character(c_char), intent(in) :: old(*), new(*)
            integer(c_int) :: code
        end function
    end interface
contains
    subroutine batch_fail(message)
        character(*), intent(in) :: message
        write(error_unit, '(a)') 'vdsolverf-batch: '//message
        if (associated(abort_callback)) call abort_callback()
        error stop 1
    end subroutine

    subroutine set_batch_abort(handler)
        procedure(abort_handler) :: handler
        abort_callback => handler
    end subroutine

    function number(value) result(s)
        integer, intent(in) :: value
        character(:), allocatable :: s
        character(32) :: buffer
        write(buffer, '(i0)') value
        s = trim(buffer)
    end function

    function parent_path(path) result(parent)
        character(*), intent(in) :: path
        character(:), allocatable :: parent
        integer :: slash
        slash = index(trim(path), '/', back=.true.)
        parent = '.'
        if (slash == 1) parent = '/'
        if (slash > 1) parent = path(:slash-1)
    end function

    function path_join(base, path) result(joined)
        character(*), intent(in) :: base, path
        character(:), allocatable :: joined
        if (len_trim(path) == 0) call batch_fail('empty path')
        joined = trim(path)
        if (path(1:1) /= '/') joined = trim(base)//'/'//trim(path)
    end function

    function shell_quote(s) result(quoted)
        character(*), intent(in) :: s
        character(:), allocatable :: quoted
        integer :: i
        quoted = "'"
        do i = 1, len_trim(s)
            if (s(i:i) == "'") then
                quoted = quoted//"'"//'"'//"'"//'"'//"'"
            else
                quoted = quoted//s(i:i)
            end if
        end do
        quoted = quoted//"'"
    end function

    subroutine make_directory(path)
        character(*), intent(in) :: path
        integer :: status, cmdstatus
        call execute_command_line('mkdir -p -- '//shell_quote(path), exitstat=status, cmdstat=cmdstatus)
        if (cmdstatus /= 0) call batch_fail('cannot launch mkdir for '//path)
        if (status /= 0) call batch_fail('cannot create directory '//path)
    end subroutine

    subroutine atomic_replace(old, new)
        character(*), intent(in) :: old, new
        if (c_rename(old//c_null_char, new//c_null_char) /= 0) &
            call batch_fail('cannot rename '//old//' to '//new)
    end subroutine

    subroutine remove_file(path)
        character(*), intent(in) :: path
        logical :: exists
        integer :: unit, ios
        inquire(file=path, exist=exists)
        if (.not. exists) return
        open(newunit=unit, file=path, status='old', iostat=ios)
        if (ios /= 0) call batch_fail('cannot remove '//path)
        close(unit, status='delete', iostat=ios)
        if (ios /= 0) call batch_fail('cannot remove '//path)
    end subroutine

    function read_text(path) result(s)
        character(*), intent(in) :: path
        character(:), allocatable :: s
        integer :: unit, ios
        integer(int64) :: length
        inquire(file=path, size=length, iostat=ios)
        if (ios /= 0 .or. length < 0 .or. length > 1048576) call batch_fail('invalid text file '//path)
        allocate(character(int(length)) :: s)
        open(newunit=unit, file=path, status='old', access='stream', form='unformatted', iostat=ios)
        if (ios /= 0) call batch_fail('cannot open '//path)
        read(unit, iostat=ios) s
        close(unit)
        if (ios /= 0) call batch_fail('cannot read '//path)
    end function

    subroutine whitespace(s, p)
        character(*), intent(in) :: s
        integer, intent(inout) :: p
        do while (p <= len(s))
            if (index(' '//achar(9)//achar(10)//achar(13), s(p:p)) == 0) exit
            p = p + 1
        end do
    end subroutine

    recursive subroutine skip_value(s, p, depth)
        character(*), intent(in) :: s
        integer, intent(inout) :: p
        integer, intent(in) :: depth
        character :: closing
        if (depth > 32) call batch_fail('bundle JSON is nested too deeply')
        call whitespace(s, p)
        if (p > len(s)) call batch_fail('truncated bundle JSON')
        select case (s(p:p))
        case ('"')
            p = p + 1
            do while (p <= len(s))
                if (s(p:p) == '"') exit
                if (iachar(s(p:p)) < 32) call batch_fail('control character in bundle JSON string')
                if (s(p:p) == achar(92)) p = p + 1
                p = p + 1
            end do
            if (p > len(s)) call batch_fail('unterminated bundle JSON string')
            p = p + 1
        case ('[', '{')
            closing = ']'
            if (s(p:p) == '{') closing = '}'
            p = p + 1
            call whitespace(s, p)
            if (p > len(s)) call batch_fail('truncated bundle JSON container')
            if (s(p:p) /= closing) then
                do
                    if (closing == '}') then
                        if (s(p:p) /= '"') call batch_fail('expected bundle JSON key')
                        call skip_value(s, p, depth+1)
                        call whitespace(s, p)
                        if (p > len(s)) call batch_fail('truncated bundle JSON object')
                        if (s(p:p) /= ':') call batch_fail('expected colon in bundle JSON')
                        p = p + 1
                    end if
                    call skip_value(s, p, depth+1)
                    call whitespace(s, p)
                    if (p > len(s)) call batch_fail('truncated bundle JSON container')
                    if (s(p:p) == closing) exit
                    if (s(p:p) /= ',') call batch_fail('expected comma in bundle JSON')
                    p = p + 1
                    call whitespace(s, p)
                    if (p > len(s)) call batch_fail('truncated bundle JSON container')
                end do
            end if
            p = p + 1
        case default
            block
                integer :: start
                start = p
                do while (p <= len(s))
                    if (index(',]} '//achar(9)//achar(10)//achar(13), s(p:p)) /= 0) exit
                    p = p + 1
                end do
                if (p == start) call batch_fail('invalid bundle JSON value')
            end block
        end select
    end subroutine

    function json_member(s, key) result(value)
        character(*), intent(in) :: s, key
        character(:), allocatable :: value
        integer :: p, start, finish
        character(:), allocatable :: name
        logical :: found
        p = 1
        found = .false.
        call whitespace(s, p)
        if (p > len(s)) call batch_fail('empty bundle JSON')
        if (s(p:p) /= '{') call batch_fail('expected bundle JSON object')
        p = p + 1
        do
            call whitespace(s, p)
            if (p > len(s)) call batch_fail('truncated bundle JSON')
            if (s(p:p) == '}') exit
            if (s(p:p) /= '"') call batch_fail('expected bundle JSON key')
            start = p
            call skip_value(s, p, 0)
            name = s(start+1:p-2)
            call whitespace(s, p)
            if (p > len(s)) call batch_fail('truncated bundle JSON')
            if (s(p:p) /= ':') call batch_fail('expected colon in bundle JSON')
            p = p + 1
            call whitespace(s, p)
            start = p
            call skip_value(s, p, 0)
            finish = p - 1
            if (name == key) then
                if (found) call batch_fail('duplicate bundle JSON key '//key)
                value = s(start:finish)
                found = .true.
            end if
            call whitespace(s, p)
            if (p > len(s)) call batch_fail('truncated bundle JSON')
            if (s(p:p) == '}') exit
            if (s(p:p) /= ',') call batch_fail('expected comma in bundle JSON')
            p = p + 1
        end do
        p = p + 1
        call whitespace(s, p)
        if (p <= len(s)) call batch_fail('trailing content in bundle JSON')
        if (.not. found) call batch_fail('missing bundle JSON key '//key)
    end function

    function json_string(s, key) result(value)
        character(*), intent(in) :: s, key
        character(:), allocatable :: value, raw
        raw = json_member(s, key)
        if (len(raw) < 2) call batch_fail('expected string for '//key)
        if (raw(1:1) /= '"' .or. raw(len(raw):) /= '"') call batch_fail('expected string for '//key)
        value = raw(2:len(raw)-1)
    end function

    function json_integer(s, key) result(value)
        character(*), intent(in) :: s, key
        integer :: value, ios, i
        character(:), allocatable :: raw
        raw = json_member(s, key)
        do i = 1, len(raw)
            if (index('0123456789', raw(i:i)) == 0) call batch_fail('expected unsigned integer for '//key)
        end do
        read(raw, *, iostat=ios) value
        if (ios /= 0) call batch_fail('invalid integer for '//key)
    end function

    function json_logical(s, key) result(value)
        character(*), intent(in) :: s, key
        logical :: value
        character(:), allocatable :: raw
        raw = json_member(s, key)
        select case (raw)
        case ('true')
            value = .true.
        case ('false')
            value = .false.
        case default
            call batch_fail('expected boolean for '//key)
        end select
    end function
end module

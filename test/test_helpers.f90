module m_test_helpers
    !! Shared assertion helpers for the Fortran test suite.
    !!
    !! Each helper prints a FAIL line with actual/expected values and calls
    !! `error stop 1` to fail the `fpm test` run. The `label` argument is
    !! printed on both success and failure so it is easy to see which
    !! assertions ran.

    implicit none

    private
    public :: assert_close
    public :: assert_close_vec
    public :: assert_equal_int
    public :: assert_true

contains

    subroutine assert_close(label, actual, expected, tolerance)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: actual, expected
        double precision, intent(in), optional :: tolerance

        double precision :: tol

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 1d-12
        end if

        if (abs(actual - expected) > tol) then
            print *, "FAILED:", trim(label)
            print *, "  actual  =", actual
            print *, "  expected=", expected
            print *, "  diff    =", actual - expected
            error stop 1
        end if
    end subroutine

    subroutine assert_close_vec(label, actual, expected, tolerance)
        character(len=*), intent(in) :: label
        double precision, intent(in) :: actual(:), expected(:)
        double precision, intent(in), optional :: tolerance

        double precision :: tol
        integer :: i

        if (size(actual) /= size(expected)) then
            print *, "FAILED:", trim(label), " (size mismatch)"
            error stop 1
        end if

        if (present(tolerance)) then
            tol = tolerance
        else
            tol = 1d-12
        end if

        do i = 1, size(actual)
            if (abs(actual(i) - expected(i)) > tol) then
                print *, "FAILED:", trim(label), " at index ", i
                print *, "  actual  =", actual(i)
                print *, "  expected=", expected(i)
                error stop 1
            end if
        end do
    end subroutine

    subroutine assert_equal_int(label, actual, expected)
        character(len=*), intent(in) :: label
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            print *, "FAILED:", trim(label)
            print *, "  actual  =", actual
            print *, "  expected=", expected
            error stop 1
        end if
    end subroutine

    subroutine assert_true(label, condition)
        character(len=*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            print *, "FAILED:", trim(label), " (expected .true.)"
            error stop 1
        end if
    end subroutine

end module

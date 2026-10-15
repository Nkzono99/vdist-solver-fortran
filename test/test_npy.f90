program test_npy
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use m_npy, only: read_npy, write_npy
    use m_batch_support, only: remove_file
    use m_test_helpers, only: assert_true
    implicit none
    real(real64) :: reals(24)
    real(real64), allocatable :: restored(:)
    integer(int32) :: ints32(24)
    integer(int32), allocatable :: restored32(:)
    integer(int64) :: ints64(24)
    integer(int64), allocatable :: restored64(:), dims(:)
    integer :: i, version
    reals = [(real(i-12, real64)/7d0, i=1, 24)]
    ints32 = [(int(i-12, int32), i=1, 24)]
    ints64 = int(ints32, int64)*2_int64**40
    do version = 1, 2
        call write_npy('test-npy-real.npy', reals, [2_int64, 3_int64, 4_int64], version)
        call read_npy('test-npy-real.npy', restored, dims)
        call assert_true('npy real shape', all(dims == [2, 3, 4]))
        call assert_true('npy real exact payload', all(restored == reals))
        call write_npy('test-npy-i32.npy', ints32, [24_int64], version)
        call read_npy('test-npy-i32.npy', restored32, dims)
        call assert_true('npy i32 exact', size(dims) == 1 .and. dims(1) == 24 .and. all(restored32 == ints32))
        call write_npy('test-npy-i64.npy', ints64, [4_int64, 6_int64], version)
        call read_npy('test-npy-i64.npy', restored64, dims)
        call assert_true('npy i64 exact (no narrowing)', all(restored64 == ints64))
    end do
    call write_npy('test-npy-real.npy', reals(:0), [0_int64])
    call read_npy('test-npy-real.npy', restored, dims)
    call assert_true('empty npy payload', size(restored) == 0 .and. dims(1) == 0)
    call remove_file('test-npy-real.npy')
    call remove_file('test-npy-i32.npy')
    call remove_file('test-npy-i64.npy')
    print *, 'test_npy: all tests passed.'
end program

program vdsolverf_batch
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use m_vdsolverf_batch
    use m_batch_support
#ifdef VDSOLVERF_MPI
    use mpi, only: MPI_STATUS_SIZE, MPI_THREAD_FUNNELED, MPI_SUCCESS, MPI_COMM_WORLD, &
        MPI_ANY_SOURCE, MPI_SOURCE, MPI_INTEGER, MPI_Init_thread, MPI_Comm_rank, MPI_Comm_size, &
        MPI_Barrier, MPI_Finalize, MPI_Abort
#endif
    implicit none
    type(t_batch) :: batch
    character(:), allocatable :: argument, input
    logical :: resume
    integer :: rank, ranks, i, chunk, computed, skipped, next
    integer, allocatable :: work(:)
    integer(int64) :: start_clock, end_clock, clock_rate
#ifdef VDSOLVERF_MPI
    integer :: ierr, provided, workers, source, request, recv_status(MPI_STATUS_SIZE)
    ! Every message is one default integer. Explicit interfaces also support
    ! legacy vendor `mpi` modules that omit choice-buffer interfaces (GCC 8).
    interface
        subroutine MPI_Send(buffer, count, datatype, destination, tag, communicator, error)
            integer, intent(in) :: buffer, count, datatype, destination, tag, communicator
            integer, intent(out) :: error
        end subroutine
        subroutine MPI_Recv(buffer, count, datatype, source, tag, communicator, status, error)
            import MPI_STATUS_SIZE
            integer, intent(out) :: buffer, status(MPI_STATUS_SIZE), error
            integer, intent(in) :: count, datatype, source, tag, communicator
        end subroutine
    end interface
#endif
    rank = 0
    ranks = 1
#ifdef VDSOLVERF_MPI
    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
    if (ierr /= MPI_SUCCESS) error stop 1
    call set_batch_abort(abort_mpi)
    if (provided < MPI_THREAD_FUNNELED) call batch_fail('MPI does not support MPI_THREAD_FUNNELED')
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, ranks, ierr)
#endif
    call system_clock(start_clock, clock_rate)
    resume = .true.
    do i = 1, command_argument_count()
        call get_argument(i, argument)
        select case (argument)
        case ('--resume')
            resume = .true.
        case ('--no-resume')
            resume = .false.
        case ('--help', '-h')
            if (rank == 0) print '(a)', 'Usage: vdsolverf-batch [--resume|--no-resume] batch.nml'
#ifdef VDSOLVERF_MPI
            call MPI_Finalize(ierr)
#endif
            stop
        case default
            if (len(argument) == 0) call batch_fail('empty command argument')
            if (argument(1:1) == '-') call batch_fail('unknown option '//argument)
            if (allocated(input)) call batch_fail('only one batch namelist may be specified')
            input = argument
        end select
    end do
    if (.not. allocated(input)) call batch_fail('Usage: vdsolverf-batch [--no-resume] batch.nml')
    ! Every rank reads and constructs its simulator once, including the dispatcher.
    call open_batch(input, batch)
    computed = 0
    skipped = 0
    if (rank == 0) then
        call prepare_work(batch, resume, work, computed, skipped)
    end if
#ifdef VDSOLVERF_MPI
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    if (ranks > 1) then
        if (rank == 0) then
            ! Rank zero never enters the costly octree: respond promptly to all workers.
            next = 1
            workers = ranks-1
            do while (workers > 0)
                call MPI_Recv(request, 1, MPI_INTEGER, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, recv_status, ierr)
                source = recv_status(MPI_SOURCE)
                chunk = -1
                if (next <= computed) then
                    chunk = work(next)
                    next = next + 1
                else
                    workers = workers-1
                end if
                call MPI_Send(chunk, 1, MPI_INTEGER, source, 2, MPI_COMM_WORLD, ierr)
            end do
        else
            do
                request = 0
                call MPI_Send(request, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, ierr)
                call MPI_Recv(chunk, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, recv_status, ierr)
                if (chunk < 0) exit
                call compute_chunk(batch, chunk)
            end do
        end if
    else
#endif
        do next = 1, computed
            call compute_chunk(batch, work(next))
        end do
#ifdef VDSOLVERF_MPI
    end if
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif
    if (rank == 0) then
        call system_clock(end_clock)
        call write_summary(batch, real(end_clock-start_clock, real64)/real(clock_rate, real64), &
            ranks, computed, skipped)
        print '(a,i0,a,i0)', 'vdsolverf-batch: chunks computed=', computed, ', skipped=', skipped
    end if
    call close_batch(batch)
#ifdef VDSOLVERF_MPI
    call MPI_Finalize(ierr)
#endif
contains
    subroutine get_argument(index, value)
        integer, intent(in) :: index
        character(:), allocatable, intent(out) :: value
        integer :: length
        call get_command_argument(index, length=length)
        allocate(character(length) :: value)
        call get_command_argument(index, value=value)
    end subroutine
#ifdef VDSOLVERF_MPI
    subroutine abort_mpi()
        integer :: error
        call MPI_Abort(MPI_COMM_WORLD, 1, error)
    end subroutine
#endif
end program

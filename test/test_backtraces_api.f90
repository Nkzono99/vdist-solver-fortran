program test_backtraces_api
    !! Exercise the C-facing `get_backtraces` entry point directly (issue #44:
    !! the optimised Intel build crashed inside this routine when called from
    !! Python). A 2x1x1 zero-field box with a Maxwellian plane at x=2.
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_double
    use m_emses_solver, only: get_backtraces
    use m_test_helpers, only: assert_close, assert_close_vec, assert_equal_int, assert_true

    implicit none

    double precision, parameter :: pi = acos(-1d0)
    character(len=*), parameter :: filename = 'test-backtraces-api.inp'
    character(kind=c_char) :: cpath(len(filename))
    integer :: unit, i

    open (newunit=unit, file=filename, status='replace')
    write (unit, '(a)') '&esorem emflag=0 /'
    write (unit, '(a)') '&plasma wp=10*0, wc=0, phiz=0, phixy=0 /'
    write (unit, '(a)') '&tmgrid nx=2, ny=1, nz=1, dt=0.1 /'
    write (unit, '(a)') '&system nspec=1, npbnd=2,0,0 /'
    write (unit, '(a)') '&intp qm=10*1, peth=1, path=1, vdri=10*0, vdthz=10*0, vdthxy=10*0 /'
    write (unit, '(a)') "&ptcond boundary_type='none' /"
    write (unit, '(a)') '&emissn use_raycast=.false., nflag_emit=10*0, nepl=10*0 /'
    close (unit)
    do i = 1, len(filename)
        cpath(i) = filename(i:i)
    end do

    call run_case(npcls=1, max_step=3, output_interval=1, n_threads=1)
    call run_case(npcls=3, max_step=3, output_interval=1, n_threads=1)
    call run_case(npcls=3, max_step=7, output_interval=2, n_threads=2)

    open (newunit=unit, file=filename, status='old')
    close (unit, status='delete')
    print *, "test_backtraces_api: all tests passed."

contains

    subroutine run_case(npcls, max_step, output_interval, n_threads)
        integer, intent(in) :: npcls, max_step, output_interval, n_threads
        integer(c_int), parameter :: lx = 2, ly = 1, lz = 1
        real(c_double) :: ebvalues(9, lx + 1, ly + 1, lz + 1)
        real(c_double) :: positions(3, npcls), velocities(3, npcls)
        real(c_double) :: ts((max_step - 1)/output_interval + 2, npcls)
        real(c_double) :: probabilities(npcls)
        real(c_double) :: traced_positions(3, (max_step - 1)/output_interval + 2, npcls)
        real(c_double) :: traced_velocities(3, (max_step - 1)/output_interval + 2, npcls)
        integer(c_int) :: last_steps(npcls)
        integer :: ipcl, nout

        ebvalues = 0d0
        do ipcl = 1, npcls
            positions(:, ipcl) = [1d0, 0.5d0, 0.5d0]
            velocities(:, ipcl) = [-1d0 - 0.5d0*(ipcl - 1), 0d0, 0d0]
        end do
        nout = (max_step - 1)/output_interval + 2

        call get_backtraces(cpath, int(len(filename), c_int), lx, ly, lz, ebvalues, 1_c_int, &
                            int(npcls, c_int), positions, velocities, 2d0, int(max_step, c_int), &
                            int(output_interval, c_int), 0_c_int, 10_c_int, &
                            ts, probabilities, traced_positions, traced_velocities, last_steps, &
                            int(n_threads, c_int))

        do ipcl = 1, npcls
            call assert_close("get_backtraces returns the Maxwellian at the plane", probabilities(ipcl), &
                              exp(-0.5d0*velocities(1, ipcl)**2)/(2d0*pi)**1.5d0)
            call assert_true("last step is within the trace buffer", &
                             last_steps(ipcl) >= 1 .and. last_steps(ipcl) <= nout)
            call assert_close_vec("final trace position is on the plane", traced_positions(:, nout, ipcl), &
                                  [2d0, 0.5d0, 0.5d0])
            call assert_close("final trace time is distance over speed", ts(nout, ipcl), &
                              1d0/abs(velocities(1, ipcl)))
        end do
    end subroutine
end program

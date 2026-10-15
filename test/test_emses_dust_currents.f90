program test_emses_dust_currents
    use, intrinsic :: iso_c_binding, only: c_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_emses_dust_currents, only: evaluate_dust_currents_impl
    use m_allcom, only: qm
    use m_simulator, only: t_ESSimulator
    use m_test_helpers, only: assert_true, assert_equal_int
    implicit none
    type(t_ESSimulator) :: simulator
    character(len=*), parameter :: filename = 'test-emses-dust-currents.inp'
    character(kind=c_char) :: cpath(len(filename))
    real(8) :: fields(9, 3, 3, 3), points(3, 2), bounds(6, 2), psi(5)
    real(8) :: flux(5, 2), error(5, 2), serial_flux(5, 2), serial_error(5, 2), expected(5), pi
    integer :: status(2), samples(2), terminals(2), serial_status(2), serial_samples(2), serial_terminals(2)
    integer :: unit, i

    open(newunit=unit, file=filename, status='replace')
    write(unit, *) '&esorem emflag=0 /'
    write(unit, *) '&plasma wp=10*0, wc=0, phiz=0, phixy=0 /'
    write(unit, *) '&tmgrid nx=2, ny=2, nz=2, dt=0.1 /'
    write(unit, *) '&system nspec=1, npbnd=30*2 /'
    ! path/peth are one-dimensional thermal standard deviations, sqrt(kT/m).
    write(unit, *) '&intp qm=10*1, spa=10*0, spe=10*0, speth=10*0, ', &
                   'vdri=10*0, vdthz=10*0, vdthxy=10*0, peth=1, path=1 /'
    write(unit, *) "&ptcond boundary_type='none' /"
    write(unit, *) '&emissn nepl=10*0, nflag_emit=10*0 /'
    close(unit)
    do i = 1, len(filename)
        cpath(i) = filename(i:i)
    end do
    fields = 0d0
    simulator = create_simulator(cpath, len(filename), 2, 2, 2, fields, 1, 100)
    points(:, 1) = [1d0, 1d0, 1d0]
    points(:, 2) = [0.75d0, 1d0, 1d0]
    bounds(:, 1) = [-4.5d0, 4.5d0, -4.5d0, 4.5d0, -4.5d0, 4.5d0]
    bounds(:, 2) = bounds(:, 1)
    psi = [-1d0, -0.25d0, 0d0, 0.25d0, 1d0]
    pi = acos(-1d0)
    expected = sqrt(8d0/pi)*merge(exp(-psi), 1d0 - psi, psi > 0d0)
    call evaluate(1)
    print *, 'dust analytic relative errors:', abs(flux(:, 1)/expected - 1d0)
    call assert_true('default octree OML Maxwellian moments within one percent', &
                     all(abs(flux(:, 1)/expected - 1d0) < 1d-2))
    call assert_true('dust errors finite and nonnegative', all(ieee_is_finite(error)) .and. all(error >= 0d0))
    call assert_true('dust flux finite and nonnegative', all(ieee_is_finite(flux)) .and. all(flux >= 0d0))
    call assert_true('dust returns samples and terminal boxes', all(samples > 0) .and. all(terminals > 0))
    serial_flux = flux
    serial_error = error
    serial_status = status
    serial_samples = samples
    serial_terminals = terminals
    call evaluate(2)
    call assert_true('OpenMP dust flux equals serial bit for bit', all(flux == serial_flux))
    call assert_true('OpenMP dust errors equal serial bit for bit', all(error == serial_error))
    call assert_true('OpenMP dust diagnostics equal serial', all(status == serial_status) .and. &
                     all(samples == serial_samples) .and. all(terminals == serial_terminals))
    ! The attractive/repulsive branches must follow the namelist qm sign.
    qm(1) = -1d0
    psi = -psi
    call evaluate(1)
    call assert_true('negative qm with opposite psi preserves flux', all(flux == serial_flux))
    qm(1) = 1d0
    psi = -psi
    bounds(2, 1) = bounds(1, 1)
    call evaluate(1)
    call assert_equal_int('dust invalid bounds', status(1), 5)
    call assert_equal_int('invalid bounds have no samples', samples(1), 0)
    call assert_equal_int('invalid bounds have no terminal boxes', terminals(1), 0)
    call assert_true('invalid bounds have zero moments', all(flux(:, 1) == 0d0))
    call destroy_simulator(simulator)
    open(newunit=unit, file=filename, status='old')
    close(unit, status='delete')
    print *, 'test_emses_dust_currents: all tests passed.'
contains
    subroutine evaluate(threads)
        integer, intent(in) :: threads
        call evaluate_dust_currents_impl(simulator, 1, 2, points, bounds, 5, psi, 2d0, 32, 0, &
            [9, 9, 9], 6, 50000, 4096, 1d-4, 1d-4, 1.5d0, 3, &
            flux, error, status, samples, terminals, threads)
    end subroutine
end program

program test_thread_determinism
    !! Nonuniform E/B, staggered E, and enough spatial work for eight workers.
    !! No namelist I/O: also usable with ifx's serial -check uninit.
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int64
!$  use omp_lib, only: omp_set_num_threads, omp_set_dynamic, omp_get_num_threads
    use finbound, only: t_BoundaryList, new_BoundaryList, t_Boundary, t_PlaneXYZ, new_PlaneXYZ
    use m_allcom, only: qm
    use m_emses_field, only: t_EMSESFieldGrid, new_EMSESFieldGrid
    use m_simulator, only: t_ESSimulator, new_ESSimulator
    use m_particle, only: new_Particle
    use m_probabilities, only: tp_Probability, new_MaxwellianProbability
    use m_solver, only: t_Solver, t_ProbabilityRecord, new_Solver
    use m_emses_octree_probabilities, only: get_probabilities_octree_impl
    use m_emses_dust_currents, only: evaluate_dust_currents_impl
    use m_test_helpers, only: assert_true
    implicit none
    integer, parameter :: n = 19, ns = 4096*n, nb = 128*n
    integer, parameter :: dimensions(3) = [3, 2, 4]
    type results
        real(c_double) :: flux(4, n), error(4, n), velocities(3, ns), probabilities(ns)
        real(c_double) :: boxes(6, nb), pmin(nb), pmax(nb), direct(7, 36*n)
        integer(c_int) :: status(n), samples(n), terminals(n), si(ns), bi(nb), depth(nb), starts(nb), counts(nb)
        integer(c_int) :: oct_status(n), oct_samples(n), oct_leaves(n), total_samples, total_leaves
        logical :: terminal(nb)
    end type
    type(results), allocatable :: reference, current
    type(t_ESSimulator) :: simulator
    type(t_Solver) :: solver
    type(t_EMSESFieldGrid) :: field
    type(t_BoundaryList) :: boundaries
    type(tp_Probability) :: probabilities(1)
    type(t_PlaneXYZ), pointer :: plane
    class(t_Boundary), pointer :: boundary
    real(c_double) :: fields(9, 0:3, 0:2, 0:4), points(3, n), bounds(6, n), psi(4), v(3, 36)
    integer :: i, axis, side, repeat, ix, iy, iz, j

    fields = reshape([(1d-5*real(i-1, c_double), i=1, size(fields))], shape(fields))
    field = new_EMSESFieldGrid(3, 2, 4, fields(:6, :, :, :), fields(7:, :, :, :))
    boundaries = new_BoundaryList()
    do axis = 1, 3
        do side = 0, 1
            allocate(plane)
            plane = new_PlaneXYZ(axis, real(side*dimensions(axis), c_double))
            boundary => plane
            boundary%material%tag = 1
            call boundaries%add_boundary(boundary)
        end do
    end do
    allocate(probabilities(1)%ref, source=new_MaxwellianProbability([0d0, 0d0, 0d0], [1d0, 1d0, 1d0]))
    simulator = new_ESSimulator(3, 2, 4, [2, 2, 2], field, boundaries, probabilities)
    solver = new_Solver(simulator)
    qm(1) = 1d0
    do i = 1, n
        points(:, i) = [0.4d0+0.1d0*(i-1), 0.8d0, 1.7d0]
        bounds(:, i) = [-4d0, 4d0, -4d0, 4d0, -4d0, 4d0]
        if (mod(i-1, 5) == 0) bounds(:, i) = [-0.5d0, 0.5d0, -0.5d0, 0.5d0, -0.5d0, 0.5d0]
    end do
    bounds(2, n) = -5d0
    psi = [-1d0, 0d0, 0.25d0, 1d0]
    j = 0
    do iz = 0, 2
        do iy = 0, 3
            do ix = 0, 2
                j = j+1
                v(:, j) = [-4d0+4d0*ix, -4d0+8d0*iy/3d0, -4d0+4d0*iz]
            end do
        end do
    end do
    allocate(reference, current)
!$  call omp_set_dynamic(.false.)
!$  call omp_set_num_threads(8)
    !$omp parallel
    !$omp single
!$  call assert_true('eight OpenMP workers active', omp_get_num_threads() == 8)
    !$omp end single
    !$omp end parallel
    call evaluate(1)
    reference = current
    ! Two independent eight-thread evaluations must match the serial baseline
    ! and therefore each other, including every probability and diagnostic.
    do repeat = 1, 2
        call evaluate(8)
        call assert_true('dust flux bits', same_bits(reshape(current%flux, [4*n]), reshape(reference%flux, [4*n])))
        call assert_true('dust error bits', same_bits(reshape(current%error, [4*n]), reshape(reference%error, [4*n])))
        call assert_true('dust diagnostics', all(current%status == reference%status) .and. &
            all(current%samples == reference%samples) .and. all(current%terminals == reference%terminals))
        call assert_true('octree probability bits', same_bits(current%probabilities, reference%probabilities))
        call assert_true('octree velocity bits', same_bits(reshape(current%velocities, [3*ns]), &
            reshape(reference%velocities, [3*ns])))
        call assert_true('octree box bits', same_bits(reshape(current%boxes, [6*nb]), reshape(reference%boxes, [6*nb])))
        call assert_true('octree extrema bits', same_bits(current%pmin, reference%pmin) .and. &
            same_bits(current%pmax, reference%pmax))
        call assert_true('octree diagnostics', all(current%si == reference%si) .and. all(current%bi == reference%bi) &
            .and. all(current%depth == reference%depth) .and. all(current%starts == reference%starts) &
            .and. all(current%counts == reference%counts) .and. all(current%terminal .eqv. reference%terminal) &
            .and. all(current%oct_status == reference%oct_status) .and. all(current%oct_samples == reference%oct_samples) &
            .and. all(current%oct_leaves == reference%oct_leaves) .and. current%total_samples == reference%total_samples &
            .and. current%total_leaves == reference%total_leaves)
        call assert_true('probability and final trajectory bits', same_bits(reshape(current%direct, [7*36*n]), &
            reshape(reference%direct, [7*36*n])))
    end do
    call boundaries%destroy()
    deallocate(probabilities(1)%ref)
    print *, 'test_thread_determinism: all tests passed.'
contains
    logical function same_bits(a, b)
        real(c_double), intent(in) :: a(:), b(:)
        same_bits = all(transfer(a, [0_int64], size(a)) == transfer(b, [0_int64], size(b)))
    end function

    subroutine evaluate(threads)
        integer, intent(in) :: threads
        integer :: k
        call evaluate_dust_currents_impl(simulator, 1, n, points, bounds, 4, psi, &
            0.7d0, 32, 0, [3, 4, 3], 2, 4096, 128, 1d-4, 1d-4, 1.5d0, 0, &
            current%flux, current%error, current%status, current%samples, current%terminals, threads)
        call get_probabilities_octree_impl(simulator, 1, n, points, bounds, 0.7d0, 32, 0, [3, 4, 3], &
            2, 4096, 128, 1d-4, 1d-4, 1.5d0, 0, current%si, current%velocities, current%probabilities, &
            current%bi, current%boxes, current%pmin, current%pmax, current%depth, current%starts, current%counts, &
            current%oct_status, current%oct_samples, current%oct_leaves, current%total_samples, current%total_leaves, &
            threads, return_leaf_is_terminal=current%terminal)
!$      call omp_set_num_threads(threads)
        !$omp parallel do schedule(dynamic, 1)
        do k = 1, 36*n
            block
                type(t_ProbabilityRecord) :: record
                record = solver%calculate_probability(new_Particle(1d0, points(:, (k-1)/36+1), &
                    v(:, mod(k-1, 36)+1)), 0.7d0, 32, .false.)
                current%direct(:, k) = [record%probability, record%particle%position, record%particle%velocity]
            end block
        end do
        !$omp end parallel do
    end subroutine
end program

module m_emses_solver
    use, intrinsic :: iso_c_binding
    use omp_lib

    use m_vector, only: rot3d_y, rot3d_z
    use finbound, only: t_Boundary, t_BoundaryList, &
                        t_PlaneXYZ, &
                        new_PlaneX, new_PlaneY, new_PlaneZ
    use forbear, only: bar_object

    use m_vdsolverf_core
    use allcom, only: npbnd, phiz, phixy, &
                      nx, ny, nz, &
                      zssurf, &
                      peth, path, vdri, vdthz, vdthxy, &
                      qm, wc
    use m_namelist
    use emses_boundaries

    implicit none

    private
    public get_probabirities

contains

    subroutine get_backtrace(inppath, &
                             length, &
                             lx, ly, lz, &
                             ebvalues, &
                             ispec, &
                             position, &
                             velocity, &
                             dt, &
                             max_step, &
                             use_adaptive_dt, &
                             max_probabirity_types, &
                             return_ts, &
                             return_positions, &
                             return_velocities, &
                             return_last_step &
                             ) bind(c)
        character(1, c_char), intent(in) :: inppath(*)
        integer(c_int), value, intent(in) :: length
        integer(c_int), value, intent(in) :: lx
        integer(c_int), value, intent(in) :: ly
        integer(c_int), value, intent(in) :: lz
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
        integer(c_int), value, intent(in) :: ispec
        real(c_double), intent(in) :: position(3)
        real(c_double), intent(in) :: velocity(3)
        real(c_double), value, intent(in) :: dt
        integer(c_int), value, intent(in) :: max_step
        integer(c_int), value, intent(in) :: use_adaptive_dt
        integer(c_int), value, intent(in) :: max_probabirity_types
        real(c_double), intent(out) :: return_ts(max_step)
        real(c_double), intent(out) :: return_positions(3, max_step)
        real(c_double), intent(out) :: return_velocities(3, max_step)
        integer(c_int), intent(out) :: return_last_step

        type(t_ESSimulator) :: simulator

        simulator = create_simulator(inppath, length, &
                                     lx, ly, lz, &
                                     ebvalues, &
                                     ispec, &
                                     max_probabirity_types)

        block
            type(t_Particle) :: particle
            type(t_BacktraceRecord) :: record
            integer :: istep
            type(t_Particle) :: trace

            particle = new_Particle(qm(ispec), position(:), velocity(:))
            record = simulator%backtrace(particle, &
                                         dt, max_step, &
                                         use_adaptive_dt == 1)

            do istep = 1, record%last_step
                trace = record%traces(istep)

                return_ts(istep) = trace%t
                return_positions(:, istep) = trace%position(:)
                return_velocities(:, istep) = trace%velocity(:)
            end do

            return_last_step = record%last_step
        end block

        call simulator%boundaries%destroy
    end subroutine

    subroutine get_probabirities(inppath, &
                                 length, &
                                 lx, ly, lz, &
                                 ebvalues, &
                                 ispec, &
                                 npcls, &
                                 positions, &
                                 velocities, &
                                 dt, &
                                 max_step, &
                                 use_adaptive_dt, &
                                 max_probabirity_types, &
                                 return_probabirities, &
                                 return_positions, &
                                 return_velocities &
                                 ) bind(c)
        character(1, c_char), intent(in) :: inppath(*)
        integer(c_int), value, intent(in) :: length
        integer(c_int), value, intent(in) :: lx
        integer(c_int), value, intent(in) :: ly
        integer(c_int), value, intent(in) :: lz
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
        integer(c_int), value, intent(in) :: ispec
        integer(c_int), value, intent(in) :: npcls
        real(c_double), intent(in) :: positions(3, npcls)
        real(c_double), intent(in) :: velocities(3, npcls)
        real(c_double), value, intent(in) :: dt
        integer(c_int), value, intent(in) :: max_step
        integer(c_int), value, intent(in) :: use_adaptive_dt
        integer(c_int), value, intent(in) :: max_probabirity_types
        real(c_double), intent(out) :: return_probabirities(npcls)
        real(c_double), intent(out) :: return_positions(3, npcls)
        real(c_double), intent(out) :: return_velocities(3, npcls)

        type(t_ESSimulator) :: simulator
        type(bar_object) :: bar
        integer :: ipcl

        integer :: ithread, nthread

        simulator = create_simulator(inppath, length, &
                                     lx, ly, lz, &
                                     ebvalues, &
                                     ispec, &
                                     max_probabirity_types)

        call bar%initialize(filled_char_string='+', &
                            prefix_string='progress |', &
                            suffix_string='| ', &
                            add_progress_percent=.true.)
        call bar%start

        !$omp parallel private(ithread)
        nthread = omp_get_num_threads()
        ithread = omp_get_thread_num()

        ! Note: Reason for not using the omp do statement
        !     Particles are often passed sorted by position and velocity.
        !     Particles with similar phase values tend to follow similar trajectories.
        !     This results in similar computation times.
        !
        !     Using the typical omp do method to divide particles evenly can cause an unbalanced load on each thread.
        !     This imbalance means parallelization does not improve speedup.
        !
        !     Therefore, this process samples from the particle list based on the maximum number of threads.
        do ipcl = ithread + 1, npcls, nthread
            if (ithread == 0) then
                ! When you print a progress bar to 100%, the opening is printed.
                ! Therefore, it is modified to print 99% until the last particle is processed.
                call bar%update(current=min(0.99d0, dble(ipcl)/dble(npcls)))
            end if

            block
                type(t_ProbabirityRecord) :: record
                type(t_Particle) :: particle

                particle = new_Particle(qm(ispec), positions(:, ipcl), velocities(:, ipcl))
                record = simulator%calculate_probabirity(particle, dt, max_step, use_adaptive_dt == 1)

                if (record%is_valid) then
                    return_probabirities(ipcl) = record%probabirity
                    return_positions(:, ipcl) = record%particle%position(:)
                    return_velocities(:, ipcl) = record%particle%velocity(:)
                else
                    return_probabirities(ipcl) = -1.0d0
                end if
            end block
        end do
        !$omp end parallel

        call bar%update(current=1d0)

        call bar%destroy
        call simulator%boundaries%destroy
    end subroutine

    function create_simulator(inppath, length, lx, ly, lz, ebvalues, ispec, max_probabirity_types) result(simulator)
        character(1, c_char), intent(in) :: inppath(*)
        integer(c_int), value, intent(in) :: length
        integer(c_int), value, intent(in) :: lx
        integer(c_int), value, intent(in) :: ly
        integer(c_int), value, intent(in) :: lz
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
        integer(c_int), value, intent(in) :: ispec
        integer(c_int), value, intent(in) :: max_probabirity_types
        type(t_ESSimulator) :: simulator

        type(t_VectorFieldGrid), target :: eb

        type(t_BoundaryList) :: boundaries
        type(tp_Probabirity), allocatable :: probabirity_functions(:)
        integer :: n_probabirity_functions = 0

        allocate (probabirity_functions(max_probabirity_types))

        block
            character(length) :: s
            integer :: i
            do i = 1, length
                s(i:i) = inppath(i)
            end do
            call read_namelist(s)
        end block

        block
            integer :: isdoms(2, 3)
            integer :: boundary_conditions(3)

            allocate (probabirity_functions(n_probabirity_functions + 1)%ref, source=new_ZeroProbabirity())
            n_probabirity_functions = n_probabirity_functions + 1

            isdoms = reshape([[0, lx], [0, ly], [0, lz]], [2, 3])
            boundaries = create_simple_collision_boundaries(isdoms, tag=n_probabirity_functions)

            eb = new_VectorFieldGrid(6, lx, ly, lz, ebvalues(1:6, 1:lx + 1, 1:ly + 1, 1:lz + 1))

            call add_probabirity_boundaries(boundaries, ispec, n_probabirity_functions, probabirity_functions)

            boundary_conditions(:) = npbnd(:, ispec)
            simulator = new_ESSimulator(lx, ly, lz, &
                                        boundary_conditions, &
                                        eb, &
                                        boundaries, &
                                        probabirity_functions)
        end block
    end function

    subroutine add_probabirity_boundaries(boundaries, &
                                          ispec, &
                                          n_probabirity_functions, &
                                          probabirity_functions)
        type(t_BoundaryList), intent(inout) :: boundaries
        integer, intent(in) :: ispec
        integer, intent(inout) :: n_probabirity_functions
        type(tp_Probabirity), intent(inout) :: probabirity_functions(:)

        call add_external_boundaries

    contains

        subroutine add_external_boundaries
            double precision :: vmean(3), vthermal(3)
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZ), pointer :: pplane

            integer :: tag_zero
            integer :: tag_vdist

            vmean(:) = vdri_vector(ispec)
            vthermal(:) = vth_vector(ispec)

            allocate (probabirity_functions(n_probabirity_functions + 1)%ref, &
                      source=new_ZeroProbabirity())
            n_probabirity_functions = n_probabirity_functions + 1
            tag_zero = n_probabirity_functions

            allocate (probabirity_functions(n_probabirity_functions + 1)%ref, &
                      source=new_MaxwellianProbabirity(vmean, vthermal))
            n_probabirity_functions = n_probabirity_functions + 1
            tag_vdist = n_probabirity_functions

            if (npbnd(1, ispec) == 2) then ! X-Boundary
                ! X lower boundary
                allocate (pplane)
                pplane = new_PlaneX(0d0)
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)

                ! X higher boundary
                allocate (pplane)
                pplane = new_PlaneX(dble(nx))
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(2, ispec) == 2) then ! Y-Boundary
                ! Y lower boundary
                allocate (pplane)
                pplane = new_PlaneY(0d0)
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)

                ! Y higher boundary
                allocate (pplane)
                pplane = new_PlaneY(dble(ny))
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(3, ispec) == 2) then ! Z-Boundary
                ! Z lower boundary
                allocate (pplane)
                pplane = new_PlaneZ(0d0)
                pbound => pplane
                if (zssurf < 0d0) then
                    pbound%material%tag = tag_vdist
                else
                    pbound%material%tag = tag_zero
                end if
                call boundaries%add_boundary(pbound)

                ! Z higher boundary
                allocate (pplane)
                pplane = new_PlaneZ(dble(nz))
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)
            end if
        end subroutine

    end subroutine

    function vdri_vector(ispec) result(ret)
        integer, intent(in) :: ispec
        double precision :: ret(3)

        ret(:) = [0d0, 0d0, vdri(ispec)]
        ret(:) = rot3d_y(ret, vdthz(ispec))
        ret(:) = rot3d_z(ret, vdthxy(ispec))
    end function

    function vth_vector(ispec) result(ret)
        integer, intent(in) :: ispec
        double precision :: ret(3)

        ret(:) = [peth(ispec), peth(ispec), path(ispec)]
        ret(:) = rot3d_y(ret, phiz)
        ret(:) = rot3d_z(ret, phixy)
    end function

end module

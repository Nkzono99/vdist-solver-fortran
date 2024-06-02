module m_emses_solver
    !! This module provides functionalities for the EMSES simulation data,
    !! including backtrace and probability calculations.

    use, intrinsic :: iso_c_binding

    ! Use OpenMP library
!$  use omp_lib

    use m_vector, only: rot3d_y, rot3d_z
    use finbound, only: t_Boundary, t_BoundaryList, &
                        t_PlaneXYZ, &
                        new_PlaneX, new_PlaneY, new_PlaneZ, &
                        t_RectangleXYZ, &
                        new_RectangleX, new_RectangleY, new_RectangleZ
    use forbear, only: bar_object

    use m_vdsolverf_core
    use allcom
    use m_namelist
    use emses_boundaries

    implicit none

    private
    public get_probabilities

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
                             max_probability_types, &
                             return_ts, &
                             return_positions, &
                             return_velocities, &
                             return_last_step &
                             ) bind(c)
        !! Perform backtrace of a particle and return the trace data.

        character(1, c_char), intent(in) :: inppath(*)
            !! Path to the input file
        integer(c_int), value, intent(in) :: length
            !! Length of the input path
        integer(c_int), value, intent(in) :: lx
            !! Number of grid cells in the x direction
        integer(c_int), value, intent(in) :: ly
            !! Number of grid cells in the y direction
        integer(c_int), value, intent(in) :: lz
            !! Number of grid cells in the z direction
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
            !! Electric and magnetic field values
        integer(c_int), value, intent(in) :: ispec
            !! Species index
        real(c_double), intent(in) :: position(3)
            !! Initial position of the particle
        real(c_double), intent(in) :: velocity(3)
            !! Initial velocity of the particle
        real(c_double), value, intent(in) :: dt
            !! Time step width (Distance moving in one step (x += v/abs(v)*dt) when use_adaptive_dt is .true.)
        integer(c_int), value, intent(in) :: max_step
            !! Maximum number of steps
        integer(c_int), value, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        integer(c_int), value, intent(in) :: max_probability_types
            !! Maximum number of probability types
        real(c_double), intent(out) :: return_ts(max_step)
            !! Array to store time steps
        real(c_double), intent(out) :: return_positions(3, max_step)
            !! Array to store positions
        real(c_double), intent(out) :: return_velocities(3, max_step)
            !! Array to store velocities
        integer(c_int), intent(out) :: return_last_step
            !! Last step index

        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver

        simulator = create_simulator(inppath, length, &
                                     lx, ly, lz, &
                                     ebvalues, &
                                     ispec, &
                                     max_probability_types)
        solver = new_Solver(simulator)

        block
            type(t_Particle) :: particle
            type(t_BacktraceRecord) :: record
            integer :: istep
            type(t_Particle) :: trace

            particle = new_Particle(qm(ispec), position(:), velocity(:))
            record = solver%backtrace(particle, &
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

    subroutine get_probabilities(inppath, &
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
                                 max_probability_types, &
                                 return_probabilities, &
                                 return_positions, &
                                 return_velocities, &
                                 n_threads &
                                 ) bind(c)
        !! Calculate probabilities for multiple particles and return the results.

        character(1, c_char), intent(in) :: inppath(*)
            !! Path to the input file
        integer(c_int), value, intent(in) :: length
            !! Length of the input path
        integer(c_int), value, intent(in) :: lx
            !! Number of grid cells in the x direction
        integer(c_int), value, intent(in) :: ly
            !! Number of grid cells in the y direction
        integer(c_int), value, intent(in) :: lz
            !! Number of grid cells in the z direction
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
            !! Electric and magnetic field values
        integer(c_int), value, intent(in) :: ispec
            !! Species index
        integer(c_int), value, intent(in) :: npcls
            !! Number of particles
        real(c_double), intent(in) :: positions(3, npcls)
            !! Initial positions of the particles
        real(c_double), intent(in) :: velocities(3, npcls)
            !! Initial velocities of the particles
        real(c_double), value, intent(in) :: dt
            !! Time step
        integer(c_int), value, intent(in) :: max_step
            !! Maximum number of steps
        integer(c_int), value, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        integer(c_int), value, intent(in) :: max_probability_types
            !! Maximum number of probability types
        real(c_double), intent(out) :: return_probabilities(npcls)
            !! Array to store calculated probabilities
        real(c_double), intent(out) :: return_positions(3, npcls)
            !! Array to store final positions
        real(c_double), intent(out) :: return_velocities(3, npcls)
            !! Array to store final velocities
        integer(c_int), optional, intent(in) :: n_threads
            !! Number of threads for parallel computation

        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver

        type(bar_object) :: bar
        integer :: ipcl

!$      integer :: ithread

        simulator = create_simulator(inppath, length, &
                                     lx, ly, lz, &
                                     ebvalues, &
                                     ispec, &
                                     max_probability_types)
        solver = new_Solver(simulator)

        call bar%initialize(filled_char_string='+', &
                            prefix_string='progress |', &
                            suffix_string='| ', &
                            add_progress_percent=.true.)
        call bar%start

!$      if (present(n_threads)) then
!$          call omp_set_num_threads(n_threads)
!$      end if

        !$omp parallel

!$      ithread = omp_get_thread_num()

        ! Note: Reason for using the "omp do schedule(static, 1)" statement
        !     Particles are often passed sorted by position and velocity.
        !     Particles with similar phase values tend to follow similar trajectories.
        !     This results in similar computation times.
        !
        !     Using the typical omp do method to divide particles evenly can cause an unbalanced load on each thread.
        !     This imbalance means parallelization does not improve speedup.
        !
        !     Therefore, this process samples every chunk size (= 1) from the particle list.
        !$omp do schedule(static, 1)
        do ipcl = 1, npcls
!$          if (ithread == 0) then
                ! When you print a progress bar to 100%, the opening is printed.
                ! Therefore, it is modified to print 99% until the last particle is processed.
                call bar%update(current=min(0.99d0, dble(ipcl)/dble(npcls)))
!$          end if

            block
                type(t_ProbabilityRecord) :: record
                type(t_Particle) :: particle

                particle = new_Particle(qm(ispec), positions(:, ipcl), velocities(:, ipcl))
                record = solver%calculate_probability(particle, dt, max_step, use_adaptive_dt == 1)

                if (record%is_valid) then
                    return_probabilities(ipcl) = record%probability
                    return_positions(:, ipcl) = record%particle%position(:)
                    return_velocities(:, ipcl) = record%particle%velocity(:)
                else
                    return_probabilities(ipcl) = -1.0d0
                end if
            end block
        end do
        !$omp end do
        !$omp end parallel

        call bar%update(current=1d0)

        call bar%destroy
        call simulator%boundaries%destroy
    end subroutine

    function create_simulator(inppath, length, lx, ly, lz, ebvalues, ispec, max_probability_types) result(simulator)
        !! Create and initialize a new ES simulator object.

        character(1, c_char), intent(in) :: inppath(*)
            !! Path to the input file
        integer(c_int), value, intent(in) :: length
            !! Length of the input path
        integer(c_int), value, intent(in) :: lx
            !! Number of grid cells in the x direction
        integer(c_int), value, intent(in) :: ly
            !! Number of grid cells in the y direction
        integer(c_int), value, intent(in) :: lz
            !! Number of grid cells in the z direction
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
            !! Electric and magnetic field values
        integer(c_int), value, intent(in) :: ispec
            !! Species index
        integer(c_int), value, intent(in) :: max_probability_types
            !! Maximum number of probability types
        type(t_ESSimulator) :: simulator
            !! Simulator object

        type(t_VectorFieldGrid), target :: eb

        type(t_BoundaryList) :: boundaries

        type(tp_Probability), allocatable :: probability_functions(:)
        integer :: n_probability_functions = 0

        allocate (probability_functions(max_probability_types))

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

            allocate (probability_functions(n_probability_functions + 1)%ref, source=new_ZeroProbability())
            n_probability_functions = n_probability_functions + 1

            isdoms = reshape([[0, lx], [0, ly], [0, lz]], [2, 3])
            boundaries = create_simple_collision_boundaries(isdoms, tag=n_probability_functions)

            eb = new_VectorFieldGrid(6, lx, ly, lz, ebvalues(1:6, 1:lx + 1, 1:ly + 1, 1:lz + 1))

            call add_probability_boundaries(boundaries, ispec, n_probability_functions, probability_functions)

            boundary_conditions(:) = npbnd(:, ispec)
            simulator = new_ESSimulator(lx, ly, lz, &
                                        boundary_conditions, &
                                        eb, &
                                        boundaries, &
                                        probability_functions)
        end block
    end function

    subroutine add_probability_boundaries(boundaries, &
                                          ispec, &
                                          n_probability_functions, &
                                          probability_functions)
        !! Add probability boundaries to the simulator.

        type(t_BoundaryList), intent(inout) :: boundaries
            !! Boundary list
        integer, intent(in) :: ispec
            !! Species index
        integer, intent(inout) :: n_probability_functions
            !! Number of probability functions
        type(tp_Probability), intent(inout) :: probability_functions(:)
            !! Array of probability functions

        call add_external_boundaries
        call add_emission_surface(priority=1)

    contains

        subroutine add_external_boundaries
            !! Add external boundaries to the simulator.

            double precision :: vmean(3), vthermal(3)
                !! Mean velocity and thermal velocity
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZ), pointer :: pplane

            integer :: tag_zero
                !! Tag for zero probability boundary
            integer :: tag_vdist
                !! Tag for maxwellian velocity distribution boundary

            vmean(:) = vdri_vector(ispec)
            vthermal(:) = vth_vector(ispec)

            allocate (probability_functions(n_probability_functions + 1)%ref, &
                      source=new_ZeroProbability())
            n_probability_functions = n_probability_functions + 1
            tag_zero = n_probability_functions

            allocate (probability_functions(n_probability_functions + 1)%ref, &
                      source=new_MaxwellianProbability(vmean, vthermal))
            n_probability_functions = n_probability_functions + 1
            tag_vdist = n_probability_functions

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

        subroutine add_emission_surface(priority)
            integer, intent(in) :: priority

            integer :: iepl_start, iepl_end
            integer :: iepl
            integer :: tag_vdist

            if (nepl(ispec) == 0) then
                return
            end if

            if (ispec == 1) then
                iepl_start = 1
            else
                iepl_start = sum(nepl(1:ispec - 1)) + 1
            end if

            iepl_end = sum(nepl(1:ispec))

            do iepl = iepl_start, iepl_end
                block
                    double precision :: vmean(3)
                    double precision :: vthermal(3)

                    vmean = emission_vdri_vector(ispec, iepl)
                    vthermal = emission_vth_vector(ispec, iepl)

                    allocate (probability_functions(n_probability_functions + 1)%ref, &
                              source=new_MaxwellianProbability(vmean, vthermal))
                    n_probability_functions = n_probability_functions + 1
                    tag_vdist = n_probability_functions
                end block

                block
                    class(t_RectangleXYZ), pointer :: prect
                    class(t_Boundary), pointer :: pbound
                    double precision :: origin(3), wx, wy, wz

                    origin(:) = [xmine(iepl), ymine(iepl), zmine(iepl)]
                    wx = xmaxe(iepl) - xmine(iepl)
                    wy = ymaxe(iepl) - ymine(iepl)
                    wz = zmaxe(iepl) - zmine(iepl)

                    if (abs(nemd(iepl)) == 1) then
                        allocate (prect, source=new_RectangleX(origin, wy, wz))
                    else if (abs(nemd(iepl)) == 2) then
                        allocate (prect, source=new_RectangleY(origin, wy, wz))
                    else if (abs(nemd(iepl)) == 3) then
                        allocate (prect, source=new_RectangleZ(origin, wy, wz))
                    end if
                    prect%priority = priority

                    pbound => prect
                    pbound%material%tag = tag_vdist
                    call boundaries%add_boundary(pbound)
                end block
            end do
        end subroutine

    end subroutine

end module

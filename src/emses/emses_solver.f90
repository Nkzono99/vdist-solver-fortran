module m_emses_solver
    !! C-facing entry points for EMSES backtrace and probability calculations.
    !!
    !! Each public subroutine is exposed through `bind(c)` so Python's ctypes
    !! wrapper (`vdsolverf/emses/wrapper.py`) can call it directly. Simulator
    !! construction is delegated to `m_emses_simulator_builder` to keep this
    !! module focused on the C bridge.

    use, intrinsic :: iso_c_binding

    ! Use OpenMP library
!$  use omp_lib

    use forbear, only: bar_object

    use m_vdsolverf_core
    use m_allcom, only: qm
    use m_emses_autorange, only: estimate_velocity_range_map_impl
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator

    implicit none

    private
    public get_probabilities
    public get_backtraces
    public estimate_velocity_range_map

contains

    subroutine estimate_velocity_range_map( &
        inppath, &
        length, &
        lx, ly, lz, &
        ebvalues, &
        ispec, &
        dt, &
        coverage_sigma, &
        safety_factor, &
        max_step, &
        use_adaptive_dt, &
        max_probability_types, &
        source_samples_per_cell, &
        velocity_sample_mode, &
        minimum_count, &
        collect_moments, &
        show_progress, &
        return_vx_min, &
        return_vx_max, &
        return_vy_min, &
        return_vy_max, &
        return_vz_min, &
        return_vz_max, &
        return_count, &
        return_weight_sum, &
        return_mean_v, &
        return_cov_v, &
        return_status, &
        return_confidence, &
        n_threads &
        ) bind(c)
        !! Estimate per-cell velocity ranges by forward-tracing source envelopes.

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
        real(c_double), intent(in) :: ebvalues(9, lx + 1, ly + 1, lz + 1)
            !! Relocated E/B values plus staggered accumulated-charge E values
        integer(c_int), value, intent(in) :: ispec
            !! Species index
        real(c_double), value, intent(in) :: dt
            !! Forward trace step width
        real(c_double), value, intent(in) :: coverage_sigma
            !! Source support radius in units of the diagonal thermal scales
        real(c_double), value, intent(in) :: safety_factor
            !! Multiplicative expansion applied to deposited min/max ranges
        integer(c_int), value, intent(in) :: max_step
            !! Maximum number of forward steps
        integer(c_int), value, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        integer(c_int), value, intent(in) :: max_probability_types
            !! Maximum number of probability types
        integer(c_int), value, intent(in) :: source_samples_per_cell
            !! Number of samples per tangential source cell direction
        integer(c_int), value, intent(in) :: velocity_sample_mode
            !! Reserved velocity support selector. 0 is ellipsoid support.
        integer(c_int), value, intent(in) :: minimum_count
            !! Count threshold used for LOW_COUNT diagnostics
        integer(c_int), value, intent(in) :: collect_moments
            !! Flag to collect mean/cov velocity diagnostics
        integer(c_int), value, intent(in) :: show_progress
            !! Flag to show progress bar logging
        real(c_double), intent(out) :: return_vx_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vx_max(lx, ly, lz)
        real(c_double), intent(out) :: return_vy_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vy_max(lx, ly, lz)
        real(c_double), intent(out) :: return_vz_min(lx, ly, lz)
        real(c_double), intent(out) :: return_vz_max(lx, ly, lz)
        integer(c_int), intent(out) :: return_count(lx, ly, lz)
        real(c_double), intent(out) :: return_weight_sum(lx, ly, lz)
        real(c_double), intent(out) :: return_mean_v(3, lx, ly, lz)
        real(c_double), intent(out) :: return_cov_v(3, 3, lx, ly, lz)
        integer(c_int), intent(out) :: return_status(lx, ly, lz)
        real(c_double), intent(out) :: return_confidence(lx, ly, lz)
        integer(c_int), optional, intent(in) :: n_threads
            !! Number of OpenMP threads for source-envelope tracing

        type(t_ESSimulator) :: simulator

        simulator = create_simulator(inppath, length, &
                                     lx, ly, lz, &
                                     ebvalues, &
                                     ispec, &
                                     max_probability_types)

        call estimate_velocity_range_map_impl( &
            simulator, &
            ispec, &
            lx, ly, lz, &
            dt, &
            coverage_sigma, &
            safety_factor, &
            max_step, &
            use_adaptive_dt, &
            source_samples_per_cell, &
            velocity_sample_mode, &
            minimum_count, &
            collect_moments, &
            show_progress, &
            return_vx_min, &
            return_vx_max, &
            return_vy_min, &
            return_vy_max, &
            return_vz_min, &
            return_vz_max, &
            return_count, &
            return_weight_sum, &
            return_mean_v, &
            return_cov_v, &
            return_status, &
            return_confidence, &
            n_threads)

        call destroy_simulator(simulator)
    end subroutine

    subroutine get_backtraces( &
        inppath, &
        length, &
        lx, ly, lz, &
        ebvalues, &
        ispec, &
        npcls, &
        positions, &
        velocities, &
        dt, &
        max_step, &
        output_interval, &
        use_adaptive_dt, &
        max_probability_types, &
        return_ts, &
        return_probabilities, &
        return_positions, &
        return_velocities, &
        return_last_steps, &
        n_threads &
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
        real(c_double), intent(in) :: ebvalues(9, lx + 1, ly + 1, lz + 1)
            !! Relocated E/B values plus staggered accumulated-charge E values
        integer(c_int), value, intent(in) :: ispec
            !! Species index
        integer(c_int), value, intent(in) :: npcls
            !! Number of particles
        real(c_double), intent(in) :: positions(3, npcls)
            !! Initial position of the particle
        real(c_double), intent(in) :: velocities(3, npcls)
            !! Initial velocity of the particle
        real(c_double), value, intent(in) :: dt
            !! Time step width (Distance moving in one step (x += v/abs(v)*dt) when use_adaptive_dt is .true.)
        integer(c_int), value, intent(in) :: max_step
            !! Maximum number of steps
        integer(c_int), value, intent(in) :: output_interval
            !! Output Interval
        integer(c_int), value, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        integer(c_int), value, intent(in) :: max_probability_types
            !! Maximum number of probability types
        real(c_double), intent(out) :: return_ts((max_step - 1)/output_interval + 2, npcls)
            !! Array to store time steps
        real(c_double), intent(out) :: return_probabilities(npcls)
        real(c_double), intent(out) :: return_positions(3, (max_step - 1)/output_interval + 2, npcls)
            !! Array to store positions
        real(c_double), intent(out) :: return_velocities(3, (max_step - 1)/output_interval + 2, npcls)
            !! Array to store velocities
        integer(c_int), intent(out) :: return_last_steps(npcls)
            !! Last step index
        integer(c_int), optional, intent(in) :: n_threads
        !! Number of threads for parallel computation

        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver

        type(bar_object) :: bar
        integer :: ipcl
        integer :: max_output_steps

        max_output_steps = (max_step - 1)/output_interval + 2

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

        ! Note: Reason for using the "omp do schedule(static, 1)" statement
        !     Particles are often passed sorted by position and velocity.
        !     Particles with similar phase values tend to follow similar trajectories.
        !     This results in similar computation times.
        !
        !     Using the typical omp do method to divide particles evenly can cause an unbalanced load on each thread.
        !     This imbalance means parallelization does not improve speedup.
        !
        !     Therefore, this process samples every chunk size (= 1) from the particle list.
        !$omp parallel do schedule(dynamic, 1)
        do ipcl = 1, npcls
!$          if (omp_get_thread_num() == 0) then
                ! When you print a progress bar to 100%, the opening is printed.
                ! Therefore, it is modified to print 99% until the last particle is processed.
                call bar%update(current=min(0.99d0, dble(ipcl)/dble(npcls)))
!$          end if

            block
                type(t_BacktraceRecord) :: record
                type(t_Particle) :: particle

                type(t_Particle) :: trace
                integer :: istep

                particle = new_Particle(qm(ispec), positions(:, ipcl), velocities(:, ipcl))
                record = solver%backtrace(particle, &
                                          dt, &
                                          max_step, &
                                          output_interval, &
                                          use_adaptive_dt == 1)

                do istep = 1, record%last_step
                    trace = record%traces(istep)

                    return_ts(istep, ipcl) = trace%t
                    return_positions(:, istep, ipcl) = trace%position(:)
                    return_velocities(:, istep, ipcl) = trace%velocity(:)
                end do

                istep = max_output_steps
                trace = record%traces(istep)
                return_ts(istep, ipcl) = trace%t
                return_positions(:, istep, ipcl) = trace%position(:)
                return_velocities(:, istep, ipcl) = trace%velocity(:)

                return_last_steps(ipcl) = record%last_step
                return_probabilities(ipcl) = record%probability
            end block
        end do
        !$omp end parallel do

        call bar%update(current=1d0)
        call bar%destroy
        call destroy_simulator(simulator)
    end subroutine

    subroutine get_probabilities( &
        inppath, &
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
        real(c_double), intent(in) :: ebvalues(9, lx + 1, ly + 1, lz + 1)
            !! Relocated E/B values plus staggered accumulated-charge E values
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

        !$omp parallel do schedule(dynamic, 1)
        do ipcl = 1, npcls
!$          if (omp_get_thread_num() == 0) then
                call bar%update(current=min(0.99d0, dble(ipcl)/dble(npcls)))
!$          end if

            block
                type(t_ProbabilityRecord) :: record
                type(t_Particle) :: particle

                particle = new_Particle(qm(ispec), positions(:, ipcl), velocities(:, ipcl))
                record = solver%calculate_probability(particle, &
                                                      dt, &
                                                      max_step, &
                                                      use_adaptive_dt == 1)

                if (record%is_valid) then
                    return_probabilities(ipcl) = record%probability
                    return_positions(:, ipcl) = record%particle%position(:)
                    return_velocities(:, ipcl) = record%particle%velocity(:)
                else
                    return_probabilities(ipcl) = -1.0d0
                    return_positions(:, ipcl) = record%particle%position(:)
                    return_velocities(:, ipcl) = record%particle%velocity(:)
                end if
            end block
        end do
        !$omp end parallel do

        call bar%update(current=1d0)

        call bar%destroy
        call destroy_simulator(simulator)
    end subroutine

end module

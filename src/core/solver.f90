module m_solver
    use finbound

    use m_particle
    use m_simulator
    use m_Probabilities
    use m_dust_charge_simulator

    implicit none

    private
    public t_Solver
    public new_Solver
    public t_BacktraceRecord
    public t_ProbabilityRecord
    public t_DustBacktraceRecord

    type t_BacktraceRecord
        !! Record for storing backtrace information
        type(t_Particle), allocatable :: traces(:)
            !! Array of particle traces
        double precision, allocatable :: ts(:)
            !! Array of time steps
        integer :: last_step
            !! Last step in the backtrace
    end type

    type t_DustBacktraceRecord
        !! Record for storing backtrace information
        type(t_DustParticle), allocatable :: traces(:)
            !! Array of particle traces
        double precision, allocatable :: ts(:)
            !! Array of time steps
        integer :: last_step
            !! Last step in the backtrace
    end type

    type t_ProbabilityRecord
        !! Record for storing probability calculation results
        logical :: is_valid = .false.
            !! Flag indicating if the record is valid (If .false., probability calculation is not completed within max_step)
        double precision :: t
            !! Time at which the probability is calculated since the start of the simulation
        double precision :: probability
            !! Calculated probability
        type(t_Particle) :: particle
            !! Particle state at the time of calculation probability
    end type

    type t_Solver

        type(t_ESSimulator) :: simulator
        type(t_DustChargeSimulator) :: charging_simulator

    contains

        procedure :: backtrace => solver_backtrace
        procedure :: calculate_probability => solver_calculate_probability
        procedure :: backtrace_dust => solver_backtrace_dust

    end type

contains

    function new_Solver(simulator, charging_simulator) result(obj)
        type(t_ESSimulator), intent(in) :: simulator
        type(t_DustChargeSimulator), intent(in), optional :: charging_simulator
        type(t_Solver) :: obj

        obj%simulator = simulator

        if (present(charging_simulator)) then
            obj%charging_simulator = charging_simulator
        end if
    end function

    function solver_backtrace(self, particle, dt, max_step, use_adaptive_dt) result(ret)
        !! Perform a backtrace of a particle.

        class(t_Solver), intent(in) :: self
            !! Instance of the ES simulator
        class(t_Particle), intent(in) :: particle
            !! Particle to backtrace
        double precision, intent(in) :: dt
            !! Time step width for backtracing
        integer, intent(in) :: max_step
            !! Maximum number of steps
        logical, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        type(t_BacktraceRecord) :: ret
            !! Backtrace record

        integer :: istep
        type(t_Particle) :: pcl

        pcl = particle
        call self%simulator%apply_boundary_condition(pcl)

        allocate (ret%traces(max_step))
        allocate (ret%ts(max_step))
        ret%traces(1) = pcl

        do istep = 2, max_step
            block
                type(t_CollisionRecord) :: record

                double precision :: tmp_dt

                if (use_adaptive_dt) then
                    tmp_dt = dt/sqrt(sum(pcl%velocity*pcl%velocity))
                else
                    tmp_dt = dt
                end if

                pcl = self%simulator%update(pcl, tmp_dt, record)

                ret%traces(istep) = pcl
                ret%ts(istep) = pcl%t

                if (record%is_collided) then
                    ret%last_step = istep
                    return
                end if
            end block
        end do

        ret%last_step = max_step
    end function

    function solver_calculate_probability(self, particle, dt, max_step, use_adaptive_dt) result(ret)
        !! Calculate the probability for a particle.

        class(t_Solver), intent(in) :: self
            !! Instance of the ES simulator
        class(t_Particle), intent(in) :: particle
            !! Particle for which to calculate the probability
        double precision, intent(in) :: dt
            !! Time step for the calculation
        integer, intent(in) :: max_step
            !! Maximum number of steps
        logical, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        type(t_ProbabilityRecord) :: ret
            !! Probability record

        integer :: i
        type(t_Particle) :: pcl

        pcl = particle
        call self%simulator%apply_boundary_condition(pcl)

        do i = 1, max_step
            block
                type(t_CollisionRecord) :: record
                type(tp_Probability), allocatable :: probability_function

                double precision :: tmp_dt

                if (use_adaptive_dt) then
                    tmp_dt = dt/sqrt(sum(pcl%velocity*pcl%velocity))
                else
                    tmp_dt = dt
                end if

                pcl = self%simulator%update(pcl, tmp_dt, record)

                if (record%is_collided) then
                    ret%is_valid = .true.
                    ret%particle = pcl

                    probability_function = self%simulator%probability_functions(record%material%tag)
                    ret%probability = probability_function%at(pcl%position, pcl%velocity)
                    return
                end if
            end block
        end do

        ret%is_valid = .false.
    end function

    function solver_backtrace_dust(self, dust, dt, max_step, use_adaptive_dt) result(ret)
        !! Perform a backtrace of a particle.

        class(t_Solver), intent(in) :: self
            !! Instance of the ES simulator
        class(t_DustParticle), intent(in) :: dust
            !! Particle to backtrace
        double precision, intent(in) :: dt
            !! Time step width for backtracing
        integer, intent(in) :: max_step
            !! Maximum number of steps
        logical, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        type(t_DustBacktraceRecord) :: ret
            !! Backtrace record

        integer :: istep
        type(t_DustParticle) :: dust_tmp

        dust_tmp = dust
        call self%simulator%apply_boundary_condition(dust_tmp%particle)

        allocate (ret%traces(max_step))
        allocate (ret%ts(max_step))
        ret%traces(1) = dust_tmp

        do istep = 2, max_step
            block
                type(t_CollisionRecord) :: record

                double precision :: tmp_dt
                type(t_Particle) :: pcl

                pcl = dust_tmp%particle

                if (use_adaptive_dt) then
                    tmp_dt = dt/sqrt(sum(pcl%velocity*pcl%velocity))
                else
                    tmp_dt = dt
                end if

                pcl = self%simulator%update(pcl, tmp_dt, record)

                dust_tmp = self%charging_simulator%update(dust_tmp, tmp_dt)
                dust_tmp%particle = pcl

                ret%traces(istep) = dust_tmp
                ret%ts(istep) = pcl%t

                if (record%is_collided) then
                    ret%last_step = istep
                    return
                end if
            end block
        end do

        ret%last_step = max_step
    end function

end module

module m_solver
    use finbound

    use m_particle
    use m_simulator
    use m_Probabilities

    implicit none

    private
    public t_Solver
    public new_Solver

    type t_Solver

        type(t_ESSimulator) :: simulator

    contains

        procedure :: backtrace => solver_backtrace
        procedure :: calculate_probability => solver_calculate_probability

    end type

contains

    function new_Solver(simulator) result(obj)
        type(t_ESSimulator), intent(in) :: simulator
        type(t_Solver) :: obj

        obj%simulator = simulator
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
end module

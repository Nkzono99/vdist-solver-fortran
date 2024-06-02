module m_simulator
    !! This module defines the ES simulator which simulates particle dynamics
    !! in a electrostatic field with specified boundary conditions.

    use m_vector, only: cross
    use finbound, only: t_BoundaryList, t_CollisionRecord
    use m_particle
    use m_Probabilities
    use m_field

    implicit none

    private
    public t_ESSimulator
    public new_ESSimulator
    public t_BacktraceRecord
    public t_ProbabilityRecord

    type t_BacktraceRecord
        !! Record for storing backtrace information
        type(t_Particle), allocatable :: traces(:)
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

    type t_ESSimulator
        !! Electric static simulator for particle dynamics
        integer :: nx
            !! Number of grid cells in the x direction
        integer :: ny
            !! Number of grid cells in the y direction
        integer :: nz
            !! Number of grid cells in the z direction
        integer :: boundary_conditions(3)
            !! Boundary conditions for each axis (0: periodic, 1: reflective)
        class(t_VectorField), allocatable :: eb
            !! Electric and magnetic fields
        type(t_BoundaryList) :: boundaries
            !! List of boundary conditions
        type(tp_Probability), allocatable :: probability_functions(:)
            !! Array of probability functions
    contains
        procedure :: update => esSimulator_update
            !! Update the state of a particle
        procedure :: backward => esSimulator_backward
            !! Perform a backward step for a particle
        procedure :: apply_boundary_condition => esSimulator_apply_boundary_condition
            !! Apply boundary conditions to a particle
    end type

contains

    function new_ESSimulator(nx, ny, nz, &
                             boundary_conditions, &
                             eb, &
                             boundaries, &
                             probability_functions) result(obj)
        !! Create a new ES simulator object with specified properties.

        integer, intent(in) :: nx
            !! Number of grid points in the x direction
        integer, intent(in) :: ny
            !! Number of grid points in the y direction
        integer, intent(in) :: nz
            !! Number of grid points in the z direction
        integer, intent(in) :: boundary_conditions(3)
            !! Boundary conditions for each axis (0: periodic, 1: reflective, 2: free)
        class(t_VectorField), intent(in) :: eb
            !! Electric and magnetic fields
        type(t_BoundaryList), intent(in) :: boundaries
            !! List of boundary conditions
        type(tp_Probability), intent(in) :: probability_functions(:)
            !! Array of probability functions

        type(t_ESSimulator) :: obj
            !! A new ES simulator object with specified properties

        obj%nx = nx
        obj%ny = ny
        obj%nz = nz
        obj%boundary_conditions(:) = boundary_conditions(:)
        obj%eb = eb
        obj%boundaries = boundaries

        allocate (obj%probability_functions(size(probability_functions)))

        block
            integer :: i
            do i = 1, size(probability_functions)
                obj%probability_functions(i)%ref => probability_functions(i)%ref
            end do
        end block
    end function

    function esSimulator_update(self, pcl, dt, use_adaptive_dt, record) result(pcl_new)
        !! Update the state of a particle.

        class(t_ESSimulator), intent(in) :: self
            !! Instance of the ES simulator
        type(t_Particle), intent(in) :: pcl
            !! Particle to update
        double precision, intent(in) :: dt
            !! Time step for the update
        logical, intent(in) :: use_adaptive_dt
            !! Flag to use adaptive time step
        type(t_CollisionRecord), intent(out) :: record
            !! Collision record
        type(t_Particle) :: pcl_new
            !! Updated particle

        double precision :: tmp_dt

        if (use_adaptive_dt) then
            tmp_dt = dt/sqrt(sum(pcl%velocity*pcl%velocity))
        else
            tmp_dt = dt
        end if

        pcl_new = self%backward(pcl, tmp_dt)

        record = self%boundaries%check_collision(pcl%position, pcl_new%position)

        if (record%is_collided) then
            block
                double precision :: r

                r = record%t
                pcl_new%q_m = pcl%q_m
                pcl_new%position = pcl%position*(1d0 - r) + pcl_new%position*r
                pcl_new%velocity = pcl%velocity*(1d0 - r) + pcl_new%velocity*r
                pcl_new%t = pcl%t + r
            end block

            return
        end if

        call self%apply_boundary_condition(pcl_new)
    end function

    function esSimulator_backward(self, particle, dt) result(ret)
        !! Perform a backward step for a particle.

        class(t_ESSimulator), intent(in) :: self
            !! Instance of the ES simulator
        class(t_Particle), intent(in) :: particle
            !! Particle to perform the backward step
        double precision, intent(in) :: dt
            !! Time step for the backward step
        type(t_Particle) :: ret
            !! Particle after the backward step

        double precision :: position_new(3)
        double precision :: velocity_new(3)

        position_new = particle%position - dt*particle%velocity

        block
            double precision :: dt2
            double precision :: eb(6)
            double precision :: ef(3)
            double precision :: bf(3)
            double precision :: s(3)
            double precision :: t(3)
            double precision :: upm(3)
            double precision :: upa(3)
            double precision :: upp(3)

            dt2 = -0.5d0*dt

            eb = self%eb%at(particle%position)
            ef(:) = eb(1:3)
            bf(:) = eb(4:6)

            t = bf*particle%q_m*dt2
            s = 2*t/(1 + t*t)

            upm = particle%velocity + particle%q_m*ef*dt2

            upa = upm + cross(upm, t)
            upp = upm + cross(upa, s)

            velocity_new = upp + particle%q_m*ef*dt2
        end block

        ret%t = particle%t + dt
        ret%q_m = particle%q_m
        ret%position = position_new
        ret%velocity = velocity_new
    end function

    subroutine esSimulator_apply_boundary_condition(self, particle)
        !! Apply boundary conditions to a particle.

        class(t_ESSimulator), intent(in) :: self
            !! Instance of the ES simulator
        class(t_Particle), intent(inout) :: particle
            !! Particle to apply boundary conditions to

        integer :: i
        double precision :: nxyz(3)

        nxyz(1) = self%nx
        nxyz(2) = self%ny
        nxyz(3) = self%nz

        do i = 1, 3
            if (self%boundary_conditions(i) == 0) then
                ! Periodic boundary condition
                particle%position(i) = modulo(particle%position(i), nxyz(i))
            else if (self%boundary_conditions(i) == 1) then
                ! Reflective boundary condition
                particle%position(i) = modulo(2*nxyz(i) - particle%position(i), nxyz(i))
            end if
        end do
    end subroutine

end module

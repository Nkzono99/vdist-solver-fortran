module m_simulator
    use m_field
    use m_particle
    use m_vector
    use m_probabirity_record
    use m_Probabirities
    use finbound
    implicit none

    type t_ESSimulator
        integer :: nx
        integer :: ny
        integer :: nz
        double precision :: q_m
        integer :: boundary_conditions(3)
        class(t_VectorField), pointer :: ef
        class(t_VectorField), pointer :: bf
        type(t_BoundaryList) :: boundaries
        type(tp_Probabirity), allocatable :: probabirity_functions(:)

    contains
        procedure :: calculate_probabirity => esSimulator_calculate_probabirity
        procedure :: backward => esSimulator_backward
        procedure :: apply_boundary_condition => esSimulator_apply_boundary_condition
    end type

contains

    function new_ESSimulator(nx, ny, nz, &
                             q_m, &
                             boundary_conditions, &
                             ef, bf, &
                             boundaries, &
                             probabirity_functions) result(obj)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        double precision, intent(in) :: q_m
        integer, intent(in) :: boundary_conditions(3)
        class(t_VectorField), pointer, intent(in) :: ef
        class(t_VectorField), pointer, intent(in) :: bf
        type(t_BoundaryList), intent(in) :: boundaries
        type(tp_Probabirity), intent(in) :: probabirity_functions(:)

        type(t_ESSimulator) :: obj

        obj%nx = nx
        obj%ny = ny
        obj%nz = nz
        obj%q_m = q_m
        obj%boundary_conditions(:) = boundary_conditions(:)
        obj%ef => ef
        obj%bf => bf
        obj%boundaries = boundaries
        obj%probabirity_functions(:) = probabirity_functions(:)
    end function

    function esSimulator_calculate_probabirity(self, particle, dt, max_step) result(ret)
        class(t_ESSimulator), intent(in) :: self
        class(t_Particle), intent(in) :: particle
        double precision, intent(in) :: dt
        integer, intent(in) :: max_step
        type(t_ProbabirityRecord) :: ret

        integer :: i
        type(t_Particle) :: pcl

        pcl = particle
        call self%apply_boundary_condition(pcl)

        do i = 1, max_step
            block
                type(t_Particle) :: pcl_new
                type(t_CollisionRecord) :: record
                double precision :: r
                type(tp_Probabirity), allocatable :: probabirity_function

                pcl_new = self%backward(pcl, dt)

                ! Detect collision (pcl, pcl_new)
                record = self%boundaries%check_collision(pcl%position, pcl_new%position)

                ! if detected => return prob
                if (record%is_collided) then
                    r = (record%t - particle%t)/dt
                    pcl_new%position = particle%position*(1d0 - r) + pcl_new%position*r
                    pcl_new%velocity = particle%velocity*(1d0 - r) + pcl_new%velocity*r
                    pcl_new%t = record%t

                    ret%is_valid = .true.
                    ret%particle = pcl_new

                    probabirity_function = self%probabirity_functions(record%material%tag)
                    ret%probabirity = probabirity_function%at(pcl_new%position, pcl_new%velocity)
                    return
                end if

                pcl = pcl_new
            end block

            ! Apply boundary condition
            call self%apply_boundary_condition(pcl)
        end do

        ret%is_valid = .false.
    end function

    subroutine esSimulator_apply_boundary_condition(self, particle)
        class(t_ESSimulator), intent(in) :: self
        class(t_Particle), intent(inout) :: particle

        integer :: i
        double precision :: nxyz(3)

        nxyz(1) = self%nx
        nxyz(2) = self%ny
        nxyz(3) = self%nz

        do i = 1, 3
            if (self%boundary_conditions(i) == 0) then ! if periodic
                particle%position(i) = modulo(particle%position(i), nxyz(i))
            else if (self%boundary_conditions(i) == 1) then ! if refrect
                particle%position(i) = modulo(2*nxyz(i) - particle%position(i), nxyz(i))
            end if
        end do
    end subroutine

    function esSimulator_backward(self, particle, dt) result(ret)
        class(t_ESSimulator), intent(in) :: self
        class(t_Particle), intent(in) :: particle
        double precision, intent(in) :: dt
        type(t_Particle) :: ret

        double precision :: position_new(3)
        double precision :: velocity_new(3)

        position_new = particle%position - dt*particle%velocity

        block
            double precision :: dt2
            double precision :: ef(3)
            double precision :: bf(3)
            double precision :: s(3)
            double precision :: t(3)
            double precision :: upm(3)
            double precision :: upa(3)
            double precision :: upp(3)
            dt2 = 0.5d0*dt

            ef = self%ef%at(particle%position)
            bf = self%bf%at(particle%position)

            t = bf*self%q_m*dt2
            s = 2*t/(1 + t*t)

            upm = particle%velocity + self%q_m*ef*dt2

            upa = upm + cross(upm, t)
            upp = upm + cross(upa, s)

            velocity_new = upp + self%q_m*ef*dt2
        end block

        ret%t = particle%t + dt
        ret%position = position_new
        ret%velocity = velocity_new
    end function
end module

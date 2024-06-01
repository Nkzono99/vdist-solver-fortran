module m_simulator
    use m_vector, only: cross
    use finbound, only: t_BoundaryList, t_CollisionRecord

    use m_particle
    use m_Probabirities
    use m_field

    implicit none

    private
    public t_ESSimulator
    public new_ESSimulator
    public t_BacktraceRecord
    public t_ProbabirityRecord

    type t_BacktraceRecord
        type(t_Particle), allocatable :: traces(:)
        double precision, allocatable :: ts(:)
        integer :: last_step
    end type

    type t_ProbabirityRecord
        logical :: is_valid = .false.
        double precision :: t
        double precision :: probabirity
        type(t_Particle) :: particle
    end type

    type t_ESSimulator
        integer :: nx
        integer :: ny
        integer :: nz
        integer :: boundary_conditions(3)
        class(t_VectorField), allocatable :: eb
        type(t_BoundaryList) :: boundaries
        type(tp_Probabirity), allocatable :: probabirity_functions(:)
    contains
        procedure :: backtrace => esSimulator_backtrace
        procedure :: calculate_probabirity => esSimulator_calculate_probabirity
        procedure :: backward => esSimulator_backward
        procedure :: update => esSimulator_update
        procedure :: apply_boundary_condition => esSimulator_apply_boundary_condition
    end type

contains

    function new_ESSimulator(nx, ny, nz, &
                             boundary_conditions, &
                             eb, &
                             boundaries, &
                             probabirity_functions) result(obj)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        integer, intent(in) :: boundary_conditions(3)
        class(t_VectorField), intent(in) :: eb
        type(t_BoundaryList), intent(in) :: boundaries
        type(tp_Probabirity), intent(in) :: probabirity_functions(:)

        type(t_ESSimulator) :: obj

        obj%nx = nx
        obj%ny = ny
        obj%nz = nz
        obj%boundary_conditions(:) = boundary_conditions(:)
        obj%eb = eb
        obj%boundaries = boundaries

        allocate (obj%probabirity_functions(size(probabirity_functions)))

        block
            integer :: i
            do i = 1, size(probabirity_functions)
                obj%probabirity_functions(i)%ref => probabirity_functions(i)%ref
            end do
        end block
    end function

    function esSimulator_backtrace(self, particle, dt, max_step, use_adaptive_dt) result(ret)
        class(t_ESSimulator), intent(in) :: self
        class(t_Particle), intent(in) :: particle
        double precision, intent(in) :: dt
        integer, intent(in) :: max_step
        logical, intent(in) :: use_adaptive_dt
        type(t_BacktraceRecord) :: ret

        integer :: istep
        type(t_Particle) :: pcl

        pcl = particle
        call self%apply_boundary_condition(pcl)

        allocate (ret%traces(max_step))
        allocate (ret%ts(max_step))
        ret%traces(1) = pcl

        do istep = 2, max_step
            block
                type(t_CollisionRecord) :: record
                pcl = self%update(pcl, dt, use_adaptive_dt, record)

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

    function esSimulator_calculate_probabirity(self, particle, dt, max_step, use_adaptive_dt) result(ret)
        class(t_ESSimulator), intent(in) :: self
        class(t_Particle), intent(in) :: particle
        double precision, intent(in) :: dt
        integer, intent(in) :: max_step
        logical, intent(in) :: use_adaptive_dt
        type(t_ProbabirityRecord) :: ret

        integer :: i
        type(t_Particle) :: pcl

        pcl = particle
        call self%apply_boundary_condition(pcl)

        do i = 1, max_step
            block
                type(t_CollisionRecord) :: record
                type(tp_Probabirity), allocatable :: probabirity_function
                pcl = self%update(pcl, dt, use_adaptive_dt, record)

                if (record%is_collided) then
                    ret%is_valid = .true.
                    ret%particle = pcl

                    probabirity_function = self%probabirity_functions(record%material%tag)
                    ret%probabirity = probabirity_function%at(pcl%position, pcl%velocity)
                    return
                end if
            end block
        end do

        ret%is_valid = .false.
    end function

    function esSimulator_update(self, pcl, dt, use_adaptive_dt, record) result(pcl_new)
        class(t_ESSimulator), intent(in) :: self
        type(t_Particle), intent(in) :: pcl
        double precision, intent(in) :: dt
        logical, intent(in) :: use_adaptive_dt
        type(t_CollisionRecord), intent(out) :: record
        type(t_Particle) :: pcl_new

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
                pcl_new%position = pcl%position*(1d0 - r) + pcl_new%position*r
                pcl_new%velocity = pcl%velocity*(1d0 - r) + pcl_new%velocity*r
                pcl_new%t = pcl%t + r
            end block

            return
        end if

        call self%apply_boundary_condition(pcl_new)
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
        ret%position = position_new
        ret%velocity = velocity_new
    end function

end module

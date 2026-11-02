module m_simulator
    !! This module defines the ES simulator which simulates particle dynamics
    !! in a electrostatic field with specified boundary conditions.

    use m_vector, only: cross
    use finbound, only: t_BoundaryList, t_CollisionRecord
    use m_particle
    use m_probabilities
    use m_field

    implicit none

    private
    public t_ESSimulator
    public new_ESSimulator

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
        procedure :: adaptive_dt => esSimulator_adaptive_dt
            !! Time step that advances a particle by about |dt| grid units
        procedure :: check_collision => esSimulator_check_collision
            !! Collision test that follows the segment through periodic faces
        procedure, private :: domain_length => esSimulator_domain_length
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

    function esSimulator_update(self, pcl, dt, record, unwrapped_position) result(pcl_new)
        !! Update the state of a particle.

        class(t_ESSimulator), intent(in) :: self
            !! Instance of the ES simulator
        type(t_Particle), intent(in) :: pcl
            !! Particle to update
        double precision, intent(in) :: dt
            !! Time step for the update
        type(t_CollisionRecord), intent(out) :: record
            !! Collision record
        double precision, intent(out), optional :: unwrapped_position(3)
            !! Collision-clipped endpoint before periodic/reflection mapping.
        type(t_Particle) :: pcl_new
            !! Updated particle

        pcl_new = self%backward(pcl, dt)

        record = self%check_collision(pcl%position, pcl_new%position)

        if (record%is_collided) then
            block
                double precision :: r
                integer :: i

                r = record%t
                pcl_new%q_m = pcl%q_m
                pcl_new%position = pcl%position*(1d0 - r) + pcl_new%position*r
                pcl_new%velocity = pcl%velocity*(1d0 - r) + pcl_new%velocity*r
                pcl_new%t = pcl%t + r*dt

                if (present(unwrapped_position)) unwrapped_position = pcl_new%position

                ! A collision found after a periodic crossing lies outside the
                ! fundamental cell; report it at the same physical point inside.
                do i = 1, 3
                    if (self%boundary_conditions(i) == 0) then
                        pcl_new%position(i) = modulo(pcl_new%position(i), self%domain_length(i))
                    end if
                end do
            end block
            return
        end if

        if (present(unwrapped_position)) unwrapped_position = pcl_new%position
        call self%apply_boundary_condition(pcl_new)
    end function

    function esSimulator_backward(self, particle, dt) result(ret)
        !! Perform an inverse Boris step for dt >= 0, or a forward step for dt < 0.

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

        ! Invert the forward drift before undoing the velocity update.
        if (dt >= 0d0) then
            position_new = particle%position - dt*particle%velocity
        end if

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
            type(t_Particle) :: field_particle

            dt2 = -0.5d0*dt

            ! Copy the base components explicitly. ifx 2023.2 lowers intrinsic
            ! assignment from this polymorphic dummy to a shared temporary
            ! descriptor with -qopenmp. Avoid that concurrent runtime assignment
            ! on the Boris hot path as well as in the probability solver.
            field_particle%t = particle%t
            field_particle%q_m = particle%q_m
            field_particle%position = particle%position
            field_particle%velocity = particle%velocity
            if (dt >= 0d0) field_particle%position = position_new
            ! Sample the same physical point after a periodic boundary crossing.
            ! Keep the unwrapped drift endpoint for collision detection.
            call self%apply_boundary_condition(field_particle)
            eb = self%eb%at(field_particle%position)
            ef(:) = eb(1:3)
            bf(:) = eb(4:6)

            t = bf*particle%q_m*dt2
            s = 2*t/(1d0 + sum(t*t))

            upm = particle%velocity + particle%q_m*ef*dt2

            upa = upm + cross(upm, t)
            upp = upm + cross(upa, s)

            velocity_new = upp + particle%q_m*ef*dt2
        end block

        ! Forward Boris updates velocity before drifting with the new velocity.
        if (dt < 0d0) then
            position_new = particle%position - dt*velocity_new
        end if

        ret%t = particle%t + dt
        ret%q_m = particle%q_m
        ret%position = position_new
        ret%velocity = velocity_new
    end function

    pure function esSimulator_domain_length(self, axis) result(length)
        !! Domain length along `axis` in grid units.
        class(t_ESSimulator), intent(in) :: self
        integer, intent(in) :: axis
        double precision :: length

        select case (axis)
        case (1)
            length = self%nx
        case (2)
            length = self%ny
        case default
            length = self%nz
        end select
    end function

    function esSimulator_check_collision(self, p1, p2) result(record)
        !! First collision of the straight segment p1 -> p2 with the boundaries.
        !!
        !! `p2` is the unwrapped endpoint. Where the segment leaves the
        !! fundamental cell through a periodic face, the remaining part is
        !! physically inside the cell again, so each piece between crossings is
        !! shifted back by whole periods before it is tested. Without this the
        !! part after the crossing was never tested and internal objects there
        !! were passed through (issue #39). `record%t` is the fraction of the
        !! whole segment, and `record%position` lies inside the cell.

        class(t_ESSimulator), intent(in) :: self
            !! Instance of the ES simulator
        double precision, intent(in) :: p1(3)
            !! Segment start (inside the fundamental cell on periodic axes)
        double precision, intent(in) :: p2(3)
            !! Unwrapped segment end
        type(t_CollisionRecord) :: record
            !! First collision along the segment

        double precision :: d(3), length
        double precision, allocatable :: fractions(:)
        integer :: axis, ncross, i, j, k, kfirst, klast, kstep

        d = p2 - p1

        ! Count the periodic planes k*L crossed strictly between p1 and p2.
        ncross = 0
        do axis = 1, 3
            if (self%boundary_conditions(axis) /= 0 .or. d(axis) == 0d0) cycle
            length = self%domain_length(axis)
            call crossing_planes(p1(axis), p2(axis), length, kfirst, klast, kstep)
            ncross = ncross + max(0, (klast - kfirst)/kstep + 1)
        end do

        if (ncross == 0) then
            record = self%boundaries%check_collision(p1, p2)
            return
        end if

        allocate (fractions(0:ncross + 1))
        fractions(0) = 0d0
        fractions(ncross + 1) = 1d0
        i = 0
        do axis = 1, 3
            if (self%boundary_conditions(axis) /= 0 .or. d(axis) == 0d0) cycle
            length = self%domain_length(axis)
            call crossing_planes(p1(axis), p2(axis), length, kfirst, klast, kstep)
            do k = kfirst, klast, kstep
                i = i + 1
                fractions(i) = (dble(k)*length - p1(axis))/d(axis)
            end do
        end do

        ! Insertion sort: crossings are few.
        do i = 2, ncross
            do j = i, 2, -1
                if (fractions(j - 1) <= fractions(j)) exit
                length = fractions(j - 1)
                fractions(j - 1) = fractions(j)
                fractions(j) = length
            end do
        end do

        do i = 0, ncross
            block
                double precision :: a, b, mid(3), shift(3), q1(3), q2(3)

                a = fractions(i)
                b = fractions(i + 1)
                if (b <= a) cycle

                ! Shift this piece back into the fundamental cell.
                mid = p1 + 0.5d0*(a + b)*d
                shift = 0d0
                do axis = 1, 3
                    if (self%boundary_conditions(axis) /= 0) cycle
                    length = self%domain_length(axis)
                    shift(axis) = -floor(mid(axis)/length)*length
                end do
                q1 = p1 + a*d + shift
                q2 = p1 + b*d + shift

                record = self%boundaries%check_collision(q1, q2)
                if (record%is_collided) then
                    record%t = a + record%t*(b - a)
                    return
                end if
            end block
        end do

        record%is_collided = .false.

    contains

        pure subroutine crossing_planes(x1, x2, period, first, last, step)
            !! Integer multiples k of `period` with k*period strictly between x1 and x2,
            !! ordered from x1 toward x2 (empty when first and last straddle step).
            double precision, intent(in) :: x1, x2, period
            integer, intent(out) :: first, last, step

            if (x2 > x1) then
                first = floor(x1/period) + 1
                last = ceiling(x2/period) - 1
                step = 1
            else
                first = ceiling(x1/period) - 1
                last = floor(x2/period) + 1
                step = -1
            end if
        end subroutine
    end function

    function esSimulator_adaptive_dt(self, particle, dt) result(step)
        !! Time step that moves `particle` by about |dt| grid units.
        !!
        !! `dt/|v|` alone diverges at zero velocity and near turning points,
        !! where a finite electric field still moves the particle. Solving
        !! 0.5*a*s**2 + |v|*s = |dt| for s with the local acceleration a keeps
        !! the step finite there and reduces to `dt/|v|` when a is negligible.
        !! When neither velocity nor acceleration can move the particle the
        !! nominal step is returned. The sign of `dt` is preserved.

        class(t_ESSimulator), intent(in) :: self
            !! Instance of the ES simulator
        type(t_Particle), intent(in) :: particle
            !! Particle at the start of the step
        double precision, intent(in) :: dt
            !! Nominal step; |dt| is the target displacement in grid units
        double precision :: step
            !! Signed time step

        double precision :: ds, speed, accel, eb(6)
        type(t_Particle) :: probe

        ds = abs(dt)
        speed = norm2(particle%velocity)

        probe%t = particle%t
        probe%q_m = particle%q_m
        probe%position = particle%position
        probe%velocity = particle%velocity
        call self%apply_boundary_condition(probe)
        eb = self%eb%at(probe%position)
        accel = abs(particle%q_m)*norm2(eb(1:3))

        if (accel*ds > 0d0) then
            ! Positive root of 0.5*a*s**2 + |v|*s - ds = 0, written without
            ! cancellation so that |v| -> 0 gives sqrt(2*ds/a).
            step = 2d0*ds/(speed + sqrt(speed*speed + 2d0*accel*ds))
        else if (speed > 0d0) then
            step = ds/speed
        else
            step = ds
        end if

        step = sign(step, dt)
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
                if (0 <= particle%position(i) .and. particle%position(i) <= nxyz(i)) then
                    cycle
                end if

                particle%position(i) = modulo(2*nxyz(i) - particle%position(i), nxyz(i))
                particle%velocity(i) = -particle%velocity(i)
            end if
        end do
    end subroutine

end module

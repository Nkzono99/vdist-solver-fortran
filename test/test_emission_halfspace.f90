program test_emission_halfspace
    !! Regression for issue #41: an explicit emission patch (nepl > 0) must
    !! emit only toward the side given by the sign of `nemd`. The rear side of
    !! a thin emitting panel sees no source.
    use, intrinsic :: iso_c_binding, only: c_char
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_simulator, only: t_ESSimulator
    use m_solver, only: t_Solver, new_Solver, t_ProbabilityRecord
    use m_particle, only: new_Particle
    use m_probabilities, only: new_HalfSpaceMaxwellianProbability, t_HalfSpaceMaxwellianProbability
    use m_test_helpers, only: assert_close, assert_true

    implicit none

    double precision, parameter :: pi = acos(-1d0)

    call test_half_space_source_is_zero_behind_its_normal()
    call test_patch_emits_only_toward_signed_nemd()

    print *, "test_emission_halfspace: all tests passed."

contains

    subroutine test_half_space_source_is_zero_behind_its_normal()
        type(t_HalfSpaceMaxwellianProbability) :: source
        double precision :: expected

        source = new_HalfSpaceMaxwellianProbability([0d0, 0d0, 0d0], [1d0, 1d0, 1d0], [0d0, -1d0, 0d0], 2d0)
        expected = 2d0*exp(-0.5d0)/(2d0*pi)**1.5d0
        call assert_close("emitted velocity keeps the Maxwellian density", &
                          source%at([0d0, 0d0, 0d0], [0d0, -1d0, 0d0]), expected)
        call assert_close("velocity into the surface has zero density", &
                          source%at([0d0, 0d0, 0d0], [0d0, 1d0, 0d0]), 0d0)
        call assert_close("tangential velocity has zero density", &
                          source%at([0d0, 0d0, 0d0], [1d0, 0d0, 0d0]), 0d0)
    end subroutine

    function build(nemd_value) result(simulator)
        !! Thin absorbing panel at the centre of a [0,4]^3 box with an
        !! emission patch on it whose signed normal is `nemd_value`.
        integer, intent(in) :: nemd_value
        type(t_ESSimulator) :: simulator
        character(len=*), parameter :: filename = 'test-emission-halfspace.inp'
        character(kind=c_char) :: cpath(len(filename))
        double precision :: fields(9, 5, 5, 5), shape(6), box(2, 3)
        integer :: unit, i, axis

        axis = abs(nemd_value)
        box(1, :) = 1d0
        box(2, :) = 3d0
        box(:, axis) = 2d0
        shape = [box(1, 1), box(2, 1), box(1, 2), box(2, 2), box(1, 3), box(2, 3)]

        open (newunit=unit, file=filename, status='replace')
        write (unit, *) '&esorem emflag=0 /'
        write (unit, *) '&plasma wp=10*0, wc=0, phiz=0, phixy=0 /'
        write (unit, *) '&tmgrid nx=4, ny=4, nz=4, dt=0.2 /'
        write (unit, *) '&system nspec=1, npbnd=30*2 /'
        write (unit, *) '&intp qm=10*1, spa=10*0, spe=10*0, speth=10*0,', &
            ' vdri=10*0, vdthz=10*0, vdthxy=10*0, peth=10*1, path=10*1 /'
        write (unit, *) "&ptcond boundary_type='complex', zssurf=-9999,", &
            " boundary_types='rectangle',9*'none', rectangle_shape(:,1)=", shape, ' /'
        write (unit, *) '&emissn nflag_emit=10*1, nepl=1,9*0, nemd=', nemd_value, ',', &
            ' curf=10*1, curfs=15*-1, use_raycast=.false.,', &
            ' xmine=', box(1, 1), ', xmaxe=', box(2, 1), ', ymine=', box(1, 2), ', ymaxe=', box(2, 2), &
            ', zmine=', box(1, 3), ', zmaxe=', box(2, 3), ' /'
        close (unit)
        do i = 1, len(filename)
            cpath(i) = filename(i:i)
        end do
        fields = 0d0
        simulator = create_simulator(cpath, len(filename), 4, 4, 4, fields, 1, 10)
        open (newunit=unit, file=filename, status='old')
        close (unit, status='delete')
    end function

    subroutine test_patch_emits_only_toward_signed_nemd()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_ProbabilityRecord) :: front, rear
        double precision :: expected, position(3), velocity(3), direction
        integer :: nemd_value, axis

        ! curf=1, path=1, zero drift: n = sqrt(2*pi), f = n*exp(-1/2)/(2*pi)**1.5
        expected = exp(-0.5d0)/(2d0*pi)
        do nemd_value = -3, 3
            if (nemd_value == 0) cycle
            axis = abs(nemd_value)
            direction = sign(1d0, dble(nemd_value))
            simulator = build(nemd_value)
            solver = new_Solver(simulator)

            ! Front: observation point on the emitting side, moving away from the panel.
            position = 2d0
            position(axis) = 2d0 + 0.5d0*direction
            velocity = 0d0
            velocity(axis) = direction
            front = solver%calculate_probability(new_Particle(1d0, position, velocity), 0.2d0, 32, .false.)
            call assert_true("front trace reaches the panel", front%is_valid)
            call assert_close("front trace hits the panel at t=0.5", front%t, 0.5d0, 1d-12)
            call assert_close("front side sees the normalised emission source", front%probability, expected, 1d-12)

            ! Rear: mirror image through the panel, moving away from its back face.
            position(axis) = 2d0 - 0.5d0*direction
            velocity(axis) = -direction
            rear = solver%calculate_probability(new_Particle(1d0, position, velocity), 0.2d0, 32, .false.)
            call assert_true("rear trace reaches the panel", rear%is_valid)
            call assert_close("rear trace hits the panel at t=0.5", rear%t, 0.5d0, 1d-12)
            call assert_close("rear side of the patch has no source", rear%probability, 0d0)

            call destroy_simulator(simulator)
        end do
    end subroutine
end program

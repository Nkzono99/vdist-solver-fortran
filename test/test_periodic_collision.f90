program test_periodic_collision
    !! Regression for issue #39: a step that crosses a periodic face must still
    !! collide with internal objects on the part of the segment that lies after
    !! the crossing (physically back inside the fundamental cell).
    use m_particle, only: t_Particle, new_Particle
    use m_field, only: t_VectorFieldGrid, new_VectorFieldGrid
    use m_probabilities, only: tp_Probability, new_ZeroProbability, new_MaxwellianProbability
    use m_simulator, only: t_ESSimulator, new_ESSimulator
    use m_solver, only: t_Solver, new_Solver, t_ProbabilityRecord
    use finbound, only: t_BoundaryList, new_BoundaryList, t_Boundary, &
                        t_PlaneXYZ, new_PlaneZ, t_RectangleXYZ, new_RectangleX, t_CollisionRecord
    use m_test_helpers, only: assert_close, assert_close_vec, assert_true

    implicit none

    integer, parameter :: nx = 8, ny = 8, nz = 4
    double precision, parameter :: pi = acos(-1d0)

    call test_issue_configuration_absorbs_after_positive_crossing()
    call test_negative_crossing_is_absorbed()
    call test_objects_in_the_middle_are_not_hit_by_wrapped_segments()
    call test_collision_after_several_periods()

    print *, "test_periodic_collision: all tests passed."

contains

    function build(sheet_x, sheet_z_height) result(simulator)
        !! x/y periodic, z open with a Maxwellian (sigma=4) source on z=0 and
        !! z=nz, plus an absorbing rectangle x=sheet_x, 2<=y<=6, 0<=z<=height.
        double precision, intent(in) :: sheet_x, sheet_z_height
        type(t_ESSimulator) :: simulator
        type(t_BoundaryList) :: boundaries
        type(tp_Probability), allocatable :: probs(:)
        double precision :: values(6, 0:nx, 0:ny, 0:nz)
        class(t_Boundary), pointer :: pbound
        type(t_PlaneXYZ), pointer :: pplane
        class(t_RectangleXYZ), pointer :: prect

        values = 0d0
        allocate (probs(2))
        allocate (probs(1)%ref, source=new_ZeroProbability())
        allocate (probs(2)%ref, source=new_MaxwellianProbability([0d0, 0d0, 0d0], [4d0, 4d0, 4d0]))

        boundaries = new_BoundaryList()
        allocate (prect, source=new_RectangleX([sheet_x, 2d0, 0d0], 4d0, sheet_z_height))
        pbound => prect
        pbound%material%tag = 1
        call boundaries%add_boundary(pbound)
        allocate (pplane)
        pplane = new_PlaneZ(0d0)
        pbound => pplane
        pbound%material%tag = 2
        call boundaries%add_boundary(pbound)
        allocate (pplane)
        pplane = new_PlaneZ(dble(nz))
        pbound => pplane
        pbound%material%tag = 2
        call boundaries%add_boundary(pbound)

        simulator = new_ESSimulator(nx, ny, nz, [0, 0, 2], &
                                    new_VectorFieldGrid(6, nx, ny, nz, values), boundaries, probs)
    end function

    subroutine test_issue_configuration_absorbs_after_positive_crossing()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl, hit
        type(t_CollisionRecord) :: collision
        type(t_ProbabilityRecord) :: coarse, fine
        double precision :: endpoint(3)

        simulator = build(0.5d0, 4d0)
        solver = new_Solver(simulator)
        pcl = new_Particle(1d0, [7.9d0, 4d0, 2.5d0], [-4d0, 0d0, -1d0])

        ! One dt=0.2 step runs 7.9 -> 8.7 unwrapped, i.e. 7.9 -> 8 then 0 -> 0.7
        ! inside the cell, and meets the sheet at x=0.5 after 0.15.
        hit = simulator%update(pcl, 0.2d0, collision, endpoint)
        call assert_true("seam-crossing step collides with the sheet", collision%is_collided)
        call assert_close("collision fraction covers the whole step", collision%t, 0.75d0, 1d-12)
        call assert_close_vec("unwrapped endpoint is clipped at the collision", endpoint, [8.5d0, 4d0, 2.65d0], 1d-12)
        call assert_close_vec("returned position is wrapped into the cell", hit%position, [0.5d0, 4d0, 2.65d0], 1d-12)
        call assert_close("collision time", hit%t, 0.15d0, 1d-12)
        call assert_close_vec("collision record position lies inside the cell", &
                              collision%position, [0.5d0, 4d0, 2.65d0], 1d-9)

        coarse = solver%calculate_probability(pcl, 0.2d0, 1000, .false.)
        fine = solver%calculate_probability(pcl, 0.013d0, 1000, .false.)
        call assert_true("coarse trace completes", coarse%is_valid)
        call assert_true("fine trace completes", fine%is_valid)
        call assert_close("fine step absorbs on the sheet", fine%probability, 0d0)
        call assert_close("coarse step absorbs on the sheet too", coarse%probability, 0d0)
        call assert_close("coarse and fine collision times agree", coarse%t, fine%t, 1d-12)
        call assert_close_vec("coarse collision position", coarse%particle%position, [0.5d0, 4d0, 2.65d0], 1d-12)
    end subroutine

    subroutine test_negative_crossing_is_absorbed()
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_ProbabilityRecord) :: record

        simulator = build(7.5d0, 4d0)
        solver = new_Solver(simulator)
        record = solver%calculate_probability(new_Particle(1d0, [0.1d0, 4d0, 2.5d0], [4d0, 0d0, -1d0]), &
                                              0.2d0, 1000, .false.)
        call assert_true("negative crossing trace completes", record%is_valid)
        call assert_close("negative crossing absorbs on the sheet", record%probability, 0d0)
        call assert_close("negative crossing collision time", record%t, 0.15d0, 1d-12)
        call assert_close_vec("negative crossing collision position", record%particle%position, &
                              [7.5d0, 4d0, 2.65d0], 1d-12)
    end subroutine

    subroutine test_objects_in_the_middle_are_not_hit_by_wrapped_segments()
        !! Testing the wrapped endpoints as one straight line would cross the
        !! whole cell; the split must only test the pieces actually travelled.
        type(t_ESSimulator) :: simulator
        type(t_Particle) :: pcl, moved
        type(t_CollisionRecord) :: collision
        type(t_ProbabilityRecord) :: record
        type(t_Solver) :: solver
        double precision :: expected

        simulator = build(4.5d0, 4d0)
        solver = new_Solver(simulator)
        pcl = new_Particle(1d0, [7.9d0, 4d0, 2.5d0], [-4d0, 0d0, -1d0])
        moved = simulator%update(pcl, 0.2d0, collision)
        call assert_true("no collision with a sheet in the middle of the cell", .not. collision%is_collided)
        call assert_close_vec("particle wraps normally", moved%position, [0.7d0, 4d0, 2.7d0], 1d-12)

        ! The trace keeps moving toward +x inside the cell and reaches the sheet
        ! at x=4.5 (unwrapped 12.5) after 1.15, with z=3.65 still inside the box,
        ! before the z=4 source at 1.5.
        record = solver%calculate_probability(pcl, 0.2d0, 1000, .false.)
        call assert_true("trace completes", record%is_valid)
        call assert_close("sheet in the middle absorbs the later pass", record%probability, 0d0)
        call assert_close("later-pass collision time", record%t, 1.15d0, 1d-12)
        call assert_close_vec("later-pass collision position", record%particle%position, [4.5d0, 4d0, 3.65d0], 1d-12)

        ! Without the sheet the same trace reaches the z=4 source after 1.5.
        simulator = build(4.5d0, 0d0)
        solver = new_Solver(simulator)
        record = solver%calculate_probability(pcl, 0.2d0, 1000, .false.)
        expected = exp(-17d0/32d0)/((2d0*pi)**1.5d0*64d0)
        call assert_true("trace reaches the open boundary", record%is_valid)
        call assert_close("open-boundary Maxwellian is returned", record%probability, expected, 1d-15)
        call assert_close("open-boundary arrival time", record%t, 1.5d0, 1d-12)
    end subroutine

    subroutine test_collision_after_several_periods()
        !! A single large step crossing x=8 and x=16 (unwrapped) must find the
        !! sheet image that is actually hit (the second one, where z is low enough).
        type(t_ESSimulator) :: simulator
        type(t_Solver) :: solver
        type(t_Particle) :: pcl, hit
        type(t_CollisionRecord) :: collision
        type(t_ProbabilityRecord) :: record
        double precision :: endpoint(3)

        simulator = build(0.5d0, 1d0)
        solver = new_Solver(simulator)
        pcl = new_Particle(1d0, [7.9d0, 4d0, 2.5d0], [-4d0, 0d0, 1d0])

        ! dt=4 (backward: position - dt*v): x 7.9 -> 23.9 crosses 8 and 16, while
        ! z 2.5 -> -1.5 leaves through z=0 at s=0.625. The sheet images at x=8.5
        ! (s=0.0375, z=2.35: above the 1-high sheet) and x=16.5 (s=0.5375, z=0.35: hit)
        ! make the second crossing the one that matters.
        hit = simulator%update(pcl, 4d0, collision, endpoint)
        call assert_true("second periodic image is hit", collision%is_collided)
        call assert_close("fraction of the whole step", collision%t, 0.5375d0, 1d-12)
        call assert_close_vec("unwrapped clipped endpoint", endpoint, [16.5d0, 4d0, 0.35d0], 1d-12)
        call assert_close_vec("wrapped collision position", hit%position, [0.5d0, 4d0, 0.35d0], 1d-12)

        record = solver%calculate_probability(pcl, 4d0, 10, .false.)
        call assert_true("trace completes", record%is_valid)
        call assert_close("absorbed on the second image", record%probability, 0d0)
        call assert_close("collision time after two crossings", record%t, 2.15d0, 1d-12)
    end subroutine
end program

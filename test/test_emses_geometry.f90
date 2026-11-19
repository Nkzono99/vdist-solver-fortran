program test_emses_geometry
    use m_allcom, only: reset_input_parameters, nx, ny, nz, boundary_type, boundary_types, &
                       disk_origin, disk_radius, disk_inner_radius, disk_height, &
                       rectangle_corners, zssurf
    use m_collision_boundary_input_config, only: load_collision_boundary_config, t_collision_boundary_config
    use m_emses_boundaries, only: create_simple_collision_boundaries
    use m_namelist, only: read_namelist
    use finbound, only: t_BoundaryList, t_HitRecord, new_Ray
    use m_test_helpers, only: assert_true, assert_close, assert_equal_int
    implicit none
    type(t_BoundaryList) :: boundaries
    type(t_HitRecord) :: hit
    type(t_collision_boundary_config) :: cfg
    integer :: domain(2, 3), unit

    call reset_input_parameters()
    nx = 20; ny = 20; nz = 20
    domain(1, :) = 0
    domain(2, :) = [nx, ny, nz]
    boundary_type = 'complex'
    boundary_types(1) = 'diskz'
    disk_origin(:, 1) = [5d0, 5d0, 2d0]
    disk_radius(1) = 3d0
    disk_inner_radius(1) = 1d0
    disk_height(1) = 0.5d0
    boundaries = create_simple_collision_boundaries(domain, tag=7)
    call assert_equal_int('disk consists of four boundaries', boundaries%nboundaries, 4)
    hit = boundaries%hit(new_Ray([7d0, 5d0, 10d0], [0d0, 0d0, -1d0]))
    call assert_true('EMSES annular disk top cap is hit', hit%is_hit)
    call assert_close('EMSES cap height', hit%position(3), 2.5d0)
    call assert_equal_int('probability tag overrides geometry tag', hit%material%tag, 7)
    hit = boundaries%hit(new_Ray([5d0, 5d0, 10d0], [0d0, 0d0, -1d0]))
    call assert_true('disk aperture remains open', .not. hit%is_hit)
    call boundaries%destroy()

    ! New EMSES shape is passed through the adapter, including its corner array.
    boundary_types(1) = 'rectangle-corners'
    rectangle_corners(:, 1, 1) = [1d0, 1d0, 3d0]
    rectangle_corners(:, 2, 1) = [4d0, 1d0, 3d0]
    rectangle_corners(:, 3, 1) = [4d0, 4d0, 3d0]
    rectangle_corners(:, 4, 1) = [1d0, 4d0, 3d0]
    call load_collision_boundary_config(cfg)
    boundaries = create_simple_collision_boundaries(domain)
    hit = boundaries%hit(new_Ray([2d0, 2d0, 10d0], [0d0, 0d0, -1d0]))
    call assert_true('rectangle-corners reaches the EMSES builder', hit%is_hit)
    call assert_close('quad height', hit%position(3), 3d0)
    call boundaries%destroy()
    call reset_input_parameters()
    call assert_true('resolved config owns its arrays', cfg%boundary_types(1) == 'rectangle-corners')
    call assert_close('new input array resets', sum(abs(rectangle_corners)), 0d0)

    open(newunit=unit, file='test_emses_geometry.inp', status='replace')
    write(unit, '(a)') '&esorem /', '&plasma /', '&tmgrid nx=20, ny=20, nz=20 /', &
        '&system /', '&intp /', '&emissn /', &
        "&ptcond boundary_type='complex', boundary_types(1)='rectangle-corners',", &
        'rectangle_corners(:,1,1)=1,1,3, rectangle_corners(:,2,1)=4,1,3,', &
        'rectangle_corners(:,3,1)=4,4,3, rectangle_corners(:,4,1)=1,4,3 /'
    close(unit)
    call read_namelist('test_emses_geometry.inp')
    open(newunit=unit, file='test_emses_geometry.inp', status='old')
    close(unit, status='delete')
    boundaries = create_simple_collision_boundaries(domain)
    hit = boundaries%hit(new_Ray([2d0, 2d0, 10d0], [0d0, 0d0, -1d0]))
    call assert_true('rectangle-corners namelist is supported', hit%is_hit)
    call assert_close('quad loaded from namelist', hit%position(3), 3d0)
    call boundaries%destroy()

    ! A subsequent construction cannot reuse the previous object's geometry.
    nx = 20; ny = 20; nz = 20
    boundary_type = 'flat-surface'
    zssurf = 4d0
    boundaries = create_simple_collision_boundaries(domain)
    hit = boundaries%hit(new_Ray([5d0, 5d0, 10d0], [0d0, 0d0, -1d0]))
    call assert_true('surface after reset', hit%is_hit)
    call assert_close('surface height after reset', hit%position(3), 4d0)
    call boundaries%destroy()
    print *, 'test_emses_geometry: all tests passed'
end program test_emses_geometry

!> Object-shaped internal boundaries for the complex-geometry mode.
!>
!> When `boundary_type == "complex"`, each element of
!> `boundary_types(1..nboundary_types)` selects one of the `add_*`
!> contained routines below to append boundary pieces to the shared
!> list. See the long namelist comment in particle_collision.F90 for
!> the supported shape names and their parameter arrays.
module m_objects
  use iso_fortran_env, only: error_unit
  use finbound, only: new_CircleXYZ, new_DonutXYZ, new_Sphere, new_cylinderXYZ, new_rectangleX, new_rectangleY, &
                      new_planeXYZWithCircleHole, &
                      new_Rectangle, t_Rectangle, new_rectangleZ, t_Boundary, t_BoundaryList, t_CircleXYZ, &
                      t_CylinderXYZ, t_DonutXYZ, &
                      t_PlaneXYZWithCircleHole, t_RectangleXYZ, t_Sphere
  use m_collision_boundary_input_config, only: load_collision_boundary_config, t_collision_boundary_config
  use m_collision_boundary_types, only: is_collision_boundary_entry_empty, is_collision_object_boundary_type, &
                                        require_supported_collision_boundary_entry
  use m_collision_accelerator, only: aabb_from_axis_plane, aabb_from_circle_xyz, aabb_from_cylinder_xyz, &
                                     aabb_from_quad, aabb_from_rectangle_xyz, aabb_from_sphere, &
                                     t_collision_accelerator
  use m_kinds, only: dp, ip
  use m_boundary_helpers, only: try_add_boundary
  implicit none
  private
  public :: add_objects

contains

  subroutine add_objects(boundaries, isdoms, offsets, accelerator)
    type(t_BoundaryList), intent(inout) :: boundaries
    integer(kind=ip), intent(in) :: isdoms(2, 3)
    real(kind=dp), intent(in) :: offsets(3)
    type(t_collision_accelerator), intent(inout), optional :: accelerator

    type(t_collision_boundary_config) :: cfg
    real(kind=dp) :: sdoms(2, 3)
    character(len=30) :: boundary_name
    integer(kind=ip) :: itype

    call load_collision_boundary_config(cfg)

    if (cfg%boundary_type /= "complex") return
    sdoms(1:2, 1:3) = dble(isdoms(1:2, 1:3))

    do itype = 1, cfg%nboundary_types
      boundary_name = trim(cfg%boundary_types(itype))
      if (is_collision_boundary_entry_empty(boundary_name)) cycle
      call require_supported_collision_boundary_entry(boundary_name, "adding collision objects", itype)
      if (.not. is_collision_object_boundary_type(boundary_name)) cycle

      select case (boundary_name)
      case ('rectangle'); call add_rectangle
      case ('rectangle-corners'); call add_rectangle_corners
      case ('sphere'); call add_sphere
      case ('circlex'); call add_circleXYZ(1)
      case ('circley'); call add_circleXYZ(2)
      case ('circlez'); call add_circleXYZ(3)
      case ('cuboid'); call add_cuboid
      case ('cylinderx'); call add_cylinderXYZ(1)
      case ('cylindery'); call add_cylinderXYZ(2)
      case ('cylinderz'); call add_cylinderXYZ(3)
      case ('open-cylinderx'); call add_open_cylinderXYZ(1)
      case ('open-cylindery'); call add_open_cylinderXYZ(2)
      case ('open-cylinderz'); call add_open_cylinderXYZ(3)
      case ('diskx'); call add_disk(1)
      case ('disky'); call add_disk(2)
      case ('diskz'); call add_disk(3)
      case ('plane-with-circlex'); call add_plane_with_circleXYZ(1)
      case ('plane-with-circley'); call add_plane_with_circleXYZ(2)
      case ('plane-with-circlez'); call add_plane_with_circleXYZ(3)
      case default; call stop_unimplemented_boundary_type(boundary_name)
      end select
    end do

  contains

    !> Allocate a rectangle boundary and hand it to try_add_boundary
    !> (lifetime afterwards is in the boundary list or freed on overlap miss).
    subroutine add_rectangle
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect_xyz
      real(kind=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
      logical :: flat_x, flat_y, flat_z
      integer(kind=ip) :: nflat

      xmin = cfg%rectangle_shape(1, itype) + offsets(1)
      xmax = cfg%rectangle_shape(2, itype) + offsets(1)
      ymin = cfg%rectangle_shape(3, itype) + offsets(2)
      ymax = cfg%rectangle_shape(4, itype) + offsets(2)
      zmin = cfg%rectangle_shape(5, itype) + offsets(3)
      zmax = cfg%rectangle_shape(6, itype) + offsets(3)

      if (xmax < xmin .or. ymax < ymin .or. zmax < zmin) then
        call stop_invalid_boundary_parameter("rectangle_shape upper bounds must be greater than or equal to lower bounds")
      end if

      flat_x = xmin == xmax
      flat_y = ymin == ymax
      flat_z = zmin == zmax
      nflat = 0
      if (flat_x) nflat = nflat + 1
      if (flat_y) nflat = nflat + 1
      if (flat_z) nflat = nflat + 1
      if (nflat /= 1) then
        call stop_invalid_boundary_parameter("rectangle_shape must have exactly one zero-thickness axis")
      end if

      if (flat_x) then
        call require_positive_width(ymax - ymin, "rectangle y width")
        call require_positive_width(zmax - zmin, "rectangle z width")
        allocate (prect_xyz)
        prect_xyz = new_rectangleX([xmin, ymin, zmin], ymax - ymin, zmax - zmin)
        pbound => prect_xyz
      else if (flat_y) then
        call require_positive_width(zmax - zmin, "rectangle z width")
        call require_positive_width(xmax - xmin, "rectangle x width")
        allocate (prect_xyz)
        prect_xyz = new_rectangleY([xmin, ymin, zmin], zmax - zmin, xmax - xmin)
        pbound => prect_xyz
      else if (flat_z) then
        call require_positive_width(xmax - xmin, "rectangle x width")
        call require_positive_width(ymax - ymin, "rectangle y width")
        allocate (prect_xyz)
        prect_xyz = new_rectangleZ([xmin, ymin, zmin], xmax - xmin, ymax - ymin)
        pbound => prect_xyz
      end if
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_rectangle_xyz(prect_xyz%axis, prect_xyz%origin, prect_xyz%w1, prect_xyz%w2))
    end subroutine

    !> Allocate a general rectangle defined by four corner points.
    subroutine add_rectangle_corners
      class(t_Boundary), pointer :: pbound
      type(t_Rectangle), pointer :: prect
      real(kind=dp) :: corners(3, 4)
      integer(kind=ip) :: icorner

      corners(:, :) = cfg%rectangle_corners(:, :, itype)
      do icorner = 1, 4
        corners(:, icorner) = corners(:, icorner) + offsets(:)
      end do

      allocate (prect)
      prect = new_Rectangle(corners)
      pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_quad(corners))
    end subroutine

    subroutine add_circleXYZ(axis)
      integer(kind=ip), intent(in) :: axis
      class(t_Boundary), pointer :: pbound
      type(t_CircleXYZ), pointer :: pcircle

      call require_positive_width(cfg%circle_radius(itype), "circle_radius")
      allocate (pcircle)
      pcircle = new_CircleXYZ(axis, cfg%circle_origin(:, itype) + offsets(:), cfg%circle_radius(itype))
      pbound => pcircle
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_circle_xyz(axis, pcircle%origin, pcircle%radius))
    end subroutine

    !> Six rectangular faces of a cuboid. The coordinate layout is:
    !>
    !>             ------------- (max)
    !>           / |           /
    !>          /  |    6     / |
    !>         /   |      5  /  |
    !>    z ^  -------------  4 |
    !>      | | 1   --------|---/  ^ y
    !>        |   /  2      |  /  /
    !>        |  /      3   | /
    !>        | /           |/
    !>  (min)  -------------/
    !>                      -> x
    subroutine add_cuboid
      real(kind=dp) :: xmin, xmax, ymin, ymax, zmin, zmax, wx, wy, wz

      xmin = cfg%cuboid_shape(1, itype) + offsets(1)
      xmax = cfg%cuboid_shape(2, itype) + offsets(1)
      ymin = cfg%cuboid_shape(3, itype) + offsets(2)
      ymax = cfg%cuboid_shape(4, itype) + offsets(2)
      zmin = cfg%cuboid_shape(5, itype) + offsets(3)
      zmax = cfg%cuboid_shape(6, itype) + offsets(3)
      wx = xmax - xmin; wy = ymax - ymin; wz = zmax - zmin
      call require_positive_width(wx, "cuboid x width")
      call require_positive_width(wy, "cuboid y width")
      call require_positive_width(wz, "cuboid z width")

      call add_rect_X([xmin, ymin, zmin], wy, wz)   ! face 1 (x = xmin)
      call add_rect_Y([xmin, ymin, zmin], wz, wx)   ! face 2 (y = ymin)
      call add_rect_Z([xmin, ymin, zmin], wx, wy)   ! face 3 (z = zmin)
      call add_rect_X([xmax, ymin, zmin], wy, wz)   ! face 4 (x = xmax)
      call add_rect_Y([xmin, ymax, zmin], wz, wx)   ! face 5 (y = ymax)
      call add_rect_Z([xmin, ymin, zmax], wx, wy)   ! face 6 (z = zmax)
    end subroutine

    subroutine add_rect_X(origin, w1, w2)
      real(kind=dp), intent(in) :: origin(3), w1, w2
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect
      allocate (prect); prect = new_rectangleX(origin, w1, w2); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_rectangle_xyz(1_ip, origin, w1, w2))
    end subroutine

    subroutine add_rect_Y(origin, w1, w2)
      real(kind=dp), intent(in) :: origin(3), w1, w2
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect
      allocate (prect); prect = new_rectangleY(origin, w1, w2); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_rectangle_xyz(2_ip, origin, w1, w2))
    end subroutine

    subroutine add_rect_Z(origin, w1, w2)
      real(kind=dp), intent(in) :: origin(3), w1, w2
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect
      allocate (prect); prect = new_rectangleZ(origin, w1, w2); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_rectangle_xyz(3_ip, origin, w1, w2))
    end subroutine

    subroutine add_sphere
      class(t_Boundary), pointer :: pbound
      type(t_Sphere), pointer :: psphere
      call require_positive_width(cfg%sphere_radius(itype), "sphere_radius")
      allocate (psphere)
      psphere = new_Sphere(cfg%sphere_origin(:, itype) + offsets(:), cfg%sphere_radius(itype))
      pbound => psphere
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_sphere(psphere%origin, psphere%radius))
    end subroutine

    subroutine add_open_cylinderXYZ(axis)
      integer(kind=ip), intent(in) :: axis
      class(t_Boundary), pointer :: pbound
      type(t_CylinderXYZ), pointer :: pcylinder

      call require_positive_width(cfg%cylinder_radius(itype), "cylinder_radius")
      call require_positive_width(cfg%cylinder_height(itype), "cylinder_height")
      allocate (pcylinder)
      pcylinder = new_cylinderXYZ(axis, cfg%cylinder_origin(:, itype) + offsets(:), &
                                  cfg%cylinder_radius(itype), cfg%cylinder_height(itype))
      pbound => pcylinder
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_cylinder_xyz(axis, pcylinder%origin, pcylinder%radius, pcylinder%height))
    end subroutine

    !> Closed cylinder: the tube surface + two end-caps.
    subroutine add_cylinderXYZ(axis)
      integer(kind=ip), intent(in) :: axis
      class(t_Boundary), pointer :: pbound
      type(t_CylinderXYZ), pointer :: pcylinder
      type(t_CircleXYZ), pointer :: pcircle
      real(kind=dp) :: lower_origin(3), upper_origin(3), radius, height

      height = cfg%cylinder_height(itype)
      lower_origin(1:3) = cfg%cylinder_origin(1:3, itype) + offsets(:)
      upper_origin(1:3) = lower_origin(1:3)
      upper_origin(axis) = lower_origin(axis) + height
      radius = cfg%cylinder_radius(itype)
      call require_positive_width(radius, "cylinder_radius")
      call require_positive_width(height, "cylinder_height")

      allocate (pcylinder)
      pcylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
      pbound => pcylinder
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_cylinder_xyz(axis, lower_origin, radius, height))

      allocate (pcircle)
      pcircle = new_CircleXYZ(axis, lower_origin, radius)
      pbound => pcircle
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_circle_xyz(axis, lower_origin, radius))

      allocate (pcircle)
      pcircle = new_CircleXYZ(axis, upper_origin, radius)
      pbound => pcircle
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_circle_xyz(axis, upper_origin, radius))
    end subroutine

    subroutine add_plane_with_circleXYZ(axis)
      integer(kind=ip), intent(in) :: axis
      class(t_Boundary), pointer :: pbound
      type(t_PlaneXYZWithCircleHole), pointer :: pplane
      real(kind=dp) :: origin(3)

      ! A plane-with-circle is an infinite plane. Periodic copies are not
      ! supported because a shifted infinite plane would cover another tile.
      if (any(offsets /= 0.0_dp)) return

      allocate (pplane)
      origin(:) = cfg%plane_with_circle_origin(:, itype)
      call require_positive_width(cfg%plane_with_circle_radius(itype), "plane_with_circle_radius")
      pplane = new_planeXYZWithCircleHole(axis, origin, cfg%plane_with_circle_radius(itype))
      pbound => pplane
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_axis_plane(axis, origin(axis)))
    end subroutine

    !> Disk = outer cylinder + inner cylinder + two donut end-caps.
    subroutine add_disk(axis)
      integer(kind=ip), intent(in) :: axis
      class(t_Boundary), pointer :: pbound
      type(t_CylinderXYZ), pointer :: pouter_cylinder, pinner_cylinder
      type(t_DonutXYZ), pointer :: plower_donut, pupper_donut
      real(kind=dp) :: lower_origin(3), upper_origin(3), radius, inner_radius, height

      height = cfg%disk_height(itype)
      lower_origin(1:3) = cfg%disk_origin(1:3, itype) + offsets(:)
      upper_origin(1:3) = lower_origin(1:3)
      upper_origin(axis) = upper_origin(axis) + height
      radius = cfg%disk_radius(itype)
      inner_radius = cfg%disk_inner_radius(itype)
      call require_positive_width(radius, "disk_radius")
      call require_positive_width(height, "disk_height")
      call require_nonnegative_width(inner_radius, "disk_inner_radius")
      if (inner_radius >= radius) then
        call stop_invalid_boundary_parameter("disk_inner_radius must be smaller than disk_radius")
      end if

      allocate (pouter_cylinder)
      pouter_cylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
      pbound => pouter_cylinder
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_cylinder_xyz(axis, lower_origin, radius, height))

      allocate (pinner_cylinder)
      pinner_cylinder = new_cylinderXYZ(axis, lower_origin, inner_radius, height)
      pbound => pinner_cylinder
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_cylinder_xyz(axis, lower_origin, inner_radius, height))

      allocate (plower_donut)
      plower_donut = new_DonutXYZ(axis, lower_origin, radius, inner_radius)
      pbound => plower_donut
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_circle_xyz(axis, lower_origin, radius))

      allocate (pupper_donut)
      pupper_donut = new_DonutXYZ(axis, upper_origin, radius, inner_radius)
      pbound => pupper_donut
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_circle_xyz(axis, upper_origin, radius))
    end subroutine

    subroutine require_positive_width(value, name)
      real(kind=dp), intent(in) :: value
      character(len=*), intent(in) :: name

      if (value <= 0.0d0) call stop_invalid_boundary_parameter(trim(name)//" must be positive")
    end subroutine require_positive_width

    subroutine require_nonnegative_width(value, name)
      real(kind=dp), intent(in) :: value
      character(len=*), intent(in) :: name

      if (value < 0.0d0) call stop_invalid_boundary_parameter(trim(name)//" must be non-negative")
    end subroutine require_nonnegative_width

    subroutine stop_invalid_boundary_parameter(message)
      character(len=*), intent(in) :: message

      write (error_unit, '(a,a,a,i0,a)') "Error.Input: ", trim(message), " at boundary_types(", itype, ")"
      error stop 1
    end subroutine stop_invalid_boundary_parameter

    subroutine stop_unimplemented_boundary_type(name)
      character(len=*), intent(in) :: name

      write (error_unit, '(2a)') "Error.Internal: collision object boundary type is registered but not implemented: ", &
        trim(name)
      error stop 1
    end subroutine stop_unimplemented_boundary_type

  end subroutine add_objects

end module m_objects

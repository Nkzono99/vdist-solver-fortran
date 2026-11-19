!> Flat / dented / holed surface boundaries. Delegated to by
!> `m_particle_collision::create_simple_collision_boundaries` when
!> `boundary_type` (or any of `boundary_types(:)` in complex mode)
!> names one of the surface shapes handled below.
module m_surfaces
  use iso_fortran_env, only: error_unit
  use finbound, only: new_CircleZ, new_EllipsoidZ, new_HyperboloidZ, new_cylinderZ, &
                      new_planeXYZWithCircleHoleZ, new_planeZ, new_rectangleX, new_rectangleY, &
                      new_Rectangle, t_Rectangle, new_rectangleZ, t_Boundary, t_BoundaryList, t_CircleXYZ, t_CylinderXYZ, &
                      t_CollisionRecord, t_HitRecord, t_Ray, &
                      t_EllipsoidXYZ, t_HyperboloidXYZ, t_PlaneXYZ, t_PlaneXYZWithCircleHole, &
                      t_RectangleXYZ
  use m_collision_boundary_input_config, only: load_collision_boundary_config, t_collision_boundary_config
  use m_collision_boundary_types, only: is_collision_boundary_entry_empty, is_collision_surface_boundary_type, &
                                        require_supported_collision_boundary_entry, &
                                        require_supported_simple_collision_boundary_type
  use m_collision_accelerator, only: aabb_from_axis_plane, aabb_from_bounded_axis_plane, aabb_from_centered_profile, &
                                     aabb_from_circle_xyz, aabb_from_cylinder_xyz, aabb_from_quad, &
                                     aabb_from_rectangle_xyz, t_collision_accelerator
  use m_kinds, only: dp, ip
  use m_boundary_helpers, only: try_add_boundary
  implicit none
  private
  public :: add_surfaces

  !> Extend the covered planes well beyond the simulation box so that
  !> particles leaving the domain still collide before exiting. Kept
  !> at 99999 grid cells for bug-compat with the original code.
  real(kind=dp), parameter :: PLANE_PADDING = 99999

  !> Finite z-plane with a circular aperture. Used for periodic raycast
  !> copies so a shifted tile's infinite plane cannot cover another tile's hole.
  type, extends(t_Boundary) :: t_BoundedPlaneZWithCircleHole
    real(kind=dp) :: origin(3)
    real(kind=dp) :: radius
    real(kind=dp) :: xmin, xmax, ymin, ymax
  contains
    procedure :: check_collision => bounded_plane_z_circle_hole_check_collision
    procedure :: hit => bounded_plane_z_circle_hole_hit
    procedure :: is_overlap => bounded_plane_z_circle_hole_is_overlap
    procedure :: pnormal => bounded_plane_z_circle_hole_pnormal
  end type t_BoundedPlaneZWithCircleHole

contains

  pure function new_bounded_plane_z_with_circle_hole(origin, radius, xmin, xmax, ymin, ymax) result(obj)
    real(kind=dp), intent(in) :: origin(3)
    real(kind=dp), intent(in) :: radius
    real(kind=dp), intent(in) :: xmin, xmax, ymin, ymax
    type(t_BoundedPlaneZWithCircleHole) :: obj

    obj%origin(:) = origin(:)
    obj%radius = radius
    obj%xmin = xmin; obj%xmax = xmax
    obj%ymin = ymin; obj%ymax = ymax
  end function new_bounded_plane_z_with_circle_hole

  pure function bounded_plane_z_circle_hole_check_collision(self, p1, p2) result(record)
    class(t_BoundedPlaneZWithCircleHole), intent(in) :: self
    real(kind=dp), intent(in) :: p1(3), p2(3)
    type(t_CollisionRecord) :: record

    real(kind=dp) :: d1, d2, r
    real(kind=dp) :: pos_collided(3)

    d1 = p1(3) - self%origin(3)
    d2 = p2(3) - self%origin(3)
    if (d1*d2 >= 0.0d0) then
      record%is_collided = .false.
      return
    end if

    r = abs(d1)/(abs(d1) + abs(d2))
    pos_collided(:) = (p2(:) - p1(:))*r + p1(:)
    if (.not. is_solid_bounded_circle_hole_point(self, pos_collided)) then
      record%is_collided = .false.
      return
    end if

    record%is_collided = .true.
    record%t = r
    record%position(:) = pos_collided(:)
    record%material = self%material
  end function bounded_plane_z_circle_hole_check_collision

  pure function bounded_plane_z_circle_hole_hit(self, ray) result(hit_record)
    class(t_BoundedPlaneZWithCircleHole), intent(in) :: self
    type(t_Ray), intent(in) :: ray
    type(t_HitRecord) :: hit_record

    real(kind=dp) :: dist, dir, t
    real(kind=dp) :: pos_hit(3)

    dist = self%origin(3) - ray%origin(3)
    dir = ray%direction(3)
    if (dist*dir <= 0.0d0) then
      hit_record%is_hit = .false.
      return
    end if

    t = abs(dist/dir)
    pos_hit(:) = ray%origin(:) + ray%direction(:)*t
    if (.not. is_solid_bounded_circle_hole_point(self, pos_hit)) then
      hit_record%is_hit = .false.
      return
    end if

    hit_record%is_hit = .true.
    hit_record%t = t
    hit_record%position(:) = pos_hit(:)
    hit_record%n(:) = self%normal(pos_hit(:), ray%origin(:))
    hit_record%material = self%material
  end function bounded_plane_z_circle_hole_hit

  pure function bounded_plane_z_circle_hole_is_overlap(self, sdoms, extent) result(is_overlap)
    class(t_BoundedPlaneZWithCircleHole), intent(in) :: self
    real(kind=dp), intent(in) :: sdoms(2, 3)
    real(kind=dp), intent(in), optional :: extent(2, 3)
    logical :: is_overlap

    real(kind=dp) :: extent_(2, 3)
    real(kind=dp) :: sdoms_(2, 3)

    extent_(:, :) = 0.0d0
    if (present(extent)) extent_(:, :) = extent(:, :)
    sdoms_(1, :) = sdoms(1, :) - extent_(1, :)
    sdoms_(2, :) = sdoms(2, :) + extent_(2, :)

    is_overlap = sdoms_(1, 3) <= self%origin(3) .and. self%origin(3) <= sdoms_(2, 3) &
                 .and. self%xmin <= sdoms_(2, 1) .and. sdoms_(1, 1) <= self%xmax &
                 .and. self%ymin <= sdoms_(2, 2) .and. sdoms_(1, 2) <= self%ymax
  end function bounded_plane_z_circle_hole_is_overlap

  pure function bounded_plane_z_circle_hole_pnormal(self, position) result(pnormal)
    class(t_BoundedPlaneZWithCircleHole), intent(in) :: self
    real(kind=dp), intent(in) :: position(3)
    real(kind=dp) :: pnormal(3)

    pnormal(:) = 0.0d0
    pnormal(3) = 1.0d0
  end function bounded_plane_z_circle_hole_pnormal

  pure function is_solid_bounded_circle_hole_point(self, position) result(is_solid)
    class(t_BoundedPlaneZWithCircleHole), intent(in) :: self
    real(kind=dp), intent(in) :: position(3)
    logical :: is_solid

    real(kind=dp) :: dx, dy

    if (position(1) < self%xmin .or. self%xmax < position(1) &
        .or. position(2) < self%ymin .or. self%ymax < position(2)) then
      is_solid = .false.
      return
    end if

    dx = position(1) - self%origin(1)
    dy = position(2) - self%origin(2)
    is_solid = dx*dx + dy*dy >= self%radius*self%radius
  end function is_solid_bounded_circle_hole_point

  subroutine add_surfaces(boundaries, isdoms, cover_all, offsets, accelerator)
    type(t_BoundaryList), intent(inout) :: boundaries
    integer(kind=ip), intent(in) :: isdoms(2, 3)
    logical, intent(in), optional :: cover_all
    real(kind=dp), intent(in), optional :: offsets(3)
    type(t_collision_accelerator), intent(inout), optional :: accelerator

    type(t_collision_boundary_config) :: cfg
    real(kind=dp) :: boundary_offsets(3)
    real(kind=dp) :: xl, yl, zl, xu, yu, zu
    real(kind=dp) :: sdoms(2, 3)
    integer(kind=ip) :: itype
    logical :: is_possible_to_be_covered
    logical :: use_bounded_tile

    call load_collision_boundary_config(cfg)

    boundary_offsets(:) = 0.0d0
    use_bounded_tile = present(offsets)
    if (use_bounded_tile) boundary_offsets(:) = offsets(:)

    sdoms(1:2, 1:3) = dble(isdoms(1:2, 1:3))
    is_possible_to_be_covered = .false.

    ! Shorthand for the hole extent used by rectangle-/cylinder-/hyperboloid-/ellipsoid-hole.
    xl = cfg%xlrechole(1) + boundary_offsets(1)
    xu = cfg%xurechole(1) + boundary_offsets(1)
    yl = cfg%ylrechole(1) + boundary_offsets(2)
    yu = cfg%yurechole(1) + boundary_offsets(2)
    zl = cfg%zlrechole(2) + boundary_offsets(3)
    zu = cfg%zssurf + boundary_offsets(3)  ! (= zurechole(1)) - fall back to zssurf when hole params are unset.

    itype = 1
    if (cfg%boundary_type == "complex") then
      do itype = 1, cfg%nboundary_types
        call dispatch_surface_shape(cfg%boundary_types(itype), is_possible_to_be_covered)
      end do
      do itype = 1, cfg%nboundary_types
        if (cfg%boundary_types(itype) == "plane-with-circle-hole") call add_plane_with_circle_hole
      end do
    else
      call require_supported_simple_collision_boundary_type(cfg%boundary_type, "adding collision surfaces")
      if (cfg%boundary_type == "plane-with-circle-hole") then
        call add_plane_with_circle_hole
      else
        call dispatch_surface_shape(cfg%boundary_type, is_possible_to_be_covered)
      end if
    end if

    if (is_possible_to_be_covered .and. present(cover_all)) then
      if (cover_all) then
        itype = -1
        call add_cover_all
      end if
    end if

  contains

    subroutine dispatch_surface_shape(name, covered)
      character(len=*), intent(in) :: name
      logical, intent(inout) :: covered

      if (is_collision_boundary_entry_empty(name)) return
      call require_supported_collision_boundary_entry(name, "adding collision surfaces", itype)
      if (.not. is_collision_surface_boundary_type(name)) return

      select case (trim(name))
      case ("rectangle-hole"); call add_rectangle_hole_surface; covered = .true.
      case ("dent-rectangle-hole"); call add_dent_rectangle_hole_surface; covered = .true.
      case ("cylinder-hole"); call add_cylinder_hole_surface; covered = .true.
      case ("hyperboloid-hole"); call add_hyperboloid_hole_surface; covered = .true.
      case ("ellipsoid-hole"); call add_ellipsoid_hole_surface; covered = .true.
      case ("flat-surface"); call add_flat_surface; covered = .true.
      case ("plane-with-circle-hole")
        return
      case default
        call stop_unimplemented_surface_type(name)
      end select
    end subroutine

    subroutine stop_unimplemented_surface_type(name)
      character(len=*), intent(in) :: name

      write (error_unit, '(2a)') "Error.Internal: collision surface boundary type is registered but not implemented: ", &
        trim(name)
      error stop 1
    end subroutine stop_unimplemented_surface_type

    !
    ! ---- Shape-local helpers: allocate → hand to try_add_boundary. ----
    ! Each variant specialises the factory new_rectangleX/Y/Z/Rectangle
    ! so the caller writes one line per face.
    !

    subroutine add_rect_X(origin, w1, w2)
      real(kind=dp), intent(in) :: origin(3), w1, w2
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect
      call require_positive_width(w1, "rectangle first width")
      call require_positive_width(w2, "rectangle second width")
      allocate (prect); prect = new_rectangleX(origin, w1, w2); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_rectangle_xyz(1_ip, origin, w1, w2))
    end subroutine

    subroutine add_rect_Y(origin, w1, w2)
      real(kind=dp), intent(in) :: origin(3), w1, w2
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect
      call require_positive_width(w1, "rectangle first width")
      call require_positive_width(w2, "rectangle second width")
      allocate (prect); prect = new_rectangleY(origin, w1, w2); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_rectangle_xyz(2_ip, origin, w1, w2))
    end subroutine

    subroutine add_rect_Z(origin, w1, w2)
      real(kind=dp), intent(in) :: origin(3), w1, w2
      class(t_Boundary), pointer :: pbound
      type(t_RectangleXYZ), pointer :: prect
      call require_positive_width(w1, "rectangle first width")
      call require_positive_width(w2, "rectangle second width")
      allocate (prect); prect = new_rectangleZ(origin, w1, w2); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_rectangle_xyz(3_ip, origin, w1, w2))
    end subroutine

    !> Add an oblique rectangle defined by four corner points (3x4 matrix).
    subroutine add_rect_quad(corners)
      real(kind=dp), intent(in) :: corners(3, 4)
      class(t_Boundary), pointer :: pbound
      type(t_Rectangle), pointer :: prect
      allocate (prect); prect = new_Rectangle(corners); pbound => prect
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_quad(corners))
    end subroutine

    subroutine add_circle_Z(origin, radius)
      real(kind=dp), intent(in) :: origin(3), radius
      class(t_Boundary), pointer :: pbound
      type(t_CircleXYZ), pointer :: pcircle
      call require_positive_width(radius, "circle radius")
      allocate (pcircle); pcircle = new_CircleZ(origin, radius); pbound => pcircle
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_circle_xyz(3_ip, origin, radius))
    end subroutine

    subroutine add_cylinder_Z(origin, radius, height)
      real(kind=dp), intent(in) :: origin(3), radius, height
      class(t_Boundary), pointer :: pbound
      type(t_CylinderXYZ), pointer :: pcylinder
      call require_positive_width(radius, "cylinder radius")
      call require_positive_width(height, "cylinder height")
      allocate (pcylinder); pcylinder = new_cylinderZ(origin, radius, height); pbound => pcylinder
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_cylinder_xyz(3_ip, origin, radius, height))
    end subroutine

    subroutine add_plane_circle_hole_Z(origin, radius)
      real(kind=dp), intent(in) :: origin(3), radius
      class(t_Boundary), pointer :: pbound
      type(t_BoundedPlaneZWithCircleHole), pointer :: pbounded_plane
      type(t_PlaneXYZWithCircleHole), pointer :: pplane
      call require_positive_width(radius, "plane-with-circle-hole radius")
      if (use_bounded_tile) then
        call require_positive_width(dble(cfg%nx), "periodic tile x width")
        call require_positive_width(dble(cfg%ny), "periodic tile y width")
        allocate (pbounded_plane)
        pbounded_plane = new_bounded_plane_z_with_circle_hole(origin, radius, &
                                                              boundary_offsets(1), boundary_offsets(1) + cfg%nx, &
                                                              boundary_offsets(2), boundary_offsets(2) + cfg%ny)
        pbound => pbounded_plane
        call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                              aabb_from_bounded_axis_plane(3_ip, origin, boundary_offsets(1), &
                                                           boundary_offsets(1) + cfg%nx, boundary_offsets(2), &
                                                           boundary_offsets(2) + cfg%ny))
        return
      else
        allocate (pplane); pplane = new_planeXYZWithCircleHoleZ(origin, radius); pbound => pplane
      end if
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, aabb_from_axis_plane(3_ip, origin(3)))
    end subroutine

    subroutine require_hole_extent(shape_name)
      character(len=*), intent(in) :: shape_name

      call require_positive_width(xu - xl, trim(shape_name)//" x width")
      call require_positive_width(yu - yl, trim(shape_name)//" y width")
      call require_positive_width(zu - zl, trim(shape_name)//" height")
    end subroutine require_hole_extent

    subroutine require_positive_width(value, name)
      real(kind=dp), intent(in) :: value
      character(len=*), intent(in) :: name

      if (value <= 0.0d0) call stop_invalid_surface_parameter(trim(name)//" must be positive")
    end subroutine require_positive_width

    subroutine stop_invalid_surface_parameter(message)
      character(len=*), intent(in) :: message

      write (error_unit, '(a,a,a,i0,a)') "Error.Input: ", trim(message), " at boundary_types(", itype, ")"
      error stop 1
    end subroutine stop_invalid_surface_parameter

    !
    ! ---- Shape constructors ----
    !

    subroutine add_flat_surface
      class(t_Boundary), pointer :: pbound
      type(t_PlaneXYZ), pointer :: pplane
      allocate (pplane); pplane = new_planeZ(cfg%zssurf + boundary_offsets(3)); pbound => pplane
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_axis_plane(3_ip, cfg%zssurf + boundary_offsets(3)))
    end subroutine

    !  z^                   z^                 y^ _____________
    !    |        |          |        |         |      8      |
    !    |        |          |        |         |___ _____ ___|
    !  1 |        | 3      2 |        | 4       | 6 |     | 7 |
    !    |        |          |        |         |___|_____|___|
    !    |        |          |        |         |      5      |
    !    |________| -> x     |________| -> y    |_____________| -> x
    !        0
    subroutine add_rectangle_hole_surface
      real(kind=dp) :: wx, wy, wz
      real(kind=dp) :: xmin, ymin, xmax, ymax
      real(kind=dp) :: cap_padding

      call require_hole_extent("rectangle-hole")
      wx = xu - xl; wy = yu - yl; wz = zu - zl

      ! 0. Bottom surface
      call add_rect_Z([xl, yl, zl], wx, wy)

      ! Wall surfaces 1..4
      call add_rect_X([xl, yl, zl], wy, wz)
      call add_rect_Y([xl, yl, zl], wz, wx)
      call add_rect_X([xu, yl, zl], wy, wz)
      call add_rect_Y([xl, yu, zl], wz, wx)

      ! Capping plane at zssurf (large padding so particles leaving the
      ! domain still collide before they exit through the open face).
      cap_padding = PLANE_PADDING
      if (use_bounded_tile) cap_padding = 0.0d0
      xmin = -cap_padding + boundary_offsets(1)
      ymin = -cap_padding + boundary_offsets(2)
      xmax = cfg%nx + cap_padding + boundary_offsets(1)
      ymax = cfg%ny + cap_padding + boundary_offsets(2)

      call add_rect_Z([xmin, ymin, zu], xmax - xmin, yl - ymin)
      call add_rect_Z([xmin, yl, zu], xl - xmin, wy)
      call add_rect_Z([xu, yl, zu], xmax - xu, wy)
      call add_rect_Z([xmin, yu, zu], xmax - xmin, ymax - yu)
    end subroutine

    ! z^                   z^                 y^ _____________
    !    \        /          \        /         |      12     |
    !   1 \      / 3       2  \      / 4        |___ _____ ___|
    !      \    /              \    /           | 10|     | 11|
    !      /    \              /    \           |___|_____|___|
    !   5 /      \ 7       6  /      \  8       |      9      |
    !    /________\ -> x     /________\ -> y    |_____________| -> x
    !        0
    subroutine add_dent_rectangle_hole_surface
      real(kind=dp) :: wx, wy, wz
      real(kind=dp) :: xmin, ymin, xmax, ymax
      real(kind=dp) :: cap_padding
      real(kind=dp) :: xmidl, xmidu, ymidl, ymidu, zmid

      call require_hole_extent("dent-rectangle-hole")
      call require_positive_width(cfg%rcurv, "rcurv")
      wx = xu - xl; wy = yu - yl; wz = zu - zl
      xmidl = 0.5d0*(xl + xu) - 0.5d0*wx*cfg%rcurv
      xmidu = xmidl + wx*cfg%rcurv
      ymidl = 0.5d0*(yl + yu) - 0.5d0*wy*cfg%rcurv
      ymidu = ymidl + wy*cfg%rcurv
      zmid = zl + 0.5d0*wz

      ! 0. Bottom surface
      call add_rect_Z([xl, yl, zl], wx, wy)

      ! Upper sloped walls 1..4 (zu face)
      call add_rect_quad(reshape([[xl, yu, zu], [xmidl, ymidu, zmid], &
                                  [xmidl, ymidl, zmid], [xl, yl, zu]], [3, 4]))
      call add_rect_quad(reshape([[xl, yl, zu], [xmidl, ymidl, zmid], &
                                  [xmidu, ymidl, zmid], [xu, yl, zu]], [3, 4]))
      call add_rect_quad(reshape([[xu, yl, zu], [xmidu, ymidl, zmid], &
                                  [xmidu, ymidu, zmid], [xu, yu, zu]], [3, 4]))
      call add_rect_quad(reshape([[xu, yu, zu], [xmidu, ymidu, zmid], &
                                  [xmidl, ymidu, zmid], [xl, yu, zu]], [3, 4]))

      ! Lower sloped walls 5..8 (zl face)
      call add_rect_quad(reshape([[xl, yu, zl], [xmidl, ymidu, zmid], &
                                  [xmidl, ymidl, zmid], [xl, yl, zl]], [3, 4]))
      call add_rect_quad(reshape([[xl, yl, zl], [xmidl, ymidl, zmid], &
                                  [xmidu, ymidl, zmid], [xu, yl, zl]], [3, 4]))
      call add_rect_quad(reshape([[xu, yl, zl], [xmidu, ymidl, zmid], &
                                  [xmidu, ymidu, zmid], [xu, yu, zl]], [3, 4]))
      call add_rect_quad(reshape([[xu, yu, zl], [xmidu, ymidu, zmid], &
                                  [xmidl, ymidu, zmid], [xl, yu, zl]], [3, 4]))

      ! Capping plane at zssurf (large padding).
      cap_padding = PLANE_PADDING
      if (use_bounded_tile) cap_padding = 0.0d0
      xmin = -cap_padding + boundary_offsets(1)
      ymin = -cap_padding + boundary_offsets(2)
      xmax = cfg%nx + cap_padding + boundary_offsets(1)
      ymax = cfg%ny + cap_padding + boundary_offsets(2)

      call add_rect_Z([xmin, ymin, zu], xmax - xmin, yl - ymin)
      call add_rect_Z([xmin, yl, zu], xl - xmin, wy)
      call add_rect_Z([xu, yl, zu], xmax - xu, wy)
      call add_rect_Z([xmin, yu, zu], xmax - xmin, ymax - yu)
    end subroutine

    subroutine add_cylinder_hole_surface
      real(kind=dp) :: origin(3), origin_bottom(3)
      real(kind=dp) :: radius, height

      call require_hole_extent("cylinder-hole")
      origin(1:3) = [0.5d0*(xl + xu), 0.5d0*(yl + yu), zu]
      radius = 0.5d0*(xu - xl)
      height = zu - zl
      origin_bottom(1:2) = origin(1:2)
      origin_bottom(3) = origin(3) - height

      call add_plane_circle_hole_Z(origin, radius)
      call add_cylinder_Z(origin_bottom, radius, height)

      origin_bottom = [origin(1), origin(2), zl]
      call add_circle_Z(origin_bottom, radius)
    end subroutine

    subroutine add_hyperboloid_hole_surface
      class(t_Boundary), pointer :: pbound
      type(t_HyperboloidXYZ), pointer :: phyperboloid
      real(kind=dp) :: origin(3), origin_hyperboloid(3), origin_bottom(3)
      real(kind=dp) :: max_radius, min_radius, height

      call require_hole_extent("hyperboloid-hole")
      max_radius = 0.5d0*(xu - xl)
      min_radius = cfg%rcurv*max_radius
      height = zu - zl
      call require_positive_width(min_radius, "hyperboloid-hole min radius")

      origin(:) = [0.5d0*(xl + xu), 0.5d0*(yl + yu), zu]
      origin_hyperboloid(:) = [origin(1), origin(2), zu - 0.5d0*height]
      origin_bottom(:) = [origin(1), origin(2), zl]

      call add_plane_circle_hole_Z(origin, max_radius)

      allocate (phyperboloid)
      phyperboloid = new_HyperboloidZ(origin_hyperboloid, max_radius, min_radius, height)
      pbound => phyperboloid
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_centered_profile(3_ip, origin_hyperboloid, max_radius, height))

      call add_circle_Z(origin_bottom, max_radius)
    end subroutine

    subroutine add_ellipsoid_hole_surface
      class(t_Boundary), pointer :: pbound
      type(t_EllipsoidXYZ), pointer :: pellipsoid
      real(kind=dp) :: origin(3), origin_ellipsoid(3), origin_bottom(3)
      real(kind=dp) :: max_radius, min_radius, height

      call require_hole_extent("ellipsoid-hole")
      min_radius = 0.5d0*(xu - xl)
      max_radius = cfg%rcurv*min_radius
      height = zu - zl
      call require_positive_width(max_radius, "ellipsoid-hole max radius")

      origin(:) = [0.5d0*(xl + xu), 0.5d0*(yl + yu), zu]
      origin_ellipsoid(:) = [origin(1), origin(2), zu - 0.5d0*height]
      origin_bottom(:) = [origin(1), origin(2), zl]

      call add_plane_circle_hole_Z(origin, min_radius)

      allocate (pellipsoid)
      pellipsoid = new_EllipsoidZ(origin_ellipsoid, min_radius, max_radius, height)
      pbound => pellipsoid
      call try_add_boundary(boundaries, pbound, itype, sdoms, accelerator, &
                            aabb_from_centered_profile(3_ip, origin_ellipsoid, max_radius, height))

      call add_circle_Z(origin_bottom, min_radius)
    end subroutine

    subroutine add_cover_all
      real(kind=dp) :: rnx, rny
      integer(kind=ip) :: iboundary, nboundaries_init

      nboundaries_init = boundaries%nboundaries
      rnx = cfg%nx; rny = cfg%ny

      call add_rect_Z([boundary_offsets(1), boundary_offsets(2), boundary_offsets(3)], rnx, rny)

      ! Tag freshly added boundaries to indicate no photoelectron emission.
      do iboundary = nboundaries_init + 1, boundaries%nboundaries
        boundaries%boundaries(iboundary)%ref%material%tag = -1
      end do
    end subroutine

    subroutine add_plane_with_circle_hole
      real(kind=dp) :: origin(3), origin_bottom(3)
      real(kind=dp) :: radius, height

      height = cfg%plane_with_circle_hole_height(itype)
      origin_bottom(1:3) = [0.5d0*cfg%nx + boundary_offsets(1), &
                            0.5d0*cfg%ny + boundary_offsets(2), &
                            cfg%plane_with_circle_hole_zlower(itype) + boundary_offsets(3)]
      origin(1:2) = origin_bottom(1:2)
      origin(3) = origin_bottom(3) + height
      radius = cfg%plane_with_circle_hole_radius(itype)
      call require_positive_width(height, "plane_with_circle_hole_height")
      call require_positive_width(radius, "plane_with_circle_hole_radius")

      call add_plane_circle_hole_Z(origin, radius)
      call add_cylinder_Z(origin_bottom, radius, height)
      call add_plane_circle_hole_Z(origin_bottom, radius)
    end subroutine

  end subroutine add_surfaces

end module m_surfaces

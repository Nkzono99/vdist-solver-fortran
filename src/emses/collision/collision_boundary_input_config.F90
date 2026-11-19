!> vdist adapter for the EMSES collision-geometry input view.
!>
!> The neighbouring objects/surfaces/boundary_* /collision_accelerator files are
!> verbatim copies of MPIEMSES3D `src/physics/collision/`. EMSES fills this type
!> from `m_input_parameters`; vdist fills it from the namelist variables in
!> `m_allcom`. The type layout follows MPIEMSES3D
!> `src/main/config/collision_boundary_input_config.F90`.
!> See docs/emses_sync.md before editing.
module m_collision_boundary_input_config
  use m_allcom, only: iboundary_type => nboundary_types
  use m_kinds, only: dp, ip
  implicit none
  private

  public :: t_collision_boundary_config
  public :: load_collision_boundary_config

  type :: t_collision_boundary_config
    integer(kind=ip) :: nboundary_types = iboundary_type
    integer(kind=ip) :: nx = 0
    integer(kind=ip) :: ny = 0
    integer(kind=ip) :: nz = 0
    character(len=30) :: boundary_type = ""
    character(len=30) :: boundary_types(iboundary_type) = ""
    real(kind=dp) :: rectangle_shape(6, iboundary_type) = 0.0d0
    real(kind=dp) :: rectangle_corners(3, 4, iboundary_type) = 0.0d0
    real(kind=dp) :: sphere_origin(3, iboundary_type) = 0.0d0
    real(kind=dp) :: sphere_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: circle_origin(3, iboundary_type) = 0.0d0
    real(kind=dp) :: circle_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: cuboid_shape(6, iboundary_type) = 0.0d0
    real(kind=dp) :: cylinder_origin(3, iboundary_type) = 0.0d0
    real(kind=dp) :: cylinder_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: cylinder_height(iboundary_type) = 0.0d0
    real(kind=dp) :: disk_origin(3, iboundary_type) = 0.0d0
    real(kind=dp) :: disk_height(iboundary_type) = 0.0d0
    real(kind=dp) :: disk_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: disk_inner_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: plane_with_circle_origin(3, iboundary_type) = 0.0d0
    real(kind=dp) :: plane_with_circle_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: plane_with_circle_hole_zlower(iboundary_type) = 0.0d0
    real(kind=dp) :: plane_with_circle_hole_height(iboundary_type) = 0.0d0
    real(kind=dp) :: plane_with_circle_hole_radius(iboundary_type) = 0.0d0
    real(kind=dp) :: rcurv = 0.0d0
    real(kind=dp) :: xlrechole(2) = 0.0d0
    real(kind=dp) :: xurechole(2) = 0.0d0
    real(kind=dp) :: ylrechole(2) = 0.0d0
    real(kind=dp) :: yurechole(2) = 0.0d0
    real(kind=dp) :: zlrechole(2) = 0.0d0
    real(kind=dp) :: zurechole(2) = 0.0d0
    real(kind=dp) :: zssurf = 0.0d0
  end type t_collision_boundary_config

contains

  subroutine load_collision_boundary_config(cfg)
    use m_allcom, only: boundary_type, boundary_types, circle_origin, circle_radius, cuboid_shape, &
                        cylinder_height, cylinder_origin, cylinder_radius, disk_height, &
                        disk_inner_radius, disk_origin, disk_radius, nboundary_types, nx, ny, nz, &
                        plane_with_circle_hole_height, plane_with_circle_hole_radius, &
                        plane_with_circle_hole_zlower, plane_with_circle_origin, plane_with_circle_radius, &
                        rcurv, rectangle_corners, rectangle_shape, sphere_origin, sphere_radius, &
                        xlrechole, xurechole, ylrechole, yurechole, zlrechole, zssurf, zurechole
    type(t_collision_boundary_config), intent(out) :: cfg

    cfg%nboundary_types = nboundary_types
    cfg%nx = nx
    cfg%ny = ny
    cfg%nz = nz
    cfg%boundary_type = boundary_type
    cfg%boundary_types = boundary_types
    cfg%rectangle_shape = rectangle_shape
    cfg%rectangle_corners = rectangle_corners
    cfg%sphere_origin = sphere_origin
    cfg%sphere_radius = sphere_radius
    cfg%circle_origin = circle_origin
    cfg%circle_radius = circle_radius
    cfg%cuboid_shape = cuboid_shape
    cfg%cylinder_origin = cylinder_origin
    cfg%cylinder_radius = cylinder_radius
    cfg%cylinder_height = cylinder_height
    cfg%disk_origin = disk_origin
    cfg%disk_height = disk_height
    cfg%disk_radius = disk_radius
    cfg%disk_inner_radius = disk_inner_radius
    cfg%plane_with_circle_origin = plane_with_circle_origin
    cfg%plane_with_circle_radius = plane_with_circle_radius
    cfg%plane_with_circle_hole_zlower = plane_with_circle_hole_zlower
    cfg%plane_with_circle_hole_height = plane_with_circle_hole_height
    cfg%plane_with_circle_hole_radius = plane_with_circle_hole_radius
    cfg%rcurv = rcurv
    cfg%xlrechole = xlrechole
    cfg%xurechole = xurechole
    cfg%ylrechole = ylrechole
    cfg%yurechole = yurechole
    cfg%zlrechole = zlrechole
    cfg%zurechole = zurechole
    cfg%zssurf = zssurf
  end subroutine load_collision_boundary_config

end module m_collision_boundary_input_config

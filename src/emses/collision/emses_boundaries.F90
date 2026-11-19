!> Manage internal boundaries.
!>
!> Namelist Parameters.
!> &ptcond
!>   boundary_type = 'none'|
!>                   'flat-surface'|
!>                   'rectangle-hole'|'cylinder-hole'|'hyperboloid-hole'|'ellipsoid-hole'|
!>                   'rectangle[xyz]'|'circle[x/y/z]'|'cuboid'|'disk[x/y/z]
!>                   'complex'
!>
!>   ! Use if boundary_type is '****-surface' or '****-hole'.
!>   zssurf = Surface Height [grid]
!>
!>   ! Use if boundary_type is '****-hole'.
!>   [x/y/z][l/u]pc = Hole [lower/upper] limit grid（[x/y/z] coordinate) [grid]
!>
!>   ! Use if boundary_type is 'complex'
!>   boundary_types(ntypes) = <boundary_type>|
!>
!>   ! Use if boundary_types(itype) is 'rectangle'
!>   rectangle_shape(ntypes, 6) = Rectangle location (xmin, xmax, ymin, ymax, zmin, zmax)
!>
!>   ! Use if boundary_types is 'circle[x/y/z]'
!>   circle_origin(ntypes, 3) = Circle center coordinates
!>   circle_radius(ntypes) = Circle radius
!>
!>   ! Use if boundary_types(itype) is 'cuboid'
!>   cuboid_shape(ntypes, 6) = Cuboid location (xmin, xmax, ymin, ymax, zmin, zmax)
!>
!>   ! Use if boundary_types(itype) is 'disk[x/y/z]
!>   disk_origin(ntypes, 3) = Disk center bottom coordinates
!>   disk_height(ntypes) = Disk height (= thickness)
!>   disk_radius(ntypes) = Disk outer radius
!>   disk_inner_radius(ntypes) = Disk inner radius
!>
!>   ! Rotation angle of all boundaries [deg].
!>   boundary_rotation_deg(3) = 0d0, 0d0, 0d0
!> &
!>
module m_emses_boundaries
    !! Build the EMSES internal collision geometry for vdist.
    !!
    !! `m_surfaces` and `m_objects` are verbatim copies of MPIEMSES3D; this
    !! module only adds the vdist probability tag (docs/emses_sync.md).
    use finbound
    use m_surfaces, only: add_surfaces
    use m_objects, only: add_objects
    use m_kinds, only: ip, dp
    implicit none

    private
    public create_simple_collision_boundaries

contains

    function create_simple_collision_boundaries(isdoms, cover_all, tag) result(boundaries)
        integer(kind=ip), intent(in) :: isdoms(2, 3)
        logical, intent(in), optional :: cover_all
        integer(kind=ip), intent(in), optional :: tag
        type(t_BoundaryList) :: boundaries

        integer(kind=ip) :: iboundary

        boundaries = new_BoundaryList()

        call add_surfaces(boundaries, isdoms, cover_all=cover_all)

        call add_objects(boundaries, isdoms, (/0d0, 0d0, 0d0/))

        if (present(tag)) then
            do iboundary = 1, boundaries%nboundaries
                boundaries%boundaries(iboundary)%ref%material%tag = tag
            end do
        end if
    end function

end module

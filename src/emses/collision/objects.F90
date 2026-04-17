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
module m_objects
    use finbound
    use m_allcom, only: nx, ny, nz, &
                      boundary_type, nboundary_types, boundary_types, &
                      cylinder_origin, cylinder_radius, cylinder_height, &
                      rectangle_shape, &
                      sphere_origin, sphere_radius, &
                      circle_origin, circle_radius, &
                      cuboid_shape, &
                      disk_origin, disk_height, disk_radius, disk_inner_radius, &
                      plane_with_circle_radius, plane_with_circle_origin
    use m_kinds, only: ip, lp, sp, dp
    implicit none

    real(kind=dp), parameter :: extent(2, 3) = &
                                reshape([[1.0d0, 1.0d0], [1.0d0, 1.0d0], [1.0d0, 1.0d0]]*2, [2, 3])

    private
    public add_objects

contains

    subroutine add_objects(boundaries, isdoms, offsets)
        type(t_BoundaryList), intent(inout) :: boundaries
        integer(kind=ip), intent(in) :: isdoms(2, 3)
        real(kind=dp), intent(in) :: offsets(3)

        real(kind=dp) :: xl, yl, zl
        real(kind=dp) :: xu, yu, zu

        real(kind=dp) :: sdoms(2, 3)
        integer(kind=ip) :: itype

        if (boundary_type /= "complex") then
            return
        end if

        sdoms(1:2, 1:3) = dble(isdoms(1:2, 1:3))

        if (boundary_type == "complex") then
            do itype = 1, nboundary_types
                if (boundary_types(itype) == 'rectangle') then
                    call add_rectangle
                else if (boundary_types(itype) == 'sphere') then
                    call add_sphere
                else if (boundary_types(itype) == 'circlex') then
                    call add_circleXYZ(1)
                else if (boundary_types(itype) == 'circley') then
                    call add_circleXYZ(2)
                else if (boundary_types(itype) == 'circlez') then
                    call add_circleXYZ(3)
                else if (boundary_types(itype) == 'cuboid') then
                    call add_cuboid
                else if (boundary_types(itype) == 'cylinderx') then
                    call add_cylinderXYZ(1)
                else if (boundary_types(itype) == 'cylindery') then
                    call add_cylinderXYZ(2)
                else if (boundary_types(itype) == 'cylinderz') then
                    call add_cylinderXYZ(3)
                else if (boundary_types(itype) == 'open-cylinderx') then
                    call add_open_cylinderXYZ(1)
                else if (boundary_types(itype) == 'open-cylindery') then
                    call add_open_cylinderXYZ(2)
                else if (boundary_types(itype) == 'open-cylinderz') then
                    call add_open_cylinderXYZ(3)
                else if (boundary_types(itype) == 'diskx') then
                    call add_disk(1)
                else if (boundary_types(itype) == 'disky') then
                    call add_disk(2)
                else if (boundary_types(itype) == 'diskz') then
                    call add_disk(3)
                else if (boundary_types(itype) == 'plane-with-circlex') then
                    call add_plane_with_circleXYZ(1)
                else if (boundary_types(itype) == 'plane-with-circley') then
                    call add_plane_with_circleXYZ(2)
                else if (boundary_types(itype) == 'plane-with-circlez') then
                    call add_plane_with_circleXYZ(3)
                end if
            end do
        end if

    contains

        subroutine add_rectangle
            real(kind=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
            class(t_Boundary), pointer :: pbound
            type(t_RectangleXYZ), pointer :: prect

            xmin = rectangle_shape(1, itype) + offsets(1)
            xmax = rectangle_shape(2, itype) + offsets(1)
            ymin = rectangle_shape(3, itype) + offsets(2)
            ymax = rectangle_shape(4, itype) + offsets(2)
            zmin = rectangle_shape(5, itype) + offsets(3)
            zmax = rectangle_shape(6, itype) + offsets(3)

            allocate (prect)
            if (xmin == xmax) then
                prect = new_rectangleX([xmin, ymin, zmin], ymax - ymin, zmax - zmin)
            else if (ymin == ymax) then
                prect = new_rectangleY([xmin, ymin, zmin], zmax - zmin, xmax - xmin)
            else if (zmin == zmax) then
                prect = new_rectangleZ([xmin, ymin, zmin], xmax - xmin, ymax - ymin)
            end if
            pbound => prect

            pbound%material%tag = itype

            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if
        end subroutine

        subroutine add_circleXYZ(axis)
            integer(kind=ip), intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CircleXYZ), pointer :: pcircle

            allocate (pcircle)
            pcircle = new_CircleXYZ(axis, circle_origin(:, itype) + offsets(:), circle_radius(itype))
            pbound => pcircle

            pbound%material%tag = itype

            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if
        end subroutine

        !             ------------- (max)
        !           / |           /
        !          /  |    6     / |
        !         /   |      5  /  |
        !    z ^  -------------  4 |
        !      | | 1   --------|---/  ^ y
        !        |   /  2      |  /  /
        !        |  /      3   | /
        !        | /           |/
        !  (min)  -------------/
        !                      -> x
        subroutine add_cuboid
            class(t_Boundary), pointer :: pbound
            type(t_RectangleXYZ), pointer :: prect
            real(kind=dp) :: xmin, xmax, ymin, ymax, zmin, zmax
            real(kind=dp) :: wx, wy, wz

            xmin = cuboid_shape(1, itype) + offsets(1)
            xmax = cuboid_shape(2, itype) + offsets(1)
            ymin = cuboid_shape(3, itype) + offsets(2)
            ymax = cuboid_shape(4, itype) + offsets(2)
            zmin = cuboid_shape(5, itype) + offsets(3)
            zmax = cuboid_shape(6, itype) + offsets(3)

            wx = xmax - xmin
            wy = ymax - ymin
            wz = zmax - zmin

            ! 1.
            allocate (prect)
            prect = new_rectangleX([xmin, ymin, zmin], wy, wz)
            pbound => prect
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 2.
            allocate (prect)
            prect = new_rectangleY([xmin, ymin, zmin], wz, wx)
            pbound => prect
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 3.
            allocate (prect)
            prect = new_rectangleZ([xmin, ymin, zmin], wx, wy)
            pbound => prect
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 4.
            allocate (prect)
            prect = new_rectangleX([xmax, ymin, zmin], wy, wz)
            pbound => prect
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 5.
            allocate (prect)
            prect = new_rectangleY([xmin, ymax, zmin], wz, wx)
            pbound => prect
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 6.
            allocate (prect)
            prect = new_rectangleZ([xmin, ymin, zmax], wx, wy)
            pbound => prect
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if
        end subroutine

        subroutine add_sphere
            real(kind=dp) :: origin(3), radius
            class(t_Boundary), pointer :: pbound
            type(t_Sphere), pointer :: pshere

            origin(:) = sphere_origin(:, itype) + offsets(:)
            radius = sphere_radius(itype)

            allocate (pshere)
            pshere = new_Sphere(origin, radius)
            pbound => pshere

            pbound%material%tag = itype

            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pshere)
            end if
        end subroutine

        subroutine add_open_cylinderXYZ(axis)
            integer(kind=ip), intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CylinderXYZ), pointer :: pcylinder

            real(kind=dp) :: height
            real(kind=dp) :: lower_origin(3)
            real(kind=dp) :: upper_origin(3)
            real(kind=dp) :: radius

            height = cylinder_height(itype)
            lower_origin(1:3) = cylinder_origin(1:3, itype) + offsets(:)
            upper_origin(1:3) = cylinder_origin(1:3, itype) + offsets(:)
            upper_origin(axis) = cylinder_origin(axis, itype) + height
            radius = cylinder_radius(itype)

            ! Outer cylinder
            allocate (pcylinder)
            pcylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
            pbound => pcylinder
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcylinder)
            end if
        end subroutine

        subroutine add_cylinderXYZ(axis)
            integer(kind=ip), intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CylinderXYZ), pointer :: pcylinder
            type(t_CircleXYZ), pointer :: pcircle

            real(kind=dp) :: height
            real(kind=dp) :: lower_origin(3)
            real(kind=dp) :: upper_origin(3)
            real(kind=dp) :: radius

            height = cylinder_height(itype)
            lower_origin(1:3) = cylinder_origin(1:3, itype) + offsets(:)
            upper_origin(1:3) = cylinder_origin(1:3, itype) + offsets(:)
            upper_origin(axis) = cylinder_origin(axis, itype) + height
            radius = cylinder_radius(itype)

            ! Outer cylinder
            allocate (pcylinder)
            pcylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
            pbound => pcylinder
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcylinder)
            end if

            ! Lower circle
            allocate (pcircle)
            pcircle = new_CircleXYZ(axis, lower_origin(:), radius)
            pbound => pcircle
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if

            ! Upper circle
            allocate (pcircle)
            pcircle = new_CircleXYZ(axis, upper_origin(:), radius)
            pbound => pcircle
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if
        end subroutine

        subroutine add_plane_with_circleXYZ(axis)
            integer(kind=ip), intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZWithCircleHole), pointer :: pplane
            type(t_CylinderXYZ), pointer :: pcylinder

            real(kind=dp) :: origin(3), origin_bottom(3)
            real(kind=dp) :: radius, height

            origin(:) = plane_with_circle_origin(:, itype)
            radius = plane_with_circle_radius(itype)

            ! Upper surface
            allocate (pplane)
            if (axis == 1) then
                pplane = new_planeXYZWithCircleHoleX(origin, radius)
            else if (axis == 2) then
                pplane = new_planeXYZWithCircleHoleY(origin, radius)
            else if (axis == 3) then
                pplane = new_planeXYZWithCircleHoleZ(origin, radius)
            end if
            pbound => pplane
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pplane)
            end if
        end subroutine

        subroutine add_disk(axis)
            integer(kind=ip), intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CylinderXYZ), pointer :: pouter_cylinder
            type(t_CylinderXYZ), pointer :: pinner_cylinder
            type(t_DonutXYZ), pointer :: plower_donut
            type(t_DonutXYZ), pointer :: pupper_donut

            real(kind=dp) :: height
            real(kind=dp) :: lower_origin(3)
            real(kind=dp) :: upper_origin(3)
            real(kind=dp) :: radius
            real(kind=dp) :: inner_radius

            height = disk_height(itype)
            lower_origin(1:3) = disk_origin(1:3, itype) + offsets(:)
            upper_origin(1:3) = disk_origin(1:3, itype) + offsets(:)
            upper_origin(axis) = upper_origin(axis) + height
            radius = disk_radius(itype)
            inner_radius = disk_inner_radius(itype)

            ! Outer cylinder
            allocate (pouter_cylinder)
            pouter_cylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
            pbound => pouter_cylinder
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pouter_cylinder)
            end if

            ! Inner cylinder
            allocate (pinner_cylinder)
            pinner_cylinder = new_cylinderXYZ(axis, lower_origin, inner_radius, height)
            pbound => pinner_cylinder
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pinner_cylinder)
            end if

            ! Lower donut
            allocate (plower_donut)
            plower_donut = new_DonutXYZ(axis, lower_origin, inner_radius, height)
            pbound => plower_donut
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (plower_donut)
            end if

            ! Upper donut
            allocate (pupper_donut)
            pupper_donut = new_DonutXYZ(axis, upper_origin, inner_radius, height)
            pbound => pupper_donut
            pbound%material%tag = itype
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pupper_donut)
            end if
        end subroutine

    end subroutine

end module

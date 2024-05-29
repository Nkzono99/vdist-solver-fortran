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
module emses_boundaries
    use finbound
    use allcom, only: xlrechole, ylrechole, zlrechole, &
                      xurechole, yurechole, zurechole, &
                      zssurf, &
                      nx, ny, nz, &
                      boundary_type, nboundary_types, boundary_types, &
                      cylinder_origin, cylinder_radius, cylinder_height, &
                      rcurv, &
                      rectangle_shape, &
                      sphere_origin, sphere_radius, &
                      circle_origin, circle_radius, &
                      cuboid_shape, &
                      disk_origin, disk_height, disk_radius, disk_inner_radius
    implicit none

    double precision, parameter :: extent(2, 3) = &
        reshape([[1.0d0, 1.0d0], [1.0d0, 1.0d0], [1.0d0, 1.0d0]]*2, [2, 3])

    private
    public create_simple_collision_boundaries

contains

    function create_simple_collision_boundaries(isdoms, cover_all, tag) result(boundaries)
        integer, intent(in) :: isdoms(2, 3)
        logical, intent(in), optional :: cover_all
        integer, intent(in), optional :: tag
        type(t_BoundaryList) :: boundaries

        double precision :: xl, yl, zl
        double precision :: xu, yu, zu

        double precision :: sdoms(2, 3)
        integer :: itype

        logical :: is_possible_to_be_covered

        sdoms(1:2, 1:3) = dble(isdoms(1:2, 1:3))
        is_possible_to_be_covered = .false.

        ! simplified parameter names
        xl = xlrechole(1)
        xu = xurechole(1)
        yl = ylrechole(1)
        yu = yurechole(1)
        zl = zlrechole(2)
        zu = zssurf ! (= zurechole(1)) ! Use 'zssurf' instead of 'zurechole(1)' in case the hole parameter is not used.

        boundaries = new_BoundaryList()
        if (boundary_type == "rectangle-hole") then
            call add_rectangle_hole_surface
            is_possible_to_be_covered = .true.
        else if (boundary_type == "dent-rectangle-hole") then
            call add_dent_rectangle_hole_surface
            is_possible_to_be_covered = .true.
        else if (boundary_type == "cylinder-hole") then
            call add_cylinder_hole_surface
            is_possible_to_be_covered = .true.
        else if (boundary_type == "hyperboloid-hole") then
            call add_hyperboloid_hole_surface
            is_possible_to_be_covered = .true.
        else if (boundary_type == "ellipsoid-hole") then
            call add_ellipsoid_hole_surface
            is_possible_to_be_covered = .true.
        else if (boundary_type == "flat-surface") then
            call add_flat_surface
            is_possible_to_be_covered = .true.
        else if (boundary_type == "complex") then
            do itype = 1, nboundary_types
                if (boundary_types(itype) == "rectangle-hole") then
                    call add_rectangle_hole_surface
                    is_possible_to_be_covered = .true.
                else if (boundary_types(itype) == "dent-rectangle-hole") then
                    call add_dent_rectangle_hole_surface
                    is_possible_to_be_covered = .true.
                else if (boundary_types(itype) == "cylinder-hole") then
                    call add_cylinder_hole_surface
                    is_possible_to_be_covered = .true.
                else if (boundary_types(itype) == "hyperboloid-hole") then
                    call add_hyperboloid_hole_surface
                    is_possible_to_be_covered = .true.
                else if (boundary_types(itype) == "ellipsoid-hole") then
                    call add_ellipsoid_hole_surface
                    is_possible_to_be_covered = .true.
                else if (boundary_types(itype) == "flat-surface") then
                    call add_flat_surface
                    is_possible_to_be_covered = .true.
                else if (boundary_types(itype) == 'rectangle') then
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
                else if (boundary_types(itype) == 'diskx') then
                    call add_disk(1)
                else if (boundary_types(itype) == 'disky') then
                    call add_disk(2)
                else if (boundary_types(itype) == 'diskz') then
                    call add_disk(3)
                end if
            end do
        end if

        if (is_possible_to_be_covered .and. present(cover_all) .and. cover_all) then
            call add_cover_all
        end if

        if (present(tag)) then
            block
                integer :: iboundary

                do iboundary = 1, boundaries%nboundaries
                    boundaries%boundaries(iboundary)%ref%material%tag = tag
                end do
            end block
        end if

    contains

        subroutine add_rectangle
            double precision :: xmin, xmax, ymin, ymax, zmin, zmax
            class(t_Boundary), pointer :: pbound
            type(t_RectangleXYZ), pointer :: prect

            xmin = rectangle_shape(itype, 1)
            xmax = rectangle_shape(itype, 2)
            ymin = rectangle_shape(itype, 3)
            ymax = rectangle_shape(itype, 4)
            zmin = rectangle_shape(itype, 5)
            zmax = rectangle_shape(itype, 6)

            allocate (prect)
            if (xmin == xmax) then
                prect = new_rectangleX([xmin, ymin, zmin], ymax - ymin, zmax - zmin)
            else if (ymin == ymax) then
                prect = new_rectangleY([xmin, ymin, zmin], zmax - zmin, xmax - xmin)
            else if (zmin == zmax) then
                prect = new_rectangleZ([xmin, ymin, zmin], xmax - xmin, ymax - ymin)
            end if
            pbound => prect

            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if
        end subroutine

        subroutine add_circleXYZ(axis)
            integer, intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CircleXYZ), pointer :: pcircle

            allocate (pcircle)
            pcircle = new_CircleXYZ(axis, circle_origin(itype, :), circle_radius(itype))
            pbound => pcircle

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
            double precision :: xmin, xmax, ymin, ymax, zmin, zmax
            double precision :: wx, wy, wz

            xmin = cuboid_shape(itype, 1)
            xmax = cuboid_shape(itype, 2)
            ymin = cuboid_shape(itype, 3)
            ymax = cuboid_shape(itype, 4)
            zmin = cuboid_shape(itype, 5)
            zmax = cuboid_shape(itype, 6)

            wx = xmax - xmin
            wy = ymax - ymin
            wz = zmax - zmin

            ! 1.
            allocate (prect)
            prect = new_rectangleX([xmin, ymin, zmin], wy, wz)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 2.
            allocate (prect)
            prect = new_rectangleY([xmin, ymin, zmin], wz, wx)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 3.
            allocate (prect)
            prect = new_rectangleZ([xmin, ymin, zmin], wx, wy)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 4.
            allocate (prect)
            prect = new_rectangleX([xmax, ymin, zmin], wy, wz)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 5.
            allocate (prect)
            prect = new_rectangleY([xmin, ymax, zmin], wz, wx)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 6.
            allocate (prect)
            prect = new_rectangleZ([xmin, ymin, zmax], wx, wy)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if
        end subroutine

        subroutine add_sphere
            double precision :: origin(3), radius
            class(t_Boundary), pointer :: pbound
            type(t_Sphere), pointer :: pshere

            origin(:) = sphere_origin(itype, :)
            radius = sphere_radius(itype)

            allocate (pshere)
            pshere = new_Sphere(origin, radius)
            pbound => pshere

            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pshere)
            end if
        end subroutine

        subroutine add_cylinderXYZ(axis)
            integer, intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CylinderXYZ), pointer :: pcylinder
            type(t_CircleXYZ), pointer :: pcircle

            double precision :: height
            double precision :: lower_origin(3)
            double precision :: upper_origin(3)
            double precision :: radius

            height = cylinder_height(itype)
            lower_origin(1:3) = cylinder_origin(itype, 1:3)
            upper_origin(1:3) = cylinder_origin(itype, 1:3)
            upper_origin(axis) = cylinder_origin(itype, axis) + height
            radius = cylinder_radius(itype)

            ! Outer cylinder
            allocate (pcylinder)
            pcylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
            pbound => pcylinder
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcylinder)
            end if

            ! Lower circle
            allocate (pcircle)
            pcircle = new_CircleXYZ(axis, lower_origin(:), radius)
            pbound => pcircle
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if

            ! Upper circle
            allocate (pcircle)
            pcircle = new_CircleXYZ(axis, upper_origin(:), radius)
            pbound => pcircle
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if
        end subroutine

        subroutine add_flat_surface
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZ), pointer :: pplane

            allocate (pplane)
            pplane = new_planeZ(zssurf)
            pbound => pplane

            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pplane)
            end if
        end subroutine

        ! z^                   z^                 y^ _____________
        !    |        |          |        |         |      8      |
        !    |        |          |        |         |___ _____ ___|
        !  1 |        | 3      2 |        | 4       | 6 |     | 7 |
        !    |        |          |        |         |___|_____|___|
        !    |        |          |        |         |      5      |
        !    |________| -> x     |________| -> y    |_____________| -> x
        !        0
        subroutine add_rectangle_hole_surface
            class(t_Boundary), pointer :: pbound
            type(t_RectangleXYZ), pointer :: prect

            double precision :: wx, wy, wz
            double precision :: xmin, ymin, xmax, ymax

            double precision :: origin(3)

            wx = xu - xl
            wy = yu - yl
            wz = zu - zl

            ! 0. Bottom surface
            origin = [xl, yl, zl]
            allocate (prect)
            prect = new_rectangleZ(origin, wx, wy)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! Wall surfaces
            ! 1.
            origin = [xl, yl, zl]
            allocate (prect)
            prect = new_rectangleX(origin, wy, wz)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 2.
            origin = [xl, yl, zl]
            allocate (prect)
            prect = new_rectangleY(origin, wz, wx)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 3.
            origin = [xu, yl, zl]
            allocate (prect)
            prect = new_rectangleX(origin, wy, wz)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 4.
            origin = [xl, yu, zl]
            allocate (prect)
            prect = new_rectangleY(origin, wz, wx)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! Add flat plane at zssurf(=zu1)
            ! Note: create the plane large enough to allow for particles to fly out of the simulation space
            xmin = -10
            ymin = -10
            xmax = nx + 10
            ymax = ny + 10

            ! 5.
            origin = [xmin, ymin, zu]
            allocate (prect)
            prect = new_rectangleZ(origin, xmax - xmin, yl - ymin)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 6.
            origin = [xmin, yl, zu]
            allocate (prect)
            prect = new_rectangleZ(origin, xl - xmin, wy)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 7.
            origin = [xu, yl, zu]
            allocate (prect)
            prect = new_rectangleZ(origin, xmax - xu, wy)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 8.
            origin = [xmin, yu, zu]
            allocate (prect)
            prect = new_rectangleZ(origin, xmax - xmin, ymax - yu)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if
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
            class(t_Boundary), pointer :: pbound
            type(t_RectangleXYZ), pointer :: prectxyz
            type(t_Rectangle), pointer :: prect

            double precision :: wx, wy, wz
            double precision :: xmin, ymin, xmax, ymax
            double precision :: xmidl, xmidu, ymidl, ymidu, zmid

            double precision :: origin(3)

            wx = xu - xl
            wy = yu - yl
            wz = zu - zl

            xmidl = 0.5d0*(xl + xu) - 0.5d0*wx*rcurv
            xmidu = xmidl + wx*rcurv
            ymidl = 0.5d0*(yl + yu) - 0.5d0*wy*rcurv
            ymidu = ymidl + wy*rcurv
            zmid = zl + 0.5d0*wz

            ! 0. Bottom surface
            origin = [xl, yl, zl]
            allocate (prectxyz)
            prectxyz = new_rectangleZ(origin, wx, wy)
            pbound => prectxyz
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prectxyz)
            end if

            ! Wall surfaces
            ! 1.
            allocate (prect)
            prect = new_Rectangle(reshape([[xl, yu, zu], &
                                           [xmidl, ymidu, zmid], &
                                           [xmidl, ymidl, zmid], &
                                           [xl, yl, zu]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 2.
            allocate (prect)
            prect = new_Rectangle(reshape([[xl, yl, zu], &
                                           [xmidl, ymidl, zmid], &
                                           [xmidu, ymidl, zmid], &
                                           [xu, yl, zu]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 3.
            allocate (prect)
            prect = new_Rectangle(reshape([[xu, yl, zu], &
                                           [xmidu, ymidl, zmid], &
                                           [xmidu, ymidu, zmid], &
                                           [xu, yu, zu]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 4.
            allocate (prect)
            prect = new_Rectangle(reshape([[xu, yu, zu], &
                                           [xmidu, ymidu, zmid], &
                                           [xmidl, ymidu, zmid], &
                                           [xl, yu, zu]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! Wall surfaces
            ! 5.
            allocate (prect)
            prect = new_Rectangle(reshape([[xl, yu, zl], &
                                           [xmidl, ymidu, zmid], &
                                           [xmidl, ymidl, zmid], &
                                           [xl, yl, zl]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 6.
            allocate (prect)
            prect = new_Rectangle(reshape([[xl, yl, zl], &
                                           [xmidl, ymidl, zmid], &
                                           [xmidu, ymidl, zmid], &
                                           [xu, yl, zl]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 7.
            allocate (prect)
            prect = new_Rectangle(reshape([[xu, yl, zl], &
                                           [xmidu, ymidl, zmid], &
                                           [xmidu, ymidu, zmid], &
                                           [xu, yu, zl]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 8.
            allocate (prect)
            prect = new_Rectangle(reshape([[xu, yu, zl], &
                                           [xmidu, ymidu, zmid], &
                                           [xmidl, ymidu, zmid], &
                                           [xl, yu, zl]], [3, 4]))
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! Add flat plane at zssurf(=zu1)
            ! Note: create the plane large enough to allow for particles to fly out of the simulation space
            xmin = -10
            ymin = -10
            xmax = nx + 10
            ymax = ny + 10

            ! 9.
            origin = [xmin, ymin, zu]
            allocate (prectxyz)
            prectxyz = new_rectangleZ(origin, xmax - xmin, yl - ymin)
            pbound => prectxyz
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prectxyz)
            end if

            ! 10.
            origin = [xmin, yl, zu]
            allocate (prectxyz)
            prectxyz = new_rectangleZ(origin, xl - xmin, wy)
            pbound => prectxyz
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prectxyz)
            end if

            ! 11.
            origin = [xu, yl, zu]
            allocate (prectxyz)
            prectxyz = new_rectangleZ(origin, xmax - xu, wy)
            pbound => prectxyz
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prectxyz)
            end if

            ! 12.
            origin = [xmin, yu, zu]
            allocate (prectxyz)
            prectxyz = new_rectangleZ(origin, xmax - xmin, ymax - yu)
            pbound => prectxyz
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prectxyz)
            end if
        end subroutine

        subroutine add_cylinder_hole_surface
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZWithCircleHole), pointer :: pplane
            type(t_CylinderXYZ), pointer :: pcylinder
            type(t_CircleXYZ), pointer :: pcircle

            double precision :: origin(3), origin_bottom(3)
            double precision :: radius, height

            origin(1:3) = [0.5d0*(xl + xu), 0.5d0*(yl + yu), zu]
            radius = 0.5d0*(xu - xl)
            height = zu - zl

            origin_bottom(1:2) = origin(1:2)
            origin_bottom(3) = origin(3) - height

            ! Surface
            allocate (pplane)
            pplane = new_planeXYZWithCircleHoleZ(origin, radius)
            pbound => pplane
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pplane)
            end if

            ! Wall surface
            allocate (pcylinder)
            pcylinder = new_cylinderZ(origin_bottom, radius, height)
            pbound => pcylinder
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcylinder)
            end if

            ! Bottom surface
            origin_bottom = [origin(1), origin(2), zl]
            allocate (pcircle)
            pcircle = new_CircleZ(origin_bottom, radius)
            pbound => pcircle
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if
        end subroutine

        subroutine add_hyperboloid_hole_surface
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZWithCircleHole), pointer :: pplane
            type(t_HyperboloidXYZ), pointer :: phyperboloid
            type(t_CircleXYZ), pointer :: pcircle

            double precision :: origin(3), origin_hyperboloid(3), origin_bottom(3)
            double precision :: max_radius, min_radius, height

            max_radius = 0.5d0*(xu - xl)
            min_radius = rcurv*max_radius
            height = zu - zl

            origin(:) = [0.5d0*(xl + xu), 0.5d0*(yl + yu), zu]
            origin_hyperboloid(:) = [origin(1), origin(2), zu - 0.5d0*height]
            origin_bottom(:) = [origin(1), origin(2), zl]

            ! Surface
            allocate (pplane)
            pplane = new_planeXYZWithCircleHoleZ(origin, max_radius)
            pbound => pplane
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pplane)
            end if

            ! Wall surface
            allocate (phyperboloid)
            phyperboloid = new_HyperboloidZ(origin_hyperboloid, max_radius, min_radius, height)
            pbound => phyperboloid
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (phyperboloid)
            end if

            ! Bottom surface
            allocate (pcircle)
            pcircle = new_CircleZ(origin_bottom, max_radius)
            pbound => pcircle
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if
        end subroutine

        subroutine add_ellipsoid_hole_surface
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZWithCircleHole), pointer :: pplane
            type(t_EllipsoidXYZ), pointer :: pellipsoid
            type(t_CircleXYZ), pointer :: pcircle

            double precision :: origin(3), origin_ellipsoid(3), origin_bottom(3)
            double precision :: max_radius, min_radius, height

            min_radius = 0.5d0*(xu - xl)
            max_radius = rcurv*min_radius
            height = zu - zl

            origin(:) = [0.5d0*(xl + xu), 0.5d0*(yl + yu), zu]
            origin_ellipsoid(:) = [origin(1), origin(2), zu - 0.5d0*height]
            origin_bottom(:) = [origin(1), origin(2), zl]

            ! Surface
            allocate (pplane)
            pplane = new_planeXYZWithCircleHoleZ(origin, min_radius)
            pbound => pplane
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pplane)
            end if

            ! Wall surface
            allocate (pellipsoid)
            pellipsoid = new_EllipsoidZ(origin_ellipsoid, min_radius, max_radius, height)
            pbound => pellipsoid
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pellipsoid)
            end if

            ! Bottom surface
            allocate (pcircle)
            pcircle = new_CircleZ(origin_bottom, min_radius)
            pbound => pcircle
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pcircle)
            end if
        end subroutine

        subroutine add_cover_all
            class(t_Boundary), pointer :: pbound
            type(t_RectangleXYZ), pointer :: prect
            double precision :: origin(3)
            double precision :: rnx, rny
            integer :: iboundary, nboundaries_init
            double precision :: xmin, ymin, xmax, ymax

            nboundaries_init = boundaries%nboundaries
            rnx = nx
            rny = ny

            ! 0. Bottom boundary
            origin = [0d0, 0d0, 0d0]
            allocate (prect)
            prect = new_rectangleZ(origin, rnx, rny)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 1. -X boundary
            origin = [0d0, 0d0, 0d0]
            allocate (prect)
            prect = new_rectangleX(origin, rny, zu)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 2. +X boundary
            origin = [rnx, 0d0, 0d0]
            allocate (prect)
            prect = new_rectangleX(origin, rny, zu)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 3. -Y boundary
            origin = [0d0, 0d0, 0d0]
            allocate (prect)
            prect = new_rectangleY(origin, zu, rnx)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! 4. +Y boundary
            origin = [0d0, rny, 0d0]
            allocate (prect)
            prect = new_rectangleY(origin, zu, rnx)
            pbound => prect
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (prect)
            end if

            ! Tag to indicate no photoelectron emission.
            do iboundary = nboundaries_init + 1, boundaries%nboundaries
                boundaries%boundaries(iboundary)%ref%material%tag = -1
            end do
        end subroutine

        subroutine add_disk(axis)
            integer, intent(in) :: axis

            class(t_Boundary), pointer :: pbound
            type(t_CylinderXYZ), pointer :: pouter_cylinder
            type(t_CylinderXYZ), pointer :: pinner_cylinder
            type(t_DonutXYZ), pointer :: plower_donut
            type(t_DonutXYZ), pointer :: pupper_donut

            double precision :: height
            double precision :: lower_origin(3)
            double precision :: upper_origin(3)
            double precision :: radius
            double precision :: inner_radius

            height = disk_height(itype)
            lower_origin(1:3) = disk_origin(itype, 1:3)
            upper_origin(1:3) = disk_origin(itype, 1:3)
            upper_origin(axis) = upper_origin(axis) + height
            radius = disk_radius(itype)
            inner_radius = disk_inner_radius(itype)

            ! Outer cylinder
            allocate (pouter_cylinder)
            pouter_cylinder = new_cylinderXYZ(axis, lower_origin, radius, height)
            pbound => pouter_cylinder
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pouter_cylinder)
            end if

            ! Inner cylinder
            allocate (pinner_cylinder)
            pinner_cylinder = new_cylinderXYZ(axis, lower_origin, inner_radius, height)
            pbound => pinner_cylinder
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pinner_cylinder)
            end if

            ! Lower donut
            allocate (plower_donut)
            plower_donut = new_DonutXYZ(axis, lower_origin, inner_radius, height)
            pbound => plower_donut
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (plower_donut)
            end if

            ! Upper donut
            allocate (pupper_donut)
            pupper_donut = new_DonutXYZ(axis, upper_origin, inner_radius, height)
            pbound => pupper_donut
            if (pbound%is_overlap(sdoms, extent=extent)) then
                call boundaries%add_boundary(pbound)
            else
                deallocate (pupper_donut)
            end if
        end subroutine

    end function

end module

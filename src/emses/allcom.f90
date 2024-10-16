module allcom
    !! This module defines common variables and parameters used across the application.

    use m_vector

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)
        !! The mathematical constant pi
    double precision, parameter :: DEG2RAD = pi/180d0
        !! Coefficient to convert from degrees to radians
    double precision, parameter :: RAD2DEG = 180d0/pi
        !! Coefficient to convert from radians to degrees

    integer, parameter :: max_nspec = 10
        !! Maximum number of species
    integer, parameter :: nboundary_types = 10
        !! Number of boundary types
    integer, parameter :: max_nepl = 15
        !! Maximum number of emission surfaces

    ! /esorem/
    integer :: emflag = 0
        !! Electromagnetic flag

    ! /plasma/
    double precision :: wp(max_nspec)
        !! Plasma frequency for each species
    double precision :: wc
        !! Cyclotron frequency
    double precision :: phixy
        !! Background magnetic field angle in XY plane
    double precision :: phiz
        !! Background magnetic field angle from z-axis

    ! /tmgrid/
    double precision :: dt
        !! Time step width
    integer :: nx
        !! Number of grid cells in x-axis
    integer :: ny
        !! Number of grid cells in y-axis
    integer :: nz
        !! Number of grid cells in z-axis

    ! /system/
    integer :: nspec
        !! Number of species
    integer :: npbnd(3, max_nspec)
        !! The boundary conditoins for each species and direction

    ! /intp/
    double precision :: qm(max_nspec)
        !! Charge to mass ratio for each species
    double precision :: path(max_nspec)
        !! Thermal velocity along the background magnetic field for each species
    double precision :: peth(max_nspec)
        !! Thermal velocity perpendicular to the background magnetic field for each species
    double precision :: vdri(max_nspec)
        !! Drift velocity for each species
    double precision :: vdthz(max_nspec)
        !! Drift velocity angle from z-axis
    double precision :: vdthxy(max_nspec)
        !! Drift velocity angle in XY plane
    double precision :: spa(max_nspec) = 0d0
        !! Drift velocity along the background magnetic field for each species
    double precision :: spe(max_nspec) = 0d0
        !! Drift velocity perpendicular to the background magnetic field for each species
    double precision :: speth(max_nspec) = 0d0
        !! The angle of drift velocity perpendicular to the background magnetic field for each species

    ! /ptcond/
    real(kind=8) :: zssurf = -9999.0d0
        !! Surface Z-coordinate
    real(kind=8) :: xlrechole(2), ylrechole(2), zlrechole(2)
        !! Lower rectangle hole coordinates
    real(kind=8) :: xurechole(2), yurechole(2), zurechole(2)
        !! Upper rectangle hole coordinates

    ! Boundary type
    ! "flat-surface"
    ! "rectangle-hole"
    ! "cylinder-hole"
    ! "hyperboloid-hole"
    ! "rectangle"
    ! "circle"
    ! "cuboid"
    ! "disk"
    character(len=30) :: boundary_type = "flat-surface"
        !! Type of boundary condition
    character(len=30) :: boundary_types(nboundary_types) = "none"
        !! List of boundary types (activated when boundary_type = "complex")
    double precision :: cylinder_origin(nboundary_types, 3) = 0.0d0
        !! Origins of cylinder boundaries
    double precision :: cylinder_radius(nboundary_types) = 0.0d0
        !! Radii of cylinder boundaries
    double precision :: cylinder_height(nboundary_types) = 0.0d0
        !! Heights of cylinder boundaries
    double precision :: rcurv = 1.0d0
        !! Minimum radius relative to the entrance radius of the hyperbolic hole (r_min == rcurv*r_max). <br>
        !! If rcurv == 1.0d0, then r_min == r_max. <br>
        !! Else if rcurv == 0.5d0 then r_min == 0.5d0 * r_max
    double precision :: rectangle_shape(nboundary_types, 6) = 0.0d0
        !! Shapes of rectangle boundaries
    double precision :: sphere_origin(nboundary_types, 3) = 0.0d0
        !! Origins of sphere boundaries
    double precision :: sphere_radius(nboundary_types) = 0.0d0
        !! Radius of sphere boundaries
    double precision :: circle_origin(nboundary_types, 3) = 0.0d0
        !! Origins of circle boundaries
    double precision :: circle_radius(nboundary_types) = 0.0d0
        !! Radius of circle boundaries
    double precision :: cuboid_shape(nboundary_types, 6) = 0.0d0
        !! Shapes of cuboid boundaries
    double precision :: disk_origin(nboundary_types, 3) = 0.0d0
        !! Origins of disk boundaries
    double precision :: disk_height(nboundary_types) = 0.0d0
        !! Heights of disk boundaries
    double precision :: disk_radius(nboundary_types) = 0.0d0
        !! Radius of disk boundaries
    double precision :: disk_inner_radius(nboundary_types) = 0.0d0
        !! Inner radius of disk boundaries
    double precision :: plane_with_circle_hole_zlower(nboundary_types) = 0.0d0
    double precision :: plane_with_circle_hole_height(nboundary_types) = 0.0d0
    double precision :: plane_with_circle_hole_radius(nboundary_types) = 0.0d0

    ! /emissn/
    integer :: nflag_emit(max_nspec)
        !! Emission flags for each species
    integer :: nepl(max_nspec) = 0
        !! Number of emission surfaces for each species
    double precision :: curf(max_nspec)
        !! Current from emission surfaces
    integer :: nemd(max_nepl)
        !! Emission directions
    double precision :: curfs(max_nepl)
        !! Current from each emission surface
    double precision :: xmine(max_nepl), xmaxe(max_nepl)
        !! Minimum and maximum X-coordinates for each emission surface
    double precision :: ymine(max_nepl), ymaxe(max_nepl)
        !! Minimum and maximum Y-coordinates for each emission surface
    double precision :: zmine(max_nepl), zmaxe(max_nepl)
        !! Minimum and maximum Z-coordinates for each emission surface
    double precision :: thetaz(max_nepl), thetaxy(max_nepl)
        !! Emission angles in the Z direction and XY plane for each emission surface

contains

    function vdri_vector(ispec) result(ret)
        !! Calculate the drift velocity vector.

        integer, intent(in) :: ispec
            !! Species index
        double precision :: ret(3)
            !! Resulting drift velocity vector

        double precision :: vdri_from_spa(3)
        double precision :: vdri_from_vdri(3)

        vdri_from_spa(:) = [spe(ispec), 0d0, spa(ispec)]
        vdri_from_spa(:) = rot3d_z(vdri_from_spa, speth(ispec))
        vdri_from_spa(:) = rot3d_y(vdri_from_spa, phiz*DEG2RAD)
        vdri_from_spa(:) = rot3d_z(vdri_from_spa, phixy*DEG2RAD)

        vdri_from_vdri(:) = [0d0, 0d0, vdri(ispec)]
        vdri_from_vdri(:) = rot3d_y(vdri_from_vdri, vdthz(ispec)*DEG2RAD)
        vdri_from_vdri(:) = rot3d_z(vdri_from_vdri, vdthxy(ispec)*DEG2RAD)

        ret(:) = vdri_from_spa + vdri_from_vdri
    end function

    function vth_vector(ispec) result(ret)
        !! Calculate the thermal velocity vector.

        integer, intent(in) :: ispec
            !! Species index
        double precision :: ret(3)
            !! Resulting thermal velocity vector

        ret(:) = [peth(ispec), peth(ispec), path(ispec)]
        ret(:) = rot3d_y(ret, phiz*DEG2RAD)
        ret(:) = rot3d_z(ret, phixy*DEG2RAD)
    end function

    function emission_vdri_vector(ispec, iepl) result(ret)
        !! Calculate the drift velocity vector.

        integer, intent(in) :: ispec
            !! Species index
        integer, intent(in) :: iepl
            !! Emission surface index
        double precision :: ret(3)
            !! Resulting drift velocity vector

        double precision :: vdri_from_spa(3)
        double precision :: vdri_from_vdri(3)

        vdri_from_spa(:) = [spe(ispec), 0d0, spa(ispec)]
        vdri_from_spa(:) = rot3d_z(vdri_from_spa, speth(ispec))
        vdri_from_spa(:) = rot3d_y(vdri_from_spa, phiz*DEG2RAD)
        vdri_from_spa(:) = rot3d_z(vdri_from_spa, phixy*DEG2RAD)

        block
            double precision :: thz, thxy

            call nemd2angle(iepl, thz, thxy)

            vdri_from_vdri(:) = [0d0, 0d0, vdri(ispec)]
            vdri_from_vdri(:) = rot3d_y(vdri_from_vdri, thz*DEG2RAD)
            vdri_from_vdri(:) = rot3d_z(vdri_from_vdri, thxy*DEG2RAD)
        end block

        ret(:) = vdri_from_spa + vdri_from_vdri
    end function

    function emission_vth_vector(ispec, iepl) result(ret)
        !! Calculate the thermal velocity vector.

        integer, intent(in) :: ispec
            !! Species index
        integer, intent(in) :: iepl
            !! Emission surface index
        double precision :: ret(3)
            !! Resulting thermal velocity vector

        ret(:) = [peth(ispec), peth(ispec), path(ispec)]
        block
            double precision :: thz, thxy

            call nemd2angle(iepl, thz, thxy)

            ret(:) = [peth(ispec), peth(ispec), path(ispec)]
            ret(:) = rot3d_y(ret, thz*DEG2RAD)
            ret(:) = rot3d_z(ret, thxy*DEG2RAD)
        end block
    end function

    subroutine nemd2angle(iepl, thz, thxy)
        !! Convert emission direction to angle
        integer, intent(in) :: iepl
            !! Emission surface index
        double precision, intent(out) :: thz
            !! Angle from z-axis
        double precision, intent(out) :: thxy
            !! Angle in XY plane

        if (abs(nemd(iepl)) == 1) then ! X-direction
            thz = 90
            thxy = 0
        else if (abs(nemd(iepl)) == 2) then ! Y-direction
            thz = 90
            thxy = 90
        else if (abs(nemd(iepl)) == 3) then ! Z-direction
            thz = 0
            thxy = 0
        end if

        if (nemd(iepl) < 0) then ! Minus direction
            thz = thz - 180
        end if
    end subroutine

end module

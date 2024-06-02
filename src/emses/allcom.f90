module allcom
    !! This module defines common variables and parameters used across the application.

    implicit none

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

    ! /emissn/
    integer :: nflag_emit(max_nspec)
        !! Emission flags for each species
    integer :: nepl(max_nspec)
        !! Number of emission surfaces for each species
    double precision :: curf
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

end module

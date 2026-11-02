module m_allcom
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

    ! Every namelist variable below has an explicit default. `read_namelist`
    ! restores all of them through `reset_input_parameters` before reading a
    ! file, so values omitted from one input never leak in from the previous
    ! one (issue #40). Keep the two lists in sync.

    ! /plasma/
    double precision :: wp(max_nspec) = 0d0
        !! Plasma frequency for each species
    double precision :: wc = 0d0
        !! Cyclotron frequency
    double precision :: phixy = 0d0
        !! Background magnetic field angle in XY plane
    double precision :: phiz = 0d0
        !! Background magnetic field angle from z-axis

    ! /tmgrid/
    double precision :: dt = 0d0
        !! Time step width
    integer :: nx = 0
        !! Number of grid cells in x-axis
    integer :: ny = 0
        !! Number of grid cells in y-axis
    integer :: nz = 0
        !! Number of grid cells in z-axis

    ! /system/
    integer :: nspec = 0
        !! Number of species
    integer :: npbnd(3, max_nspec) = 0
        !! The boundary conditoins for each species and direction

    ! /intp/
    double precision :: qm(max_nspec) = 0d0
        !! Charge to mass ratio for each species
    double precision :: path(max_nspec) = 0d0
        !! Thermal velocity along the background magnetic field for each species
    double precision :: peth(max_nspec) = 0d0
        !! Thermal velocity perpendicular to the background magnetic field for each species
    double precision :: vdri(max_nspec) = 0d0
        !! Drift velocity for each species
    double precision :: vdthz(max_nspec) = 0d0
        !! Drift velocity angle from z-axis
    double precision :: vdthxy(max_nspec) = 0d0
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
    real(kind=8) :: xlrechole(2) = 0d0, ylrechole(2) = 0d0, zlrechole(2) = 0d0
        !! Lower rectangle hole coordinates
    real(kind=8) :: xurechole(2) = 0d0, yurechole(2) = 0d0, zurechole(2) = 0d0
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
    integer :: boundary_conductor_id(nboundary_types) = -1
        !! Conductor ID assigned to each complex boundary
    double precision :: cylinder_origin(3, nboundary_types) = 0.0d0
        !! Origins of cylinder boundaries
    double precision :: cylinder_radius(nboundary_types) = 0.0d0
        !! Radii of cylinder boundaries
    double precision :: cylinder_height(nboundary_types) = 0.0d0
        !! Heights of cylinder boundaries
    double precision :: rcurv = 1.0d0
        !! Minimum radius relative to the entrance radius of the hyperbolic hole (r_min == rcurv*r_max). <br>
        !! If rcurv == 1.0d0, then r_min == r_max. <br>
        !! Else if rcurv == 0.5d0 then r_min == 0.5d0 * r_max
    double precision :: rectangle_shape(6, nboundary_types) = 0.0d0
        !! Shapes of rectangle boundaries
    double precision :: sphere_origin(3, nboundary_types) = 0.0d0
        !! Origins of sphere boundaries
    double precision :: sphere_radius(nboundary_types) = 0.0d0
        !! Radius of sphere boundaries
    double precision :: circle_origin(3, nboundary_types) = 0.0d0
        !! Origins of circle boundaries
    double precision :: circle_radius(nboundary_types) = 0.0d0
        !! Radius of circle boundaries
    double precision :: cuboid_shape(6, nboundary_types) = 0.0d0
        !! Shapes of cuboid boundaries
    double precision :: disk_origin(3, nboundary_types) = 0.0d0
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
    double precision :: plane_with_circle_origin(3, nboundary_types) = 0.0d0
        !! Origins of planes with circular apertures
    double precision :: plane_with_circle_radius(nboundary_types) = 0.0d0
        !! Radii of circular apertures defined on planes
    double precision :: conductivity = 0.0d0
        !! Boundary conductivity
    integer :: max_bounce_count = 3
        !! Maximum number of repeated collisions treated on a boundary
    double precision :: boundary_mirror_reflection_rate(nboundary_types, max_nspec) = 0.0d0
        !! Mirror reflection rate for each boundary and species
    double precision :: boundary_reversal_reflection_rate(nboundary_types, max_nspec) = 0.0d0
        !! Reversal reflection rate for each boundary and species
    double precision :: boundary_mirror_reflect_alpha(nboundary_types, max_nspec) = 0.0d0
        !! Mirror reflection angular parameter for each boundary and species
    double precision :: boundary_reversal_reflect_alpha(nboundary_types, max_nspec) = 0.0d0
        !! Reversal reflection angular parameter for each boundary and species
    double precision :: boundary_mirror_reflect_energy_loss_frac(nboundary_types, max_nspec) = 0.0d0
        !! Energy loss fraction for mirror reflection
    double precision :: boundary_reversal_reflect_energy_loss_frac(nboundary_types, max_nspec) = 0.0d0
        !! Energy loss fraction for reversal reflection
    logical :: enable_secondary_electron_emission(nboundary_types, max_nspec) = .false.
        !! Flag to enable secondary electron emission on each boundary and species
    integer :: boundary_se_model_type(nboundary_types, max_nspec) = 0
        !! Secondary electron emission model for each boundary and species
    double precision :: boundary_se_const_yield(nboundary_types, max_nspec) = 0.0d0
        !! Constant SEE yield for each boundary and species
    double precision :: boundary_se_yield_max(nboundary_types, max_nspec) = 0.0d0
        !! Maximum SEE yield for each boundary and species
    double precision :: boundary_se_energy_max(nboundary_types, max_nspec) = 300.0d0
        !! Energy at which the SEE yield peaks
    integer :: boundary_se_species_id(nboundary_types, max_nspec) = 1
        !! Emitted species ID for secondary electron emission

    ! /emissn/
    integer :: nflag_emit(max_nspec) = 0
        !! Emission flags for each species
    integer :: nepl(max_nspec) = 0
        !! Number of emission surfaces for each species
    double precision :: curf(max_nspec) = 0d0
        !! Current from emission surfaces
    integer :: nemd(max_nepl) = 0
        !! Emission directions
    double precision :: curfs(max_nepl) = -9999d0
        !! Current from each emission surface
    double precision :: xmine(max_nepl) = 0d0, xmaxe(max_nepl) = 0d0
        !! Minimum and maximum X-coordinates for each emission surface
    double precision :: ymine(max_nepl) = 0d0, ymaxe(max_nepl) = 0d0
        !! Minimum and maximum Y-coordinates for each emission surface
    double precision :: zmine(max_nepl) = 0d0, zmaxe(max_nepl) = 0d0
        !! Minimum and maximum Z-coordinates for each emission surface
    double precision :: thetaz(max_nepl) = 0d0, thetaxy(max_nepl) = 0d0
        !! Emission angles in the Z direction and XY plane for each emission surface

    logical :: use_raycast = .false.
        !! Enable raycast-based photoelectron probability for species with nflag_emit == 2
    double precision :: ray_zenith_angle_deg(max_nspec) = 9999d0
        !! Sentinel value 9999d0 means "fall back to vdthz(ispec)". Otherwise the
        !! ray direction follows the same rotation convention as `vdri_vector`
        !! (zenith from +z with sign -ray_zenith_angle_deg around y).
    double precision :: ray_azimuth_angle_deg(max_nspec) = 9999d0
        !! Sentinel value 9999d0 means "fall back to vdthxy(ispec)".

contains

    subroutine reset_input_parameters()
        !! Restore every namelist variable to its declared default.
        !!
        !! Namelist reads only overwrite the keys present in the file, so a
        !! value set by one input would otherwise persist into the next input
        !! read in the same process. Called by `read_namelist` before reading.

        ! /esorem/
        emflag = 0

        ! /plasma/
        wp = 0d0
        wc = 0d0
        phixy = 0d0
        phiz = 0d0

        ! /tmgrid/
        dt = 0d0
        nx = 0
        ny = 0
        nz = 0

        ! /system/
        nspec = 0
        npbnd = 0

        ! /intp/
        qm = 0d0
        path = 0d0
        peth = 0d0
        vdri = 0d0
        vdthz = 0d0
        vdthxy = 0d0
        spa = 0d0
        spe = 0d0
        speth = 0d0

        ! /ptcond/
        zssurf = -9999.0d0
        xlrechole = 0d0
        ylrechole = 0d0
        zlrechole = 0d0
        xurechole = 0d0
        yurechole = 0d0
        zurechole = 0d0
        boundary_type = "flat-surface"
        boundary_types = "none"
        boundary_conductor_id = -1
        cylinder_origin = 0.0d0
        cylinder_radius = 0.0d0
        cylinder_height = 0.0d0
        rcurv = 1.0d0
        rectangle_shape = 0.0d0
        sphere_origin = 0.0d0
        sphere_radius = 0.0d0
        circle_origin = 0.0d0
        circle_radius = 0.0d0
        cuboid_shape = 0.0d0
        disk_origin = 0.0d0
        disk_height = 0.0d0
        disk_radius = 0.0d0
        disk_inner_radius = 0.0d0
        plane_with_circle_hole_zlower = 0.0d0
        plane_with_circle_hole_height = 0.0d0
        plane_with_circle_hole_radius = 0.0d0
        plane_with_circle_origin = 0.0d0
        plane_with_circle_radius = 0.0d0
        conductivity = 0.0d0
        max_bounce_count = 3
        boundary_mirror_reflection_rate = 0.0d0
        boundary_reversal_reflection_rate = 0.0d0
        boundary_mirror_reflect_alpha = 0.0d0
        boundary_reversal_reflect_alpha = 0.0d0
        boundary_mirror_reflect_energy_loss_frac = 0.0d0
        boundary_reversal_reflect_energy_loss_frac = 0.0d0
        enable_secondary_electron_emission = .false.
        boundary_se_model_type = 0
        boundary_se_const_yield = 0.0d0
        boundary_se_yield_max = 0.0d0
        boundary_se_energy_max = 300.0d0
        boundary_se_species_id = 1

        ! /emissn/
        nflag_emit = 0
        nepl = 0
        curf = 0d0
        nemd = 0
        curfs = -9999d0
        xmine = 0d0
        xmaxe = 0d0
        ymine = 0d0
        ymaxe = 0d0
        zmine = 0d0
        zmaxe = 0d0
        thetaz = 0d0
        thetaxy = 0d0
        use_raycast = .false.
        ray_zenith_angle_deg = 9999d0
        ray_azimuth_angle_deg = 9999d0
    end subroutine

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
        vdri_from_vdri(:) = rot3d_y(vdri_from_vdri, -vdthz(ispec)*DEG2RAD)
        vdri_from_vdri(:) = rot3d_z(vdri_from_vdri, vdthxy(ispec)*DEG2RAD)

        ret(:) = vdri_from_spa + vdri_from_vdri
    end function

    function vth_vector(ispec) result(ret)
        !! Cartesian marginal standard deviations for envelope support.

        integer, intent(in) :: ispec
            !! Species index
        double precision :: ret(3)
        double precision :: basis(3, 3), scales(3)
        integer :: i

        scales = [peth(ispec), peth(ispec), path(ispec)]
        basis = thermal_basis()
        do i = 1, 3
            ret(i) = sqrt(sum((basis(i, :)*scales)**2))
        end do
    end function

    function thermal_basis() result(ret)
        !! Principal thermal axes, using the same rotation as the drift.
        double precision :: ret(3, 3), axis(3)
        integer :: i

        do i = 1, 3
            axis = 0d0
            axis(i) = 1d0
            axis = rot3d_y(axis, phiz*DEG2RAD)
            ret(:, i) = rot3d_z(axis, phixy*DEG2RAD)
        end do
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

        integer :: axis0, axis1, axis2

        vdri_from_spa(:) = [spe(ispec), 0d0, spa(ispec)]
        vdri_from_spa(:) = rot3d_z(vdri_from_spa, speth(ispec))
        vdri_from_spa(:) = rot3d_y(vdri_from_spa, phiz*DEG2RAD)
        vdri_from_spa(:) = rot3d_z(vdri_from_spa, phixy*DEG2RAD)

        axis0 = abs(nemd(iepl))
        axis1 = mod(axis0, 3) + 1
        axis2 = mod(axis0 + 1, 3) + 1

        if (nemd(iepl) > 0) then
            ret(axis0) = vdri_from_spa(3)
        else
            ret(axis0) = -vdri_from_spa(3)
        end if

        ret(axis1) = vdri_from_spa(1)
        ret(axis2) = vdri_from_spa(2)

    end function

    function emission_vth_vector(ispec, iepl) result(ret)
        !! Emission-frame Gaussian standard deviations (not signed vectors).

        integer, intent(in) :: ispec
            !! Species index
        integer, intent(in) :: iepl
            !! Emission surface index
        double precision :: ret(3)
            !! Resulting thermal velocity vector

        ! MPIEMSES3D samples the normal flux distribution with path and both
        ! tangential Gaussians with peth, then permutes axes and signs the normal.
        ! Standard deviations are positive and do not acquire that sign.
        ret(:) = peth(ispec)
        ret(abs(nemd(iepl))) = path(ispec)
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

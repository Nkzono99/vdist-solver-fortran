!
module m_namelist
    use, intrinsic :: iso_c_binding

    use allcom

    implicit none

    private
    public read_namelist

    namelist /esorem/ emflag

    namelist /plasma/ wp, wc, phixy, phiz

    namelist /tmgrid/ dt, nx, ny, nz

    namelist /system/ nspec, npbnd

    namelist /intp/ qm, path, peth, vdri, vdthz, vdthxy, &
        spa, spe, speth

    namelist /ptcond/ boundary_type, boundary_types, &
        zssurf, &
        xlrechole, ylrechole, zlrechole, xurechole, yurechole, zurechole, &
        boundary_conductor_id, &
        cylinder_origin, cylinder_radius, cylinder_height, rcurv, &
        rectangle_shape, &
        sphere_origin, sphere_radius, &
        circle_origin, circle_radius, &
        cuboid_shape, &
        disk_origin, disk_height, disk_radius, disk_inner_radius, &
        conductivity, &
        plane_with_circle_hole_zlower, &
        plane_with_circle_hole_height, &
        plane_with_circle_hole_radius, &
        plane_with_circle_origin, &
        plane_with_circle_radius, &
        max_bounce_count, &
        boundary_mirror_reflection_rate, &
        boundary_mirror_reflect_alpha, &
        boundary_reversal_reflection_rate, &
        boundary_reversal_reflect_alpha, &
        boundary_mirror_reflect_energy_loss_frac, &
        boundary_reversal_reflect_energy_loss_frac, &
        enable_secondary_electron_emission, &
        boundary_se_model_type, &
        boundary_se_const_yield, &
        boundary_se_yield_max, &
        boundary_se_energy_max, &
        boundary_se_species_id

    namelist /emissn/ nflag_emit, &
        nepl, curf, nemd, curfs, &
        xmaxe, xmine, ymaxe, ymine, zmaxe, zmine, &
        thetaz, thetaxy

contains

    subroutine read_namelist(inppath)
        !! Read namelists from a specified input file path.

        character(len=*), intent(in) :: inppath
            !! Path to the input file

        integer :: iostat
            !! I/O status variable

        open (10, file=trim(inppath), status='old', iostat=iostat)

        rewind (10); read (10, nml=esorem, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=esorem not found"
        end if

        rewind (10); read (10, nml=plasma, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=plasma not found"
        end if

        rewind (10); read (10, nml=tmgrid, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=tmgrid not found"
        end if

        rewind (10); read (10, nml=system, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=system not found"
        end if

        rewind (10); read (10, nml=intp, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=intp not found"
        end if

        rewind (10); read (10, nml=ptcond, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=ptcond not found"
        end if

        rewind (10); read (10, nml=emissn, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=emissn not found"
        end if

        close (10)
    end subroutine

end module

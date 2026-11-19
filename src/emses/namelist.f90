!
module m_namelist
    use, intrinsic :: iso_c_binding

    use m_allcom

    implicit none

    private
    public read_namelist

    namelist /esorem/ emflag

    namelist /plasma/ wp, wc, phixy, phiz

    namelist /tmgrid/ dt, nx, ny, nz

    namelist /system/ nspec, npbnd

    namelist /intp/ qm, path, peth, vdri, vdthz, vdthxy, &
        spa, spe, speth, vdx, vdy, vdz

    namelist /ptcond/ boundary_type, boundary_types, &
        zssurf, &
        xlrechole, ylrechole, zlrechole, xurechole, yurechole, zurechole, &
        boundary_conductor_id, &
        cylinder_origin, cylinder_radius, cylinder_height, rcurv, &
        rectangle_shape, rectangle_corners, &
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
        nepl, curf, nemd, curfs, abvdem, &
        xmaxe, xmine, ymaxe, ymine, zmaxe, zmine, &
        thetaz, thetaxy, decrease_current_according_to_incident_angle, &
        use_raycast, &
        ray_zenith_angle_deg, ray_azimuth_angle_deg

contains

    subroutine read_namelist(inppath, status)
        !! Read namelists from a specified input file path.
        !!
        !! Without `status`, an unreadable file or a malformed group is skipped
        !! as before. With `status`, the first such error is returned (nonzero)
        !! and reading stops; a missing group is still only a warning.

        character(len=*), intent(in) :: inppath
            !! Path to the input file
        integer, intent(out), optional :: status
            !! 0 on success, otherwise the failing iostat

        integer :: unit, iostat

        ! Namelist reads leave omitted keys untouched; start from the defaults
        ! so one input cannot inherit values from a previous one (issue #40).
        call reset_input_parameters()
        if (present(status)) status = 0

        open (newunit=unit, file=trim(inppath), status='old', iostat=iostat)
        if (iostat /= 0) then
            if (present(status)) status = iostat
            return
        end if

        rewind (unit); read (unit, nml=esorem, iostat=iostat)
        if (stop_reading('esorem')) return
        rewind (unit); read (unit, nml=plasma, iostat=iostat)
        if (stop_reading('plasma')) return
        rewind (unit); read (unit, nml=tmgrid, iostat=iostat)
        if (stop_reading('tmgrid')) return
        rewind (unit); read (unit, nml=system, iostat=iostat)
        if (stop_reading('system')) return
        rewind (unit); read (unit, nml=intp, iostat=iostat)
        if (stop_reading('intp')) return
        rewind (unit); read (unit, nml=ptcond, iostat=iostat)
        if (stop_reading('ptcond')) return
        rewind (unit); read (unit, nml=emissn, iostat=iostat)
        if (stop_reading('emissn')) return

        close (unit)

    contains

        logical function stop_reading(group)
            character(len=*), intent(in) :: group

            stop_reading = .false.
            if (iostat == -1) then
                print *, "Warning.Input: nml=", group, " not found"
            else if (iostat > 0 .and. present(status)) then
                status = iostat
                close (unit)
                stop_reading = .true.
            end if
        end function
    end subroutine

end module

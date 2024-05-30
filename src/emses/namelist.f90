!
module m_namelist
    use, intrinsic :: iso_c_binding

    use allcom

    implicit none

    namelist /esorem/ emflag

    namelist /plasma/ wp, wc, phixy, phiz

    namelist /tmgrid/ dt, nx, ny, nz

    namelist /system/ nspec, npbnd

    namelist /intp/ qm, path, peth, vdri, vdthz, vdthxy

    namelist /ptcond/ boundary_type, boundary_types, &
        zssurf, &
        xlrechole, ylrechole, zlrechole, xurechole, yurechole, zurechole, &
        cylinder_origin, cylinder_radius, cylinder_height, rcurv, &
        rectangle_shape, &
        sphere_origin, sphere_radius, &
        circle_origin, circle_radius, &
        cuboid_shape, &
        disk_origin, disk_height, disk_radius, disk_inner_radius

    namelist /emissn/ nflag_emit, &
        nepl, curf, nemd, curfs, &
        xmaxe, xmine, ymaxe, ymine, zmaxe, zmine, &
        thetaz, thetaxy

contains

    subroutine read_namelist(inppath)
        character(len=*), intent(In) :: inppath
        integer :: iostat

        ! Convert the array of characters to string.
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

!
module m_namelist
    use iso_c_binding
    !
    !   ____________________________________________________________
    !
    !                    M O D U L E   N A M E L S
    !   ____________________________________________________________
    !
    !   ............................................................
    !   .                                                          .
    !   .        module for common declaration of namelist         .
    !   ............................................................
    !
    !-------------------- common declaration of namelist
    use allcom

    namelist /esorem/ emflag

    namelist /plasma/ wp, wc, phixy, phiz

    namelist /tmgrid/ dt, nx, ny, nz

    namelist /system/ nspec, npbnd

    namelist /intp/ qm, path, peth, vdri, vdthz, vdthxy

    namelist /ptcond/ zssurf, &
        xlrechole, ylrechole, zlrechole, xurechole, yurechole, zurechole, &
        boundary_type, boundary_types, &
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
        character(1, c_char), intent(In) :: inppath(:)
        integer :: iostat
        integer :: i

        block ! Convert the array of characters to string.
            character(size(inppath)) :: s
            do i = 1, size(inppath)
                s(i:i) = inppath(i)
            end do
            open (1, file=s)
        end block

        rewind (1); read (1, nml=esorem, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=esorem not found"
        end if

        rewind (1); read (1, nml=plasma, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=plasma not found"
        end if

        rewind (1); read (1, nml=tmgrid, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=tmgrid not found"
        end if

        rewind (1); read (1, nml=system, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=system not found"
        end if

        rewind (1); read (1, nml=intp, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=intp not found"
        end if

        rewind (1); read (1, nml=ptcond, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=ptcond not found"
        end if

        rewind (1); read (1, nml=emissn, IOSTAT=iostat)
        if (iostat .eq. -1) then
            print *, "Warning.Input: nml=emissn not found"
        end if

        close (1)
    end subroutine

end module

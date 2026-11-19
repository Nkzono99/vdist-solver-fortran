module m_emses_input_validation
    !! Recoverable checks of a vdist EMSES input before building a simulator.
    !!
    !! The geometry builders copied from MPIEMSES3D stop the process on invalid
    !! shapes, and settings that EMSES uses but vdist cannot replay would be
    !! ignored silently. `check_emses_input` reports both as a message. The
    !! C entry `validate_emses_input` lets Python raise ValueError before
    !! calling the void solver routines. Settings visible only in the original
    !! EMSES namelist are checked in `vdsolverf/emses/source_support.py`.

    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use m_namelist, only: read_namelist
    use m_allcom
    use m_collision_boundary_types, only: is_collision_boundary_entry_empty, &
                                          is_collision_boundary_type_supported, &
                                          is_collision_object_boundary_type

    implicit none

    private
    public :: check_emses_input
    public :: validate_emses_input

contains

    integer(c_int) function validate_emses_input(filename, length, ispec, regular_pdf, message, capacity) &
        result(status) bind(c)
        !! Read `filename` and check it for species `ispec` (one-based).
        !!
        !! Returns 0 when valid; otherwise 1 and a NUL-terminated reason in
        !! `message(1:capacity)`. `regular_pdf /= 0` additionally requires
        !! positive thermal widths, since cold beams have no regular density.
        character(kind=c_char), intent(in) :: filename(*)
        integer(c_int), value, intent(in) :: length
        integer(c_int), value, intent(in) :: ispec
        integer(c_int), value, intent(in) :: regular_pdf
        character(kind=c_char), intent(out) :: message(*)
        integer(c_int), value, intent(in) :: capacity

        character(len=length) :: input_path
        character(:), allocatable :: reason
        integer :: i, iostat

        do i = 1, length
            input_path(i:i) = filename(i)
        end do

        call read_namelist(input_path, iostat)
        if (iostat /= 0) then
            reason = 'cannot read the EMSES namelist (missing file, malformed group or unknown key)'
        else
            call check_emses_input(ispec, regular_pdf /= 0, reason)
        end if

        status = merge(1, 0, len(reason) > 0)
        do i = 1, max(0, capacity)
            message(i) = c_null_char
        end do
        do i = 1, min(len(reason), max(0, capacity - 1))
            message(i) = reason(i:i)
        end do
    end function

    subroutine check_emses_input(ispec, regular_pdf, message)
        !! Check the namelist state in `m_allcom` for species `ispec`.
        integer, intent(in) :: ispec
        logical, intent(in) :: regular_pdf
        character(:), allocatable, intent(out) :: message

        message = ''
        if (ispec < 1 .or. ispec > min(nspec, max_nspec)) then
            message = 'ispec is outside the EMSES species range'
        else if (min(nx, ny, nz) < 1) then
            message = 'EMSES grid dimensions nx, ny, nz must be positive'
        else if (any(nepl < 0) .or. sum(nepl) > max_nepl) then
            message = 'invalid number of emission patches (nepl)'
        else if (nflag_emit(ispec) < 0 .or. nflag_emit(ispec) > 2) then
            message = 'nflag_emit must be 0, 1 or 2; secondary electron emission (3) cannot be replayed'
        else if (.not. all(ieee_is_finite([path(ispec), peth(ispec), spa(ispec), spe(ispec), speth(ispec), &
                                           abvdem(ispec), vdx(ispec), vdy(ispec), vdz(ispec), vdri(ispec), &
                                           vdthz(ispec), vdthxy(ispec), curf(ispec)]))) then
            message = 'nonfinite source parameter'
        else if (min(path(ispec), peth(ispec)) < 0d0) then
            message = 'thermal standard deviations path and peth must be nonnegative'
        else if (regular_pdf .and. min(path(ispec), peth(ispec)) == 0d0) then
            message = 'a cold beam (path or peth = 0) has no regular probability density'
        end if
        if (len(message) > 0) return

        call check_patches(ispec, message)
        if (len(message) > 0) return
        call check_geometry(message)
    end subroutine

    subroutine check_patches(ispec, message)
        integer, intent(in) :: ispec
        character(:), allocatable, intent(inout) :: message

        integer :: iepl, first

        if (nepl(ispec) == 0) return
        if (vdri(ispec) /= 0d0) then
            ! EMSES adds the vdri bulk velocity after sampling the outward flux,
            ! which shifts the half-space; vdist's patch source cannot do that.
            message = 'patch emission with a nonzero vdri bulk drift is not supported'
            return
        end if

        first = sum(nepl(1:ispec - 1)) + 1
        do iepl = first, first + nepl(ispec) - 1
            if (abs(nemd(iepl)) < 1 .or. abs(nemd(iepl)) > 3) then
                message = 'nemd of an emission patch must be +/-1, +/-2 or +/-3'
                return
            end if
            if (.not. (curfs(iepl) > 0d0 .or. curf(ispec) > 0d0)) then
                ! EMSES then uses flpf/flpfs or the npin density, which vdist does not read.
                message = 'patch emission needs curfs > 0 or curf > 0'
                return
            end if
        end do
    end subroutine

    subroutine check_geometry(message)
        !! Mirror the fail-fast checks of the EMSES geometry builders.
        character(:), allocatable, intent(inout) :: message

        integer :: itype, ntypes
        character(len=30) :: name
        double precision :: width(3), edge1(3), edge2(3), normal(3)

        if (.not. all(ieee_is_finite([rectangle_shape, rectangle_corners, sphere_origin, sphere_radius, &
                                      circle_origin, circle_radius, cuboid_shape, cylinder_origin, &
                                      cylinder_radius, cylinder_height, disk_origin, disk_height, &
                                      disk_radius, disk_inner_radius, plane_with_circle_origin, &
                                      plane_with_circle_radius, plane_with_circle_hole_zlower, &
                                      plane_with_circle_hole_height, plane_with_circle_hole_radius, &
                                      zssurf, rcurv, xlrechole, xurechole, ylrechole, yurechole, &
                                      zlrechole, zurechole]))) then
            message = 'nonfinite geometry parameter'
            return
        end if

        ntypes = 1
        if (boundary_type == 'complex') ntypes = nboundary_types
        do itype = 1, ntypes
            name = boundary_type
            if (boundary_type == 'complex') name = boundary_types(itype)
            if (is_collision_boundary_entry_empty(name)) cycle
            if (.not. is_collision_boundary_type_supported(name)) then
                message = 'unsupported boundary type: '//trim(name)
                return
            end if
            if (boundary_type /= 'complex' .and. is_collision_object_boundary_type(name)) then
                message = 'object boundary type '//trim(name)//' requires boundary_type = complex'
                return
            end if

            select case (trim(name))
            case ('sphere')
                if (sphere_radius(itype) <= 0d0) message = 'sphere_radius must be positive'
            case ('circlex', 'circley', 'circlez')
                if (circle_radius(itype) <= 0d0) message = 'circle_radius must be positive'
            case ('cylinderx', 'cylindery', 'cylinderz', 'open-cylinderx', 'open-cylindery', 'open-cylinderz')
                if (cylinder_radius(itype) <= 0d0 .or. cylinder_height(itype) <= 0d0) &
                    message = 'cylinder_radius and cylinder_height must be positive'
            case ('diskx', 'disky', 'diskz')
                if (disk_radius(itype) <= 0d0 .or. disk_height(itype) <= 0d0 .or. disk_inner_radius(itype) < 0d0 &
                    .or. disk_inner_radius(itype) >= disk_radius(itype)) &
                    message = 'disk needs radius > inner radius >= 0 and a positive height'
            case ('plane-with-circlex', 'plane-with-circley', 'plane-with-circlez')
                if (plane_with_circle_radius(itype) <= 0d0) message = 'plane_with_circle_radius must be positive'
            case ('rectangle')
                width = rectangle_shape(2:6:2, itype) - rectangle_shape(1:5:2, itype)
                if (any(width < 0d0) .or. count(width == 0d0) /= 1) &
                    message = 'rectangle_shape must span two axes and be flat along one'
            case ('cuboid')
                width = cuboid_shape(2:6:2, itype) - cuboid_shape(1:5:2, itype)
                if (any(width <= 0d0)) message = 'cuboid_shape extents must be positive'
            case ('rectangle-corners')
                edge1 = rectangle_corners(:, 2, itype) - rectangle_corners(:, 1, itype)
                edge2 = rectangle_corners(:, 4, itype) - rectangle_corners(:, 1, itype)
                normal = [edge1(2)*edge2(3) - edge1(3)*edge2(2), edge1(3)*edge2(1) - edge1(1)*edge2(3), &
                          edge1(1)*edge2(2) - edge1(2)*edge2(1)]
                if (norm2(normal) <= tiny(1d0)) message = 'rectangle_corners are degenerate'
            case ('plane-with-circle-hole')
                if (plane_with_circle_hole_radius(itype) <= 0d0 .or. plane_with_circle_hole_height(itype) <= 0d0) &
                    message = 'plane_with_circle_hole radius and height must be positive'
            case ('rectangle-hole', 'dent-rectangle-hole', 'cylinder-hole', 'hyperboloid-hole', 'ellipsoid-hole')
                width = [xurechole(1) - xlrechole(1), yurechole(1) - ylrechole(1), zssurf - zlrechole(2)]
                if (any(width <= 0d0)) then
                    message = 'hole extents (x/y/zlrechole to x/yurechole, zssurf) must be positive'
                else if ((name == 'dent-rectangle-hole' .or. name == 'hyperboloid-hole' .or. &
                          name == 'ellipsoid-hole') .and. rcurv <= 0d0) then
                    message = 'rcurv must be positive'
                end if
            end select
            if (len(message) > 0) return
        end do
    end subroutine

end module

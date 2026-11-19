program test_emses_input_validation
    !! `check_emses_input` reports settings that vdist cannot replay or that the
    !! EMSES geometry builders would stop on, instead of failing silently or
    !! aborting the Python process.
    use m_allcom
    use m_emses_input_validation, only: check_emses_input
    use m_test_helpers, only: assert_true

    implicit none

    call test_valid_baseline()
    call test_source_rules()
    call test_patch_rules()
    call test_geometry_rules()

    print *, "test_emses_input_validation: all tests passed."

contains

    subroutine baseline()
        call reset_input_parameters()
        nx = 8; ny = 8; nz = 8
        nspec = 1
        path(1) = 1d0
        peth(1) = 1d0
        boundary_type = 'none'
    end subroutine

    subroutine expect(label, ispec, regular_pdf, valid)
        character(len=*), intent(in) :: label
        integer, intent(in) :: ispec
        logical, intent(in) :: regular_pdf, valid
        character(:), allocatable :: message

        call check_emses_input(ispec, regular_pdf, message)
        if (valid .and. len(message) > 0) print *, 'unexpected message: ', message
        call assert_true(label, (len(message) == 0) .eqv. valid)
    end subroutine

    subroutine test_valid_baseline()
        call baseline()
        call expect('Maxwellian outer source is valid', 1, .true., .true.)
        boundary_type = 'complex'
        boundary_types(1) = 'sphere'
        sphere_radius(1) = 1d0
        call expect('complex sphere is valid', 1, .true., .true.)
    end subroutine

    subroutine test_source_rules()
        call baseline()
        call expect('species outside nspec', 2, .true., .false.)
        nflag_emit(1) = 3
        call expect('secondary electron emission is rejected', 1, .true., .false.)
        call baseline()
        peth(1) = 0d0
        call expect('cold beam has no regular density', 1, .true., .false.)
        call expect('cold beam is allowed for range envelopes', 1, .false., .true.)
        peth(1) = -1d0
        call expect('negative thermal width', 1, .false., .false.)
        call baseline()
        nz = 0
        call expect('nonpositive grid', 1, .true., .false.)
    end subroutine

    subroutine test_patch_rules()
        call baseline()
        nflag_emit(1) = 1
        nepl(1) = 1
        nemd(1) = -2
        curf(1) = 1d0
        call expect('curf patch is valid', 1, .true., .true.)
        curf(1) = -1d0
        call expect('patch without positive curf/curfs is rejected', 1, .true., .false.)
        curfs(1) = 0.5d0
        call expect('positive curfs is enough', 1, .true., .true.)
        vdri(1) = 0.3d0
        call expect('bulk-shifted patch emission is rejected', 1, .true., .false.)
        vdri(1) = 0d0
        nemd(1) = 0
        call expect('patch normal must be an axis', 1, .true., .false.)
    end subroutine

    subroutine test_geometry_rules()
        call baseline()
        boundary_type = 'typo'
        call expect('unknown boundary type', 1, .true., .false.)
        boundary_type = 'sphere'
        call expect('object type outside complex mode', 1, .true., .false.)
        call baseline()
        boundary_type = 'complex'
        boundary_types(1) = 'diskz'
        disk_radius(1) = 1d0
        disk_inner_radius(1) = 1d0
        disk_height(1) = 0.5d0
        call expect('disk inner radius must be smaller', 1, .true., .false.)
        disk_inner_radius(1) = 0.5d0
        call expect('annular disk is valid', 1, .true., .true.)
        boundary_types(1) = 'rectangle'
        rectangle_shape(:, 1) = [1d0, 2d0, 1d0, 2d0, 1d0, 2d0]
        call expect('rectangle needs one flat axis', 1, .true., .false.)
        rectangle_shape(:, 1) = [1d0, 2d0, 1d0, 2d0, 3d0, 3d0]
        call expect('flat rectangle is valid', 1, .true., .true.)
        call baseline()
        boundary_type = 'rectangle-hole'
        zssurf = 4d0
        xlrechole(1) = 0d0; xurechole(1) = 3d0
        ylrechole(1) = 2d0; yurechole(1) = 3d0
        zlrechole(2) = 1d0
        call expect('hole touching the domain edge is valid (EMSES pads the caps)', 1, .true., .true.)
        zlrechole(2) = 4d0
        call expect('hole depth must be positive', 1, .true., .false.)
    end subroutine

end program

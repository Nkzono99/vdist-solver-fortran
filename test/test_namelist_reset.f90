program test_namelist_reset
    !! Regression for issue #40: values omitted from a namelist file must not
    !! survive from the previous file read in the same process.
    use, intrinsic :: iso_c_binding, only: c_char
    use m_allcom
    use m_namelist, only: read_namelist
    use m_emses_simulator_builder, only: create_simulator, destroy_simulator
    use m_simulator, only: t_ESSimulator
    use m_test_helpers, only: assert_close, assert_close_vec, assert_equal_int, assert_true

    implicit none

    character(len=*), parameter :: minimal_file = 'test-namelist-minimal.inp'
    character(len=*), parameter :: full_file = 'test-namelist-full.inp'

    call write_minimal(minimal_file, extra='')
    call write_full(full_file)

    call test_omitted_values_return_to_defaults()
    call test_simulator_source_is_independent_of_call_order()

    call delete_file(minimal_file)
    call delete_file(full_file)
    print *, "test_namelist_reset: all tests passed."

contains

    subroutine write_minimal(filename, extra)
        character(len=*), intent(in) :: filename, extra
        integer :: unit

        open (newunit=unit, file=filename, status='replace')
        write (unit, *) '&esorem emflag=0 /'
        write (unit, *) '&plasma wp=10*0, wc=0, phiz=0, phixy=0 /'
        write (unit, *) '&tmgrid nx=2, ny=1, nz=1, dt=0.1 /'
        write (unit, *) '&system nspec=1, npbnd=2,0,0 /'
        write (unit, *) '&intp qm=10*1, peth=1, path=1, vdri=10*0, vdthz=10*0, vdthxy=10*0 '//trim(extra)//' /'
        write (unit, *) "&ptcond boundary_type='none' /"
        write (unit, *) '&emissn use_raycast=.false., nflag_emit=10*0, nepl=10*0 /'
        close (unit)
    end subroutine

    subroutine write_full(filename)
        !! Sets a non-default value in every namelist group.
        character(len=*), intent(in) :: filename
        integer :: unit

        open (newunit=unit, file=filename, status='replace')
        write (unit, *) '&esorem emflag=1 /'
        write (unit, *) '&plasma wp=10*2, wc=3, phiz=45, phixy=31 /'
        write (unit, *) '&tmgrid nx=8, ny=8, nz=4, dt=0.2 /'
        write (unit, *) '&system nspec=2, npbnd=0,0,2,27*2 /'
        write (unit, *) '&intp qm=10*-1, peth=10*4, path=10*4, vdri=10*1, vdthz=10*30, vdthxy=10*60,', &
            ' spa=10*3, spe=10*2, speth=10*90 /'
        write (unit, *) "&ptcond boundary_type='complex', zssurf=1.5, boundary_types='rectangle',9*'none',", &
            ' rectangle_shape(:,1)=0.5,0.5,2,6,0,4, boundary_conductor_id=10*1, rcurv=0.5,', &
            ' cylinder_origin=30*1, cylinder_radius=10*1, cylinder_height=10*1,', &
            ' sphere_origin=30*1, sphere_radius=10*1, circle_origin=30*1, circle_radius=10*1,', &
            ' cuboid_shape=60*1, disk_origin=30*1, disk_height=10*1, disk_radius=10*1,', &
            ' disk_inner_radius=10*1, conductivity=2, max_bounce_count=7,', &
            ' plane_with_circle_hole_zlower=10*1, plane_with_circle_hole_height=10*1,', &
            ' plane_with_circle_hole_radius=10*1, plane_with_circle_origin=30*1, plane_with_circle_radius=10*1,', &
            ' xlrechole=1,2, ylrechole=1,2, zlrechole=1,2, xurechole=1,2, yurechole=1,2, zurechole=1,2,', &
            ' boundary_mirror_reflection_rate=100*0.5, boundary_reversal_reflection_rate=100*0.5,', &
            ' boundary_mirror_reflect_alpha=100*0.5, boundary_reversal_reflect_alpha=100*0.5,', &
            ' boundary_mirror_reflect_energy_loss_frac=100*0.5, boundary_reversal_reflect_energy_loss_frac=100*0.5,', &
            ' enable_secondary_electron_emission=100*.true., boundary_se_model_type=100*2,', &
            ' boundary_se_const_yield=100*0.5, boundary_se_yield_max=100*0.5, boundary_se_energy_max=100*5,', &
            ' boundary_se_species_id=100*2 /'
        write (unit, *) '&emissn nflag_emit=10*1, nepl=1,9*0, nemd=15*-3, curf=10*1, curfs=15*2,', &
            ' xmine=15*1, xmaxe=15*3, ymine=15*1, ymaxe=15*3, zmine=15*2, zmaxe=15*2,', &
            ' thetaz=15*10, thetaxy=15*20, use_raycast=.true.,', &
            ' ray_zenith_angle_deg=10*180, ray_azimuth_angle_deg=10*90 /'
        close (unit)
    end subroutine

    subroutine delete_file(filename)
        character(len=*), intent(in) :: filename
        integer :: unit

        open (newunit=unit, file=filename, status='old')
        close (unit, status='delete')
    end subroutine

    subroutine test_omitted_values_return_to_defaults()
        call read_namelist(full_file)
        call assert_close("full input sets spa", spa(1), 3d0)
        call assert_true("full input sets raycast", use_raycast)
        call assert_equal_int("full input sets nepl", nepl(1), 1)

        call read_namelist(minimal_file)
        call assert_equal_int("emflag returns to default", emflag, 0)
        call assert_close_vec("npbnd returns to the minimal input", dble(npbnd(:, 1)), [2d0, 0d0, 0d0])
        call assert_equal_int("nspec returns to the minimal input", nspec, 1)
        call assert_equal_int("omitted second-species npbnd returns to default", npbnd(1, 2), 0)
        call assert_close("omitted spa returns to default", spa(1), 0d0)
        call assert_close("omitted spe returns to default", spe(1), 0d0)
        call assert_close("omitted speth returns to default", speth(1), 0d0)
        call assert_close("omitted vdri is read from the minimal input", vdri(1), 0d0)
        call assert_close("omitted zssurf returns to default", zssurf, -9999d0)
        call assert_true("omitted boundary_types return to default", all(boundary_types == "none"))
        call assert_equal_int("omitted boundary_conductor_id returns to default", boundary_conductor_id(1), -1)
        call assert_close("omitted rectangle_shape returns to default", maxval(abs(rectangle_shape)), 0d0)
        call assert_close("omitted rcurv returns to default", rcurv, 1d0)
        call assert_close("omitted conductivity returns to default", conductivity, 0d0)
        call assert_equal_int("omitted max_bounce_count returns to default", max_bounce_count, 3)
        call assert_close("omitted reflection rates return to default", &
                          maxval(abs(boundary_mirror_reflection_rate)), 0d0)
        call assert_true("omitted secondary emission flag returns to default", &
                         .not. any(enable_secondary_electron_emission))
        call assert_close("omitted SEE peak energy returns to default", boundary_se_energy_max(1, 1), 300d0)
        call assert_equal_int("omitted SEE species returns to default", boundary_se_species_id(1, 1), 1)
        call assert_close("omitted hole coordinates return to default", maxval(abs(xlrechole)), 0d0)
        call assert_equal_int("omitted nemd returns to default", nemd(1), 0)
        call assert_close("omitted curf returns to default", curf(1), 0d0)
        call assert_close("omitted curfs returns to default", curfs(1), -9999d0)
        call assert_close("omitted emission box returns to default", maxval(abs([xmine, xmaxe, ymine, ymaxe, zmine, zmaxe])), 0d0)
        call assert_close("omitted emission angles return to default", maxval(abs([thetaz, thetaxy])), 0d0)
        call assert_true("omitted use_raycast returns to default", .not. use_raycast)
        call assert_close("omitted ray zenith returns to sentinel", ray_zenith_angle_deg(1), 9999d0)
        call assert_close("omitted ray azimuth returns to sentinel", ray_azimuth_angle_deg(1), 9999d0)
    end subroutine

    function source_pdf(filename, extra) result(pdf)
        !! Build a simulator from `filename` (rewritten with `extra` intp keys)
        !! and evaluate the external Maxwellian at v = (-1, 0, 0).
        character(len=*), intent(in) :: filename, extra
        double precision :: pdf
        type(t_ESSimulator) :: simulator
        character(kind=c_char) :: cpath(len(filename))
        double precision :: fields(9, 3, 2, 2)
        integer :: i

        call write_minimal(filename, extra)
        do i = 1, len(filename)
            cpath(i) = filename(i:i)
        end do
        fields = 0d0
        simulator = create_simulator(cpath, len(filename), 2, 1, 1, fields, 1, 10)
        pdf = simulator%probability_functions(3)%at([2d0, 0.5d0, 0.5d0], [-1d0, 0d0, 0d0])
        call destroy_simulator(simulator)
    end function

    subroutine test_simulator_source_is_independent_of_call_order()
        double precision :: expected, fresh, after_drift, explicit_zero

        expected = exp(-0.5d0)/(2d0*pi)**1.5d0
        fresh = source_pdf(minimal_file, '')
        call assert_close("fresh input B matches the drift-free Maxwellian", fresh, expected)
        call assert_true("input A with spa=3 carries a drift", abs(source_pdf(minimal_file, ', spa=3') - expected) > 1d-3)
        after_drift = source_pdf(minimal_file, '')
        call assert_close("input B after input A is unchanged", after_drift, expected)
        explicit_zero = source_pdf(minimal_file, ', spa=0')
        call assert_close("explicit spa=0 agrees with the omitted key", explicit_zero, expected)
    end subroutine

end program

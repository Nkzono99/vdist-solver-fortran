program test_field
    use m_field, only: t_VectorFieldGrid, new_VectorFieldGrid
    use m_test_helpers, only: assert_close, assert_close_vec

    implicit none

    call test_uniform_field_returns_constant()
    call test_corner_position_returns_corner_value()
    call test_linear_gradient_in_x()
    call test_trilinear_interpolation_center()
    call test_position_outside_domain_linearly_extrapolates()
    call test_component_offset_shifts_interpolation_position()
    call test_shared_stencil_matches_reference_on_nonaffine_field()
    call test_singleton_axis_is_constant()

    print *, "test_field: all tests passed."

contains

    subroutine test_uniform_field_returns_constant()
        type(t_VectorFieldGrid) :: grid
        double precision :: values(2, 0:1, 0:1, 0:1)

        values = 0d0
        values(1, :, :, :) = 3d0
        values(2, :, :, :) = -5d0

        grid = new_VectorFieldGrid(2, 1, 1, 1, values)

        call assert_close_vec("Uniform field at (0.3, 0.7, 0.2)", &
                              grid%at([0.3d0, 0.7d0, 0.2d0]), [3d0, -5d0])
    end subroutine

    subroutine test_corner_position_returns_corner_value()
        type(t_VectorFieldGrid) :: grid
        double precision :: values(1, 0:1, 0:1, 0:1)

        values(1, :, :, :) = 0d0
        values(1, 1, 1, 1) = 9d0

        grid = new_VectorFieldGrid(1, 1, 1, 1, values)

        call assert_close_vec("grid%at at corner (1,1,1)", &
                              grid%at([1d0, 1d0, 1d0]), [9d0])
        call assert_close_vec("grid%at at corner (0,0,0)", &
                              grid%at([0d0, 0d0, 0d0]), [0d0])
    end subroutine

    subroutine test_linear_gradient_in_x()
        type(t_VectorFieldGrid) :: grid
        double precision :: values(1, 0:2, 0:1, 0:1)
        integer :: ix, iy, iz

        do ix = 0, 2
            do iy = 0, 1
                do iz = 0, 1
                    values(1, ix, iy, iz) = dble(ix)
                end do
            end do
        end do

        grid = new_VectorFieldGrid(1, 2, 1, 1, values)

        call assert_close_vec("Linear-in-x at x=0.25", grid%at([0.25d0, 0.5d0, 0.5d0]), [0.25d0])
        call assert_close_vec("Linear-in-x at x=1.75", grid%at([1.75d0, 0.5d0, 0.5d0]), [1.75d0])
    end subroutine

    subroutine test_trilinear_interpolation_center()
        !! Set the eight corners of a single cell to 0..7 and query the
        !! centroid — the expected value is the mean of the eight corners.
        type(t_VectorFieldGrid) :: grid
        double precision :: values(1, 0:1, 0:1, 0:1)
        integer :: k, iz, iy, ix
        double precision :: mean_value

        k = 0
        do iz = 0, 1
            do iy = 0, 1
                do ix = 0, 1
                    values(1, ix, iy, iz) = dble(k)
                    k = k + 1
                end do
            end do
        end do

        mean_value = sum(values(1, :, :, :))/8d0

        grid = new_VectorFieldGrid(1, 1, 1, 1, values)

        call assert_close_vec("Trilinear at (0.5, 0.5, 0.5) = mean of corners", &
                              grid%at([0.5d0, 0.5d0, 0.5d0]), [mean_value])
    end subroutine

    subroutine test_position_outside_domain_linearly_extrapolates()
        !! Beyond the grid, `at` clamps the cell index (ip) but leaves the
        !! intra-cell parameter rp untouched. That yields linear extrapolation
        !! along the nearest edge cell — for a field stepping 0 -> 4 between
        !! x=0 and x=1 at a 1-cell grid, querying x=2 yields 8 (rp=2).
        type(t_VectorFieldGrid) :: grid
        double precision :: values(1, 0:1, 0:1, 0:1)

        values(1, :, :, :) = 0d0
        values(1, 1, 0, 0) = 4d0
        values(1, 1, 1, 0) = 4d0
        values(1, 1, 0, 1) = 4d0
        values(1, 1, 1, 1) = 4d0

        grid = new_VectorFieldGrid(1, 1, 1, 1, values)

        call assert_close_vec("Position past +X linearly extrapolates", &
                              grid%at([2d0, 0.3d0, 0.7d0]), [8d0])
    end subroutine

    subroutine test_component_offset_shifts_interpolation_position()
        type(t_VectorFieldGrid) :: grid
        double precision :: values(1, 0:2, 0:1, 0:1)
        double precision :: offsets(3, 1)
        integer :: ix, iy, iz

        do ix = 0, 2
            do iy = 0, 1
                do iz = 0, 1
                    values(1, ix, iy, iz) = dble(ix)
                end do
            end do
        end do

        offsets(:, :) = 0d0
        offsets(1, 1) = 0.5d0

        grid = new_VectorFieldGrid(1, 2, 1, 1, values, offsets)

        call assert_close_vec("Offset component samples x - 0.5", &
                              grid%at([0.75d0, 0.5d0, 0.5d0]), [0.25d0])
    end subroutine

    pure function reference_trilinear(values, nx, ny, nz, component, point) result(ret)
        !! Independent per-component trilinear formula with the documented
        !! edge rule: the lower node is clamped to [0, n-1]; singleton axes
        !! (n = 0) are constant.
        integer, intent(in) :: nx, ny, nz, component
        double precision, intent(in) :: values(:, 0:, 0:, 0:)
        double precision, intent(in) :: point(3)
        double precision :: ret

        integer :: lower(3), upper(3), n(3), a, corner, c(3)
        double precision :: t(3), weight

        n = [nx, ny, nz]
        do a = 1, 3
            if (n(a) == 0) then
                lower(a) = 0; upper(a) = 0; t(a) = 0d0
            else
                lower(a) = max(0, min(floor(point(a)), n(a) - 1))
                upper(a) = lower(a) + 1
                t(a) = point(a) - lower(a)
            end if
        end do
        ret = 0d0
        do corner = 0, 7
            weight = 1d0
            do a = 1, 3
                if (btest(corner, a - 1)) then
                    c(a) = upper(a); weight = weight*t(a)
                else
                    c(a) = lower(a); weight = weight*(1d0 - t(a))
                end if
            end do
            ret = ret + weight*values(component, c(1), c(2), c(3))
        end do
    end function

    subroutine test_shared_stencil_matches_reference_on_nonaffine_field()
        type(t_VectorFieldGrid) :: grid
        double precision :: values(6, 0:2, 0:3, 0:1), point(3), expected(6)
        integer :: i, j, k, component, probe

        do k = 0, 1
            do j = 0, 3
                do i = 0, 2
                    do component = 1, 6
                        values(component, i, j, k) = sin(dble(7*i + 3*j + 11*k + component))
                    end do
                end do
            end do
        end do
        grid = new_VectorFieldGrid(6, 2, 3, 1, values)
        do probe = -4, 12
            point = [0.25d0*probe, 0.31d0*probe, 0.125d0*probe]
            do component = 1, 6
                expected(component) = reference_trilinear(values, 2, 3, 1, component, point)
            end do
            call assert_close_vec("Nodal interior and extrapolated points", grid%at(point), expected, 1d-13)
        end do
    end subroutine

    subroutine test_singleton_axis_is_constant()
        !! An axis with a single node must not read index -1 (bounds-checked
        !! builds abort) and must be constant along that axis.
        type(t_VectorFieldGrid) :: grid
        double precision :: values(2, 0:0, 0:3, 0:1), point(3), expected(2), offsets(3, 2)
        integer :: j, k, component, probe

        do k = 0, 1
            do j = 0, 3
                do component = 1, 2
                    values(component, 0, j, k) = cos(dble(5*j + 13*k + component))
                end do
            end do
        end do
        offsets = 0d0
        offsets(2, 2) = 0.5d0
        grid = new_VectorFieldGrid(2, 0, 3, 1, values, offsets)
        do probe = -4, 12
            point = [0.4d0*probe, 0.31d0*probe, 0.125d0*probe]
            expected(1) = reference_trilinear(values, 0, 3, 1, 1, point)
            expected(2) = reference_trilinear(values, 0, 3, 1, 2, point - offsets(:, 2))
            call assert_close_vec("Singleton x axis", grid%at(point), expected, 1d-13)
        end do
    end subroutine

end program

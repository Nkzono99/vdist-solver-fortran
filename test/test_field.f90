program test_field
    use m_field, only: t_VectorFieldGrid, new_VectorFieldGrid
    use m_test_helpers, only: assert_close, assert_close_vec

    implicit none

    call test_uniform_field_returns_constant()
    call test_corner_position_returns_corner_value()
    call test_linear_gradient_in_x()
    call test_trilinear_interpolation_center()
    call test_position_outside_domain_linearly_extrapolates()

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

end program

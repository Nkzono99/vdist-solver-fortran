program test_emses_field
    use m_emses_field, only: t_EMSESFieldGrid, new_EMSESFieldGrid
    use m_test_helpers, only: assert_close_vec

    implicit none

    call test_accumulated_e_uses_component_offsets()

    print *, "test_emses_field: all tests passed."

contains

    subroutine test_accumulated_e_uses_component_offsets()
        type(t_EMSESFieldGrid) :: grid
        double precision :: relocated_eb(6, 0:2, 0:1, 0:1)
        double precision :: accumulated_e(3, 0:2, 0:1, 0:1)
        integer :: ix, iy, iz

        relocated_eb(:, :, :, :) = 0d0
        relocated_eb(1, :, :, :) = 10d0
        relocated_eb(4, :, :, :) = 4d0

        accumulated_e(:, :, :, :) = 0d0
        do ix = 0, 2
            do iy = 0, 1
                do iz = 0, 1
                    accumulated_e(1, ix, iy, iz) = dble(ix)
                    accumulated_e(2, ix, iy, iz) = 2d0
                    accumulated_e(3, ix, iy, iz) = 3d0
                end do
            end do
        end do

        grid = new_EMSESFieldGrid(2, 1, 1, relocated_eb, accumulated_e)

        call assert_close_vec("EMSES field adds shifted accumulated E", &
                              grid%at([0.75d0, 0.5d0, 0.5d0]), &
                              [10.25d0, 2d0, 3d0, 4d0, 0d0, 0d0])
    end subroutine

end program

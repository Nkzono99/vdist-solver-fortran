program test_particle
    use m_particle, only: t_Particle, new_Particle
    use m_test_helpers, only: assert_close, assert_close_vec

    implicit none

    call test_new_Particle_sets_all_fields()
    call test_new_Particle_zero_initial_time()

    print *, "test_particle: all tests passed."

contains

    subroutine test_new_Particle_sets_all_fields()
        type(t_Particle) :: p

        p = new_Particle(2d0, [1d0, 2d0, 3d0], [4d0, 5d0, 6d0])

        call assert_close("new_Particle sets q_m", p%q_m, 2d0)
        call assert_close_vec("new_Particle sets position", p%position, [1d0, 2d0, 3d0])
        call assert_close_vec("new_Particle sets velocity", p%velocity, [4d0, 5d0, 6d0])
    end subroutine

    subroutine test_new_Particle_zero_initial_time()
        type(t_Particle) :: p

        p = new_Particle(1d0, [0d0, 0d0, 0d0], [0d0, 0d0, 0d0])

        call assert_close("new_Particle starts at t = 0", p%t, 0d0)
    end subroutine

end program

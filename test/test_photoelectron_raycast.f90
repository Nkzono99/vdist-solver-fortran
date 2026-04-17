program test_photoelectron_raycast
    !! Unit tests for the raycast-based photoelectron probability.
    !!
    !! The probability equals the Maxwellian PDF of the velocity when the
    !! sun ray from the collision point escapes without hitting any blocker,
    !! and zero when an internal boundary occludes the sun.
    use m_photoelectron_raycast, only: t_PhotoelectronRaycastProbability, &
                                       new_PhotoelectronRaycastProbability
    use finbound, only: t_BoundaryList, new_BoundaryList, &
                        t_Boundary, t_PlaneXYZ, new_PlaneZ
    use m_test_helpers, only: assert_close, assert_true

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)

    call test_unoccluded_ray_returns_maxwell_pdf()
    call test_occluded_ray_returns_zero()
    call test_zero_coefficient_stays_zero_even_unoccluded()

    print *, "test_photoelectron_raycast: all tests passed."

contains

    function empty_boundary_list() result(blist)
        type(t_BoundaryList) :: blist

        blist = new_BoundaryList()
    end function

    function single_plane_z_blocker(z0) result(blist)
        !! Boundary list with one planeZ occluder at z = z0.
        double precision, intent(in) :: z0
        type(t_BoundaryList) :: blist

        class(t_Boundary), pointer :: pb
        type(t_PlaneXYZ), pointer :: pplane

        blist = new_BoundaryList()

        allocate (pplane)
        pplane = new_PlaneZ(z0)
        pb => pplane
        call blist%add_boundary(pb)
    end function

    function maxwell_pdf_3d(v, loc, sigma) result(ret)
        double precision, intent(in) :: v(3)
        double precision, intent(in) :: loc(3)
        double precision, intent(in) :: sigma(3)
        double precision :: ret
        integer :: i

        ret = 1d0
        do i = 1, 3
            ret = ret*(1d0/sqrt(2d0*pi*sigma(i)*sigma(i)))* &
                  exp(-(v(i) - loc(i))*(v(i) - loc(i))/(2d0*sigma(i)*sigma(i)))
        end do
    end function

    subroutine test_unoccluded_ray_returns_maxwell_pdf()
        !! No blockers -> ray escapes -> probability = Maxwell PDF at velocity.
        type(t_PhotoelectronRaycastProbability) :: prob
        double precision :: loc(3), sigma(3), sun_dir(3)
        double precision :: position(3), velocity(3)
        double precision :: expected

        loc = [0d0, 0d0, 0d0]
        sigma = [1d0, 1d0, 1d0]
        sun_dir = [0d0, 0d0, 1d0]  ! sun straight overhead
        position = [0.5d0, 0.5d0, 0.5d0]
        velocity = [0.2d0, -0.1d0, 0.3d0]

        prob = new_PhotoelectronRaycastProbability(loc, sigma, sun_dir, &
                                                   empty_boundary_list())

        expected = maxwell_pdf_3d(velocity, loc, sigma)
        call assert_close("unoccluded ray returns Maxwell PDF", &
                          prob%at(position, velocity), expected)
    end subroutine

    subroutine test_occluded_ray_returns_zero()
        !! A plane at z=1 in the sun direction occludes the ray -> 0.
        type(t_PhotoelectronRaycastProbability) :: prob
        double precision :: loc(3), sigma(3), sun_dir(3)
        double precision :: position(3), velocity(3)

        loc = [0d0, 0d0, 0d0]
        sigma = [1d0, 1d0, 1d0]
        sun_dir = [0d0, 0d0, 1d0]
        position = [0.5d0, 0.5d0, 0.5d0]
        velocity = [0.2d0, -0.1d0, 0.3d0]

        prob = new_PhotoelectronRaycastProbability(loc, sigma, sun_dir, &
                                                   single_plane_z_blocker(1d0))

        call assert_close("occluded ray returns 0", &
                          prob%at(position, velocity), 0d0)
    end subroutine

    subroutine test_zero_coefficient_stays_zero_even_unoccluded()
        !! coefficient = 0 should yield 0 without attempting a raycast,
        !! which also means an empty blocker list cannot accidentally
        !! produce a non-zero value.
        type(t_PhotoelectronRaycastProbability) :: prob
        double precision :: loc(3), sigma(3), sun_dir(3)

        loc = [0d0, 0d0, 0d0]
        sigma = [1d0, 1d0, 1d0]
        sun_dir = [0d0, 0d0, 1d0]

        prob = new_PhotoelectronRaycastProbability(loc, sigma, sun_dir, &
                                                   empty_boundary_list(), &
                                                   coefficient=0d0)

        call assert_close("coefficient=0 forces probability to 0", &
                          prob%at([0.5d0, 0.5d0, 0.5d0], [0d0, 0d0, 0d0]), 0d0)
    end subroutine

end program

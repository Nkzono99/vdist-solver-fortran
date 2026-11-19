program test_photoelectron_raycast
    !! Unit tests for the EMSES raycast photoelectron source.
    !!
    !! The density at a lit surface point is flux*max(n.sun,0)/rate times the
    !! surface-frame Gaussian, for outgoing velocities only, and zero in shadow
    !! or on the non-emitting EMSES cover (tag -1).
    use m_photoelectron_raycast, only: t_PhotoelectronRaycastProbability, &
                                       new_PhotoelectronRaycastProbability, surface_basis
    use finbound, only: t_BoundaryList, new_BoundaryList, t_Boundary, new_PlaneZ, new_Sphere
    use m_test_helpers, only: assert_close, assert_close_vec, assert_true

    implicit none

    double precision, parameter :: pi = acos(-1.0d0)
    double precision, parameter :: mean(3) = [0.3d0, 0.2d0, -0.1d0]
    double precision, parameter :: sigma(3) = [0.7d0, 0.5d0, 0.4d0]
    double precision, parameter :: flux = 1.7d0

    call test_surface_basis_matches_emses_sampler()
    call test_lit_plane_returns_flux_normalised_density()
    call test_inward_velocity_returns_zero()
    call test_occluded_point_returns_zero()
    call test_cover_and_empty_space_do_not_emit()
    call test_zero_coefficient_returns_zero()
    call test_curved_surface_uses_local_normal_and_projected_flux()

    print *, "test_photoelectron_raycast: all tests passed."

contains

    function plane_list(z0, tag) result(blist)
        double precision, intent(in) :: z0
        integer, intent(in) :: tag
        type(t_BoundaryList) :: blist
        class(t_Boundary), pointer :: boundary

        blist = new_BoundaryList()
        allocate (boundary, source=new_PlaneZ(z0))
        boundary%material%tag = tag
        call blist%add_boundary(boundary)
    end function

    function gaussian(x, mu, s) result(ret)
        double precision, intent(in) :: x, mu, s
        double precision :: ret

        ret = exp(-0.5d0*((x - mu)/s)**2)/(sqrt(2d0*pi)*s)
    end function

    function normal_flux_rate(mu, s) result(rate)
        !! Independent midpoint quadrature of int_0^inf v G(v) dv.
        double precision, intent(in) :: mu, s
        double precision :: rate, x, dx
        integer :: i

        rate = 0d0
        dx = (abs(mu) + 12d0*s)/200000d0
        do i = 1, 200000
            x = (i - 0.5d0)*dx
            rate = rate + x*gaussian(x, mu, s)*dx
        end do
    end function

    function expected_density(local, illumination) result(ret)
        double precision, intent(in) :: local(3), illumination
        double precision :: ret
        integer :: i

        ret = flux*illumination/normal_flux_rate(mean(1), sigma(1))
        do i = 1, 3
            ret = ret*gaussian(local(i), mean(i), sigma(i))
        end do
    end function

    subroutine test_surface_basis_matches_emses_sampler()
        !! raycast_sampler.F90: t1 = cross(ref, n)/|.|, t2 = cross(n, t1) with
        !! ref = +x unless |n_x| >= 0.9, then +y.
        double precision :: basis(3, 3)

        basis = surface_basis([0d0, 0d0, 1d0])
        call assert_close_vec("z normal tangent 1", basis(:, 2), [0d0, -1d0, 0d0])
        call assert_close_vec("z normal tangent 2", basis(:, 3), [1d0, 0d0, 0d0])
        basis = surface_basis([1d0, 0d0, 0d0])
        call assert_close_vec("x normal tangent 1", basis(:, 2), [0d0, 0d0, -1d0])
        call assert_close_vec("x normal tangent 2", basis(:, 3), [0d0, 1d0, 0d0])
    end subroutine

    subroutine test_lit_plane_returns_flux_normalised_density()
        type(t_PhotoelectronRaycastProbability) :: prob
        double precision :: velocity(3)

        prob = new_PhotoelectronRaycastProbability(mean, sigma, flux, [0d0, 0d0, 1d0], plane_list(0d0, 1))
        velocity = [0.2d0, -0.1d0, 0.3d0]
        ! n = +z, t1 = -y, t2 = +x
        call assert_close("lit plane density", prob%at([0.5d0, 0.5d0, 0d0], velocity), &
                          expected_density([0.3d0, 0.1d0, 0.2d0], 1d0), 1d-9)
        call prob%blocking_boundaries%destroy()
    end subroutine

    subroutine test_inward_velocity_returns_zero()
        type(t_PhotoelectronRaycastProbability) :: prob

        prob = new_PhotoelectronRaycastProbability(mean, sigma, flux, [0d0, 0d0, 1d0], plane_list(0d0, 1))
        call assert_close("velocity into the surface returns 0", &
                          prob%at([0.5d0, 0.5d0, 0d0], [0.2d0, -0.1d0, -0.3d0]), 0d0)
        call prob%blocking_boundaries%destroy()
    end subroutine

    subroutine test_occluded_point_returns_zero()
        type(t_PhotoelectronRaycastProbability) :: prob
        class(t_Boundary), pointer :: boundary

        prob = new_PhotoelectronRaycastProbability(mean, sigma, flux, [0d0, 0d0, 1d0], plane_list(0d0, 1))
        allocate (boundary, source=new_PlaneZ(1d0))
        boundary%material%tag = 2
        call prob%blocking_boundaries%add_boundary(boundary)
        call assert_close("shadowed point returns 0", &
                          prob%at([0.5d0, 0.5d0, 0d0], [0.2d0, -0.1d0, 0.3d0]), 0d0)
        call prob%blocking_boundaries%destroy()
    end subroutine

    subroutine test_cover_and_empty_space_do_not_emit()
        type(t_PhotoelectronRaycastProbability) :: prob

        prob = new_PhotoelectronRaycastProbability(mean, sigma, flux, [0d0, 0d0, 1d0], plane_list(0d0, -1))
        call assert_close("EMSES cover (tag -1) does not emit", &
                          prob%at([0.5d0, 0.5d0, 0d0], [0.2d0, -0.1d0, 0.3d0]), 0d0)
        call assert_close("a point off every surface does not emit", &
                          prob%at([0.5d0, 0.5d0, 0.5d0], [0.2d0, -0.1d0, 0.3d0]), 0d0)
        call prob%blocking_boundaries%destroy()
    end subroutine

    subroutine test_zero_coefficient_returns_zero()
        type(t_PhotoelectronRaycastProbability) :: prob

        prob = new_PhotoelectronRaycastProbability(mean, sigma, flux, [0d0, 0d0, 1d0], plane_list(0d0, 1), &
                                                   coefficient=0d0)
        call assert_close("coefficient=0 forces probability to 0", &
                          prob%at([0.5d0, 0.5d0, 0d0], [0.2d0, -0.1d0, 0.3d0]), 0d0)
        call prob%blocking_boundaries%destroy()
    end subroutine

    subroutine test_curved_surface_uses_local_normal_and_projected_flux()
        type(t_PhotoelectronRaycastProbability) :: prob
        type(t_BoundaryList) :: boundaries
        class(t_Boundary), pointer :: boundary
        double precision :: normal(3), tangent(3), position(3), velocity(3)

        boundaries = new_BoundaryList()
        allocate (boundary, source=new_Sphere([2d0, 2d0, 2d0], 1d0))
        boundary%material%tag = 1
        call boundaries%add_boundary(boundary)
        prob = new_PhotoelectronRaycastProbability(mean, sigma, flux, [0d0, 0d0, 1d0], boundaries)

        normal = [0.6d0, 0d0, 0.8d0]
        tangent = [0.8d0, 0d0, -0.6d0]
        position = [2d0, 2d0, 2d0] + normal
        velocity = 0.2d0*normal + 2d0*tangent
        call assert_true("test velocity leaves the surface but points away from the sun", &
                         dot_product(velocity, normal) > 0d0 .and. velocity(3) < 0d0)
        ! n = (0.6, 0, 0.8): t1 = -y, t2 = (0.8, 0, -0.6); illumination n.sun = 0.8
        call assert_close("curved surface uses the hit normal and projected flux", &
                          prob%at(position, velocity), expected_density([0.2d0, 0d0, 2d0], 0.8d0), 1d-9)
        call assert_close("unlit hemisphere does not emit", &
                          prob%at([2d0, 2d0, 2d0] - normal, -velocity), 0d0)
        call prob%blocking_boundaries%destroy()
    end subroutine

end program

module m_photoelectron_raycast
    !! Raycast photoelectron source matching MPIEMSES3D raycast emission.
    !!
    !! EMSES casts rays from a parallel camera along `-sun_direction`. A ray
    !! that first hits an emitting surface (tag /= -1) emits one particle with
    !! velocity `vn*n + vt1*t1 + vt2*t2`, where `n` is the hit normal, `vn` is
    !! drawn from the one-sided flux distribution of N(mean(1), sigma(1)) and the
    !! tangents from N(mean(2:3), sigma(2:3)) (raycast_sampler.F90). Rays per
    !! projected area are fixed, so the emitted flux per surface area is
    !! `flux*max(n.sun, 0)`.
    !!
    !! The equivalent phase-space density at the surface is
    !!   f(v) = flux*max(n.sun, 0)/rate * G(v_n) G(v_t1) G(v_t2),  v_n > 0,
    !! where `rate` is the one-sided first moment of G(v_n). A point in shadow
    !! (the ray toward the sun hits another boundary) has zero density.

    use m_probabilities, only: t_Probability
    use m_maxwell_flux_erf, only: flux_integral_halfspace_erf
    use finbound, only: t_BoundaryList, t_Ray, new_Ray, t_HitRecord

    implicit none

    private
    public t_PhotoelectronRaycastProbability
    public new_PhotoelectronRaycastProbability
    public surface_basis

    double precision, parameter :: pi = acos(-1.0d0)

    double precision, parameter :: DEFAULT_RAY_ORIGIN_OFFSET = 1d-6
        !! Small offset along the sun direction used to move the ray origin
        !! just off the surface it starts from.

    type, extends(t_Probability) :: t_PhotoelectronRaycastProbability
        double precision :: mean(3)
            !! Mean velocity in the surface frame (normal, tangent 1, tangent 2).
        double precision :: sigma(3)
            !! Standard deviations in the surface frame.
        double precision :: flux
            !! Emitted particle flux of a surface facing the sun.
        double precision :: sun_direction(3)
            !! Unit vector pointing from the surface toward the sun.
        type(t_BoundaryList) :: blocking_boundaries
            !! Emitting and occluding geometry, including the EMSES cover (tag -1).
        double precision :: ray_origin_offset = DEFAULT_RAY_ORIGIN_OFFSET
            !! Offset along `sun_direction` applied to the ray origin.
    contains
        procedure :: at => photoelectronRaycast_at
    end type

contains

    function new_PhotoelectronRaycastProbability(mean, sigma, flux, sun_direction, &
                                                 blocking_boundaries, &
                                                 coefficient, &
                                                 ray_origin_offset) result(obj)
        double precision, intent(in) :: mean(3)
        double precision, intent(in) :: sigma(3)
        double precision, intent(in) :: flux
        double precision, intent(in) :: sun_direction(3)
        type(t_BoundaryList), intent(in) :: blocking_boundaries
        double precision, intent(in), optional :: coefficient
        double precision, intent(in), optional :: ray_origin_offset
        type(t_PhotoelectronRaycastProbability) :: obj

        double precision :: mag

        obj%mean = mean
        obj%sigma = sigma
        obj%flux = flux

        mag = norm2(sun_direction)
        if (mag > 0d0) then
            obj%sun_direction = sun_direction/mag
        else
            obj%sun_direction = [0d0, 0d0, 1d0]
        end if

        obj%blocking_boundaries = blocking_boundaries

        if (present(coefficient)) then
            obj%coefficient = coefficient
        else
            obj%coefficient = 1d0
        end if

        if (present(ray_origin_offset)) then
            obj%ray_origin_offset = ray_origin_offset
        end if
    end function

    function photoelectronRaycast_at(self, position, velocity) result(ret)
        !! Photoelectron phase-space density at a surface collision point.
        class(t_PhotoelectronRaycastProbability), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        double precision :: ret

        type(t_Ray) :: ray
        type(t_HitRecord) :: hit
        double precision :: origin(3), basis(3, 3), local(3)
        double precision :: illumination, rate
        integer :: i

        ret = 0d0
        if (self%coefficient == 0d0 .or. self%flux == 0d0) return

        ! Recover the face the EMSES camera ray would hit at this point: step
        ! toward the sun and cast back along the camera direction.
        origin = position + self%ray_origin_offset*self%sun_direction
        ray = new_Ray(origin, -self%sun_direction)
        hit = self%blocking_boundaries%hit(ray)
        if (.not. hit%is_hit) return
        if (abs(hit%t - self%ray_origin_offset) > 0.1d0*self%ray_origin_offset) return
        if (hit%material%tag == -1) return

        illumination = dot_product(hit%n, self%sun_direction)
        if (illumination <= 0d0) return

        basis = surface_basis(hit%n)
        local = matmul(transpose(basis), velocity)
        if (local(1) <= 0d0) return

        rate = flux_integral_halfspace_erf(sqrt(2d0)*self%sigma(1), self%mean(1))
        if (.not. (rate > 0d0)) return

        ret = self%coefficient*self%flux*illumination/rate
        do i = 1, 3
            ret = ret*gaussian_pdf(local(i), self%mean(i), self%sigma(i))
        end do
        if (ret == 0d0) return

        ray = new_Ray(origin, self%sun_direction)
        hit = self%blocking_boundaries%hit(ray)
        if (hit%is_hit) ret = 0d0
    end function

    pure function surface_basis(normal) result(basis)
        !! Columns: normal and the two tangents of EMSES raycast_sampler.F90.
        double precision, intent(in) :: normal(3)
        double precision :: basis(3, 3)
        double precision :: reference(3)

        if (abs(normal(1)) < 0.9d0) then
            reference = [1d0, 0d0, 0d0]
        else
            reference = [0d0, 1d0, 0d0]
        end if
        basis(:, 1) = normal
        basis(:, 2) = cross(reference, normal)
        basis(:, 2) = basis(:, 2)/norm2(basis(:, 2))
        basis(:, 3) = cross(normal, basis(:, 2))
    end function

    pure function cross(a, b) result(c)
        double precision, intent(in) :: a(3), b(3)
        double precision :: c(3)

        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function

    pure function gaussian_pdf(x, mu, sigma) result(ret)
        !! Normal density; a cold component (sigma <= 0) has no regular density.
        double precision, intent(in) :: x
        double precision, intent(in) :: mu
        double precision, intent(in) :: sigma
        double precision :: ret

        ret = 0d0
        if (sigma <= 0d0) return
        ret = exp(-0.5d0*((x - mu)/sigma)**2)/(sqrt(2d0*pi)*sigma)
    end function

end module

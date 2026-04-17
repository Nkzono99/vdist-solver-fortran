module m_photoelectron_raycast
    !! Raycast-based photoelectron probability.
    !!
    !! When a backtraced particle collides with an internal boundary, this
    !! probability function casts a ray from the collision point toward the
    !! sun (the direction opposite to the photoelectron drift vector derived
    !! from `vdri`, `vdthz`, `vdthxy`). If the ray is blocked by another
    !! internal boundary the probability is zero (the emission surface is
    !! in shadow); otherwise it equals the Maxwellian PDF of the photo-
    !! electron velocity distribution evaluated at the incoming particle
    !! velocity. The caller is expected to hand over a `blocking_boundaries`
    !! list containing only geometry that can occlude the sun — typically
    !! the internal surfaces/objects, not the outer domain planes.

    use m_probabilities, only: t_Probability
    use finbound, only: t_BoundaryList, t_Ray, new_Ray, t_HitRecord

    implicit none

    private
    public t_PhotoelectronRaycastProbability
    public new_PhotoelectronRaycastProbability

    double precision, parameter :: pi = acos(-1.0d0)

    double precision, parameter :: DEFAULT_RAY_ORIGIN_OFFSET = 1d-6
        !! Small offset along the ray direction used to move the ray origin
        !! just off the surface it starts from, so the ray does not
        !! immediately self-intersect the host boundary.

    type, extends(t_Probability) :: t_PhotoelectronRaycastProbability
        double precision :: locs(3)
            !! Mean velocity of the photoelectron Maxwellian.
        double precision :: scales(3)
            !! Thermal velocity (sigma) of the photoelectron Maxwellian.
        double precision :: sun_direction(3)
            !! Unit vector pointing from the emission point toward the sun.
            !! This is the direction the occlusion ray is cast along.
        type(t_BoundaryList) :: blocking_boundaries
            !! Boundary list used to test occlusion — typically populated
            !! with only the internal surfaces/objects (no outer planes).
        double precision :: ray_origin_offset = DEFAULT_RAY_ORIGIN_OFFSET
            !! Offset along `sun_direction` applied to the ray origin to
            !! avoid self-intersection at the host boundary.
    contains
        procedure :: at => photoelectronRaycast_at
    end type

contains

    function new_PhotoelectronRaycastProbability(locs, scales, sun_direction, &
                                                 blocking_boundaries, &
                                                 coefficient, &
                                                 ray_origin_offset) result(obj)
        double precision, intent(in) :: locs(3)
        double precision, intent(in) :: scales(3)
        double precision, intent(in) :: sun_direction(3)
        type(t_BoundaryList), intent(in) :: blocking_boundaries
        double precision, intent(in), optional :: coefficient
        double precision, intent(in), optional :: ray_origin_offset
        type(t_PhotoelectronRaycastProbability) :: obj

        double precision :: mag

        obj%locs = locs
        obj%scales = scales

        mag = norm2(sun_direction)
        if (mag > 0d0) then
            obj%sun_direction = sun_direction/mag
        else
            ! Degenerate (no drift) -> default to +Z (overhead sun).
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
        !! Probability density for a photoelectron emission event.
        !!
        !! Physical model:
        !!   1. Photoelectrons leave the surface with velocity whose component
        !!      along the outward normal (approximated by `sun_direction`) is
        !!      strictly positive. Particles with v . sun_direction <= 0 could
        !!      not have been emitted outward, so their density is 0.
        !!   2. The remaining density is a shifted 3D Maxwellian multiplied by
        !!      a half-space normalization factor 2 (exact for mu = 0 along the
        !!      normal; a reasonable approximation for drifting distributions
        !!      when |drift| is not dominated by the thermal spread).
        !!   3. If a blocking boundary intercepts the ray toward the sun the
        !!      surface element is in shadow and no photoelectron could have
        !!      been produced there, so the density is again 0.
        class(t_PhotoelectronRaycastProbability), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        double precision :: ret

        type(t_Ray) :: ray
        type(t_HitRecord) :: hit
        double precision :: origin(3)
        double precision :: v_parallel
        integer :: i

        ! Half-space check: only keep particles moving outward (along +sun_dir).
        v_parallel = dot_product(velocity, self%sun_direction)
        if (v_parallel <= 0d0) then
            ret = 0d0
            return
        end if

        ret = 2d0*self%coefficient
        do i = 1, 3
            ret = ret*maxwell_pdf(velocity(i), self%locs(i), self%scales(i))
        end do

        if (ret == 0d0) return

        origin = position + self%ray_origin_offset*self%sun_direction
        ray = new_Ray(origin, self%sun_direction)
        hit = self%blocking_boundaries%hit(ray)

        if (hit%is_hit .and. hit%t > 0d0) then
            ret = 0d0
        end if
    end function

    pure function maxwell_pdf(x, mu, sigma) result(ret)
        double precision, intent(in) :: x
        double precision, intent(in) :: mu
        double precision, intent(in) :: sigma
        double precision :: ret

        ret = 1d0/sqrt(2d0*pi*sigma*sigma)*exp(-(x - mu)*(x - mu)/(2d0*sigma*sigma))
    end function

end module

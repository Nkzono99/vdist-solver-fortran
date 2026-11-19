!> Shared primitives for m_objects and m_surfaces: the subdomain
!> overlap padding parameter and the "allocate → try to add, else
!> discard" boundary helper used dozens of times in both modules.
module m_boundary_helpers
  use finbound, only: t_Boundary, t_BoundaryList
  use m_collision_accelerator, only: t_collision_accelerator
  use m_kinds, only: dp, ip
  implicit none
  private
  public :: try_add_boundary
  public :: extent

  !> Extra padding applied to is_overlap(...)'s subdomain extents. Two
  !> grid cells on each face in every direction — matches the original
  !> constant used by both m_objects and m_surfaces.
  real(kind=dp), parameter :: extent(2, 3) = &
                              reshape([[1.0d0, 1.0d0], [1.0d0, 1.0d0], [1.0d0, 1.0d0]]*2, [2, 3])

contains

  !> Stamp `tag` onto `pbound%material` and add the boundary to the
  !> list iff it overlaps the subdomain (with the shared `extent`
  !> padding); otherwise deallocate the polymorphic pointer.
  !>
  !> Callers retain a specific-type local pointer (e.g. t_RectangleXYZ)
  !> that aliases the same target as `pbound`. After this routine the
  !> target is either owned by `boundaries` or freed, and the caller's
  !> specific-type pointer must not be dereferenced again.
  subroutine try_add_boundary(boundaries, pbound, tag, sdoms, accelerator, aabb)
    type(t_BoundaryList), intent(inout) :: boundaries
    class(t_Boundary), pointer, intent(inout) :: pbound
    integer(kind=ip), intent(in) :: tag
    real(kind=dp), intent(in) :: sdoms(2, 3)
    type(t_collision_accelerator), intent(inout), optional :: accelerator
    real(kind=dp), intent(in), optional :: aabb(2, 3)

    pbound%material%tag = tag
    if (pbound%is_overlap(sdoms, extent=extent)) then
      call boundaries%add_boundary(pbound)
      if (present(accelerator)) then
        if (present(aabb)) then
          call accelerator%add_aabb(aabb)
        else
          call accelerator%add_unbounded()
        end if
      end if
    else
      deallocate (pbound)
    end if
  end subroutine try_add_boundary

end module m_boundary_helpers

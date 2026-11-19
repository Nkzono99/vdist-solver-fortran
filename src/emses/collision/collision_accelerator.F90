!> Lightweight broad-phase acceleration for finbound collision checks.
module m_collision_accelerator
  use m_kinds, only: dp, ip, lp
  implicit none
  private

  real(kind=dp), parameter :: DEFAULT_AABB_PADDING = 1.0d-9
  real(kind=dp), parameter :: UNBOUNDED_LIMIT = 1.0d100
  integer(kind=ip), parameter :: DEFAULT_MAX_NBOUNDARIES = 30_ip
  integer(kind=ip), parameter :: DEFAULT_CELL_DOMAIN_PADDING = 2_ip
  integer(kind=lp), parameter :: MAX_CELL_INDEX_CELLS = 50000000_lp

  type, public :: t_collision_accelerator
    integer(kind=ip) :: nboundaries = 0_ip
    integer(kind=ip) :: max_nboundaries = 0_ip
    real(kind=dp), allocatable :: lower(:, :)
    real(kind=dp), allocatable :: upper(:, :)
    logical :: has_aggregate_bounds = .false.
    real(kind=dp) :: aggregate_lower(3) = UNBOUNDED_LIMIT
    real(kind=dp) :: aggregate_upper(3) = -UNBOUNDED_LIMIT
    logical :: cell_index_enabled = .false.
    integer(kind=ip) :: cell_lower(3) = 0_ip
    integer(kind=ip) :: cell_upper(3) = -1_ip
    integer(kind=ip) :: cell_counts(3) = 0_ip
    logical, allocatable :: active_cell(:)
  contains
    procedure :: initialize => collision_accelerator_initialize
    procedure :: add_aabb => collision_accelerator_add_aabb
    procedure :: add_unbounded => collision_accelerator_add_unbounded
    procedure :: build_cell_index => collision_accelerator_build_cell_index
    procedure :: segment_bounds => collision_accelerator_segment_bounds
    procedure :: padded_segment_bounds => collision_accelerator_padded_segment_bounds
    procedure :: segment_may_overlap_any => collision_accelerator_segment_may_overlap_any
    procedure :: bounds_may_overlap_any => collision_accelerator_bounds_may_overlap_any
    procedure :: bounds_may_touch_active_cell => collision_accelerator_bounds_may_touch_active_cell
    procedure :: bounds_may_overlap => collision_accelerator_bounds_may_overlap
    procedure :: segment_may_touch_active_cell => collision_accelerator_segment_may_touch_active_cell
    procedure :: segment_may_overlap => collision_accelerator_segment_may_overlap
    procedure :: destroy => collision_accelerator_destroy
  end type t_collision_accelerator

  public :: aabb_from_axis_plane
  public :: aabb_from_bounded_axis_plane
  public :: aabb_from_bounds
  public :: aabb_from_centered_profile
  public :: aabb_from_circle_xyz
  public :: aabb_from_cylinder_xyz
  public :: aabb_from_quad
  public :: aabb_from_rectangle_xyz
  public :: aabb_from_sphere

contains

  subroutine collision_accelerator_initialize(self, initial_capacity)
    class(t_collision_accelerator), intent(inout) :: self
    integer(kind=ip), intent(in), optional :: initial_capacity

    integer(kind=ip) :: capacity

    call self%destroy()

    capacity = DEFAULT_MAX_NBOUNDARIES
    if (present(initial_capacity)) capacity = max(1_ip, initial_capacity)

    allocate (self%lower(3, capacity))
    allocate (self%upper(3, capacity))
    self%nboundaries = 0_ip
    self%max_nboundaries = capacity
  end subroutine collision_accelerator_initialize

  subroutine collision_accelerator_add_aabb(self, aabb)
    class(t_collision_accelerator), intent(inout) :: self
    real(kind=dp), intent(in) :: aabb(2, 3)

    if (.not. allocated(self%lower)) call self%initialize()
    if (self%nboundaries >= self%max_nboundaries) call collision_accelerator_extend(self)
    if (self%cell_index_enabled) call collision_accelerator_destroy_cell_index(self)

    self%nboundaries = self%nboundaries + 1_ip
    self%lower(:, self%nboundaries) = min(aabb(1, :), aabb(2, :)) - DEFAULT_AABB_PADDING
    self%upper(:, self%nboundaries) = max(aabb(1, :), aabb(2, :)) + DEFAULT_AABB_PADDING
    call collision_accelerator_update_aggregate_bounds(self, self%nboundaries)
  end subroutine collision_accelerator_add_aabb

  subroutine collision_accelerator_add_unbounded(self)
    class(t_collision_accelerator), intent(inout) :: self

    call self%add_aabb(aabb_from_bounds([-UNBOUNDED_LIMIT, -UNBOUNDED_LIMIT, -UNBOUNDED_LIMIT], &
                                        [UNBOUNDED_LIMIT, UNBOUNDED_LIMIT, UNBOUNDED_LIMIT]))
  end subroutine collision_accelerator_add_unbounded

  subroutine collision_accelerator_build_cell_index(self, isdoms, padding)
    class(t_collision_accelerator), intent(inout) :: self
    integer(kind=ip), intent(in) :: isdoms(2, 3)
    integer(kind=ip), intent(in), optional :: padding

    integer(kind=ip) :: pad
    integer(kind=ip) :: iboundary, ix, iy, iz
    integer(kind=ip) :: lower_cell(3), upper_cell(3)
    integer(kind=lp) :: total_cells

    call collision_accelerator_destroy_cell_index(self)
    if (.not. allocated(self%lower)) return

    pad = DEFAULT_CELL_DOMAIN_PADDING
    if (present(padding)) pad = max(0_ip, padding)

    self%cell_lower(:) = isdoms(1, :) - pad
    self%cell_upper(:) = isdoms(2, :) + pad
    self%cell_counts(:) = self%cell_upper(:) - self%cell_lower(:) + 1_ip
    if (any(self%cell_counts(:) <= 0_ip)) return

    total_cells = int(self%cell_counts(1), lp)*int(self%cell_counts(2), lp)*int(self%cell_counts(3), lp)
    if (total_cells > MAX_CELL_INDEX_CELLS) return

    allocate (self%active_cell(int(total_cells, ip)))
    self%active_cell(:) = .false.

    do iboundary = 1, self%nboundaries
      call boundary_cell_range(self, iboundary, lower_cell, upper_cell)
      if (any(upper_cell(:) < lower_cell(:))) cycle

      do iz = lower_cell(3), upper_cell(3)
        do iy = lower_cell(2), upper_cell(2)
          do ix = lower_cell(1), upper_cell(1)
            self%active_cell(cell_flat_index(self, ix, iy, iz)) = .true.
          end do
        end do
      end do
    end do

    self%cell_index_enabled = .true.
  end subroutine collision_accelerator_build_cell_index

  pure subroutine collision_accelerator_segment_bounds(self, p1, p2, segment_lower, segment_upper)
    class(t_collision_accelerator), intent(in) :: self
    real(kind=dp), intent(in) :: p1(3), p2(3)
    real(kind=dp), intent(out) :: segment_lower(3), segment_upper(3)

    segment_lower(:) = min(p1(:), p2(:))
    segment_upper(:) = max(p1(:), p2(:))
  end subroutine collision_accelerator_segment_bounds

  pure subroutine collision_accelerator_padded_segment_bounds(self, p1, p2, segment_lower, segment_upper)
    class(t_collision_accelerator), intent(in) :: self
    real(kind=dp), intent(in) :: p1(3), p2(3)
    real(kind=dp), intent(out) :: segment_lower(3), segment_upper(3)

    segment_lower(:) = min(p1(:), p2(:)) - DEFAULT_AABB_PADDING
    segment_upper(:) = max(p1(:), p2(:)) + DEFAULT_AABB_PADDING
  end subroutine collision_accelerator_padded_segment_bounds

  pure function collision_accelerator_segment_may_overlap_any(self, p1, p2) result(may_overlap)
    class(t_collision_accelerator), intent(in) :: self
    real(kind=dp), intent(in) :: p1(3), p2(3)
    logical :: may_overlap

    real(kind=dp) :: segment_lower(3), segment_upper(3)

    call self%padded_segment_bounds(p1, p2, segment_lower, segment_upper)
    may_overlap = self%bounds_may_overlap_any(segment_lower, segment_upper)
  end function collision_accelerator_segment_may_overlap_any

  pure function collision_accelerator_bounds_may_overlap_any(self, segment_lower, segment_upper) result(may_overlap)
    class(t_collision_accelerator), intent(in) :: self
    real(kind=dp), intent(in) :: segment_lower(3), segment_upper(3)
    logical :: may_overlap

    if (self%nboundaries <= 0_ip .or. .not. self%has_aggregate_bounds) then
      may_overlap = .false.
      return
    end if

    may_overlap = .not. any(segment_upper(:) < self%aggregate_lower(:) &
                            .or. self%aggregate_upper(:) < segment_lower(:))
  end function collision_accelerator_bounds_may_overlap_any

  pure function collision_accelerator_bounds_may_touch_active_cell(self, segment_lower, segment_upper) result(may_touch)
    class(t_collision_accelerator), intent(in) :: self
    real(kind=dp), intent(in) :: segment_lower(3), segment_upper(3)
    logical :: may_touch

    integer(kind=ip) :: iaxis, ix, iy, iz
    integer(kind=ip) :: lower_cell(3), upper_cell(3)

    if (.not. self%cell_index_enabled .or. .not. allocated(self%active_cell)) then
      may_touch = .true.
      return
    end if

    do iaxis = 1, 3
      if (segment_lower(iaxis) < real(self%cell_lower(iaxis), dp) &
          .or. segment_upper(iaxis) > real(self%cell_upper(iaxis) + 1_ip, dp)) then
        may_touch = .true.
        return
      end if

      lower_cell(iaxis) = clamped_floor_to_cell(segment_lower(iaxis), self%cell_lower(iaxis), self%cell_upper(iaxis))
      upper_cell(iaxis) = clamped_floor_to_cell(segment_upper(iaxis), self%cell_lower(iaxis), self%cell_upper(iaxis))
    end do

    do iz = lower_cell(3), upper_cell(3)
      do iy = lower_cell(2), upper_cell(2)
        do ix = lower_cell(1), upper_cell(1)
          if (self%active_cell(cell_flat_index(self, ix, iy, iz))) then
            may_touch = .true.
            return
          end if
        end do
      end do
    end do

    may_touch = .false.
  end function collision_accelerator_bounds_may_touch_active_cell

  pure function collision_accelerator_segment_may_touch_active_cell(self, p1, p2) result(may_touch)
    class(t_collision_accelerator), intent(in) :: self
    real(kind=dp), intent(in) :: p1(3), p2(3)
    logical :: may_touch

    real(kind=dp) :: segment_lower(3), segment_upper(3)

    call self%segment_bounds(p1, p2, segment_lower, segment_upper)
    may_touch = self%bounds_may_touch_active_cell(segment_lower, segment_upper)
  end function collision_accelerator_segment_may_touch_active_cell

  pure function collision_accelerator_bounds_may_overlap(self, iboundary, segment_lower, segment_upper) result(may_overlap)
    class(t_collision_accelerator), intent(in) :: self
    integer(kind=ip), intent(in) :: iboundary
    real(kind=dp), intent(in) :: segment_lower(3), segment_upper(3)
    logical :: may_overlap

    if (.not. allocated(self%lower)) then
      may_overlap = .true.
      return
    end if

    if (iboundary < 1_ip .or. iboundary > self%nboundaries) then
      may_overlap = .true.
      return
    end if

    may_overlap = .not. any(segment_upper(:) < self%lower(:, iboundary) &
                            .or. self%upper(:, iboundary) < segment_lower(:))
  end function collision_accelerator_bounds_may_overlap

  pure function collision_accelerator_segment_may_overlap(self, iboundary, p1, p2) result(may_overlap)
    class(t_collision_accelerator), intent(in) :: self
    integer(kind=ip), intent(in) :: iboundary
    real(kind=dp), intent(in) :: p1(3), p2(3)
    logical :: may_overlap

    real(kind=dp) :: segment_lower(3), segment_upper(3)

    call self%padded_segment_bounds(p1, p2, segment_lower, segment_upper)
    may_overlap = self%bounds_may_overlap(iboundary, segment_lower, segment_upper)
  end function collision_accelerator_segment_may_overlap

  subroutine collision_accelerator_destroy(self)
    class(t_collision_accelerator), intent(inout) :: self

    if (allocated(self%lower)) deallocate (self%lower)
    if (allocated(self%upper)) deallocate (self%upper)
    call collision_accelerator_destroy_cell_index(self)
    self%nboundaries = 0_ip
    self%max_nboundaries = 0_ip
    self%has_aggregate_bounds = .false.
    self%aggregate_lower(:) = UNBOUNDED_LIMIT
    self%aggregate_upper(:) = -UNBOUNDED_LIMIT
  end subroutine collision_accelerator_destroy

  subroutine collision_accelerator_destroy_cell_index(self)
    type(t_collision_accelerator), intent(inout) :: self

    if (allocated(self%active_cell)) deallocate (self%active_cell)
    self%cell_index_enabled = .false.
    self%cell_lower(:) = 0_ip
    self%cell_upper(:) = -1_ip
    self%cell_counts(:) = 0_ip
  end subroutine collision_accelerator_destroy_cell_index

  subroutine collision_accelerator_extend(self)
    type(t_collision_accelerator), intent(inout) :: self

    integer(kind=ip) :: new_capacity
    real(kind=dp), allocatable :: new_lower(:, :)
    real(kind=dp), allocatable :: new_upper(:, :)

    new_capacity = max(self%max_nboundaries + 1_ip, int(real(self%max_nboundaries, dp)*1.5d0, ip))
    allocate (new_lower(3, new_capacity))
    allocate (new_upper(3, new_capacity))
    if (self%nboundaries > 0_ip) then
      new_lower(:, 1:self%nboundaries) = self%lower(:, 1:self%nboundaries)
      new_upper(:, 1:self%nboundaries) = self%upper(:, 1:self%nboundaries)
    end if

    call move_alloc(new_lower, self%lower)
    call move_alloc(new_upper, self%upper)
    self%max_nboundaries = new_capacity
  end subroutine collision_accelerator_extend

  subroutine collision_accelerator_update_aggregate_bounds(self, iboundary)
    type(t_collision_accelerator), intent(inout) :: self
    integer(kind=ip), intent(in) :: iboundary

    if (.not. self%has_aggregate_bounds) then
      self%aggregate_lower(:) = self%lower(:, iboundary)
      self%aggregate_upper(:) = self%upper(:, iboundary)
      self%has_aggregate_bounds = .true.
    else
      self%aggregate_lower(:) = min(self%aggregate_lower(:), self%lower(:, iboundary))
      self%aggregate_upper(:) = max(self%aggregate_upper(:), self%upper(:, iboundary))
    end if
  end subroutine collision_accelerator_update_aggregate_bounds

  subroutine boundary_cell_range(self, iboundary, lower_cell, upper_cell)
    type(t_collision_accelerator), intent(in) :: self
    integer(kind=ip), intent(in) :: iboundary
    integer(kind=ip), intent(out) :: lower_cell(3), upper_cell(3)

    integer(kind=ip) :: iaxis

    lower_cell(:) = 1_ip
    upper_cell(:) = 0_ip

    do iaxis = 1, 3
      if (self%upper(iaxis, iboundary) < real(self%cell_lower(iaxis), dp) &
          .or. self%lower(iaxis, iboundary) > real(self%cell_upper(iaxis) + 1_ip, dp)) return

      lower_cell(iaxis) = clamped_floor_to_cell(self%lower(iaxis, iboundary), &
                                                self%cell_lower(iaxis), self%cell_upper(iaxis))
      upper_cell(iaxis) = clamped_floor_to_cell(self%upper(iaxis, iboundary), &
                                                self%cell_lower(iaxis), self%cell_upper(iaxis))
    end do
  end subroutine boundary_cell_range

  pure integer(kind=ip) function clamped_floor_to_cell(value, lower_bound, upper_bound) result(cell)
    real(kind=dp), intent(in) :: value
    integer(kind=ip), intent(in) :: lower_bound, upper_bound

    if (value <= real(lower_bound, dp)) then
      cell = lower_bound
    else if (value >= real(upper_bound + 1_ip, dp)) then
      cell = upper_bound
    else
      cell = floor(value, kind=ip)
      cell = max(lower_bound, min(upper_bound, cell))
    end if
  end function clamped_floor_to_cell

  pure integer(kind=ip) function cell_flat_index(self, ix, iy, iz) result(index)
    type(t_collision_accelerator), intent(in) :: self
    integer(kind=ip), intent(in) :: ix, iy, iz

    index = 1_ip + (ix - self%cell_lower(1)) &
            + self%cell_counts(1)*(iy - self%cell_lower(2)) &
            + self%cell_counts(1)*self%cell_counts(2)*(iz - self%cell_lower(3))
  end function cell_flat_index

  pure function aabb_from_bounds(lower, upper) result(aabb)
    real(kind=dp), intent(in) :: lower(3), upper(3)
    real(kind=dp) :: aabb(2, 3)

    aabb(1, :) = min(lower(:), upper(:))
    aabb(2, :) = max(lower(:), upper(:))
  end function aabb_from_bounds

  pure function aabb_from_axis_plane(axis, pos) result(aabb)
    integer(kind=ip), intent(in) :: axis
    real(kind=dp), intent(in) :: pos
    real(kind=dp) :: aabb(2, 3)

    aabb(1, :) = -UNBOUNDED_LIMIT
    aabb(2, :) = UNBOUNDED_LIMIT
    aabb(1, axis) = pos
    aabb(2, axis) = pos
  end function aabb_from_axis_plane

  pure function aabb_from_bounded_axis_plane(axis, origin, lower_axis1, upper_axis1, lower_axis2, upper_axis2) result(aabb)
    integer(kind=ip), intent(in) :: axis
    real(kind=dp), intent(in) :: origin(3)
    real(kind=dp), intent(in) :: lower_axis1, upper_axis1, lower_axis2, upper_axis2
    real(kind=dp) :: aabb(2, 3)

    integer(kind=ip) :: axis1, axis2

    axis1 = mod(axis, 3) + 1
    axis2 = mod(axis + 1, 3) + 1

    aabb(1, :) = -UNBOUNDED_LIMIT
    aabb(2, :) = UNBOUNDED_LIMIT
    aabb(1, axis) = origin(axis)
    aabb(2, axis) = origin(axis)
    aabb(1, axis1) = min(lower_axis1, upper_axis1)
    aabb(2, axis1) = max(lower_axis1, upper_axis1)
    aabb(1, axis2) = min(lower_axis2, upper_axis2)
    aabb(2, axis2) = max(lower_axis2, upper_axis2)
  end function aabb_from_bounded_axis_plane

  pure function aabb_from_rectangle_xyz(axis, origin, w1, w2) result(aabb)
    integer(kind=ip), intent(in) :: axis
    real(kind=dp), intent(in) :: origin(3), w1, w2
    real(kind=dp) :: aabb(2, 3)

    integer(kind=ip) :: axis1, axis2

    axis1 = mod(axis, 3) + 1
    axis2 = mod(axis + 1, 3) + 1

    aabb(1, :) = origin(:)
    aabb(2, :) = origin(:)
    aabb(2, axis1) = origin(axis1) + w1
    aabb(2, axis2) = origin(axis2) + w2
    aabb = aabb_from_bounds(aabb(1, :), aabb(2, :))
  end function aabb_from_rectangle_xyz

  pure function aabb_from_circle_xyz(axis, origin, radius) result(aabb)
    integer(kind=ip), intent(in) :: axis
    real(kind=dp), intent(in) :: origin(3), radius
    real(kind=dp) :: aabb(2, 3)

    integer(kind=ip) :: axis1, axis2

    axis1 = mod(axis, 3) + 1
    axis2 = mod(axis + 1, 3) + 1

    aabb(1, :) = origin(:)
    aabb(2, :) = origin(:)
    aabb(1, axis1) = origin(axis1) - radius
    aabb(2, axis1) = origin(axis1) + radius
    aabb(1, axis2) = origin(axis2) - radius
    aabb(2, axis2) = origin(axis2) + radius
  end function aabb_from_circle_xyz

  pure function aabb_from_cylinder_xyz(axis, origin, radius, height) result(aabb)
    integer(kind=ip), intent(in) :: axis
    real(kind=dp), intent(in) :: origin(3), radius, height
    real(kind=dp) :: aabb(2, 3)

    integer(kind=ip) :: axis1, axis2

    axis1 = mod(axis, 3) + 1
    axis2 = mod(axis + 1, 3) + 1

    aabb(1, :) = origin(:)
    aabb(2, :) = origin(:)
    aabb(1, axis1) = origin(axis1) - radius
    aabb(2, axis1) = origin(axis1) + radius
    aabb(1, axis2) = origin(axis2) - radius
    aabb(2, axis2) = origin(axis2) + radius
    aabb(2, axis) = origin(axis) + height
    aabb = aabb_from_bounds(aabb(1, :), aabb(2, :))
  end function aabb_from_cylinder_xyz

  pure function aabb_from_centered_profile(axis, origin, radius, height) result(aabb)
    integer(kind=ip), intent(in) :: axis
    real(kind=dp), intent(in) :: origin(3), radius, height
    real(kind=dp) :: aabb(2, 3)

    integer(kind=ip) :: axis1, axis2

    axis1 = mod(axis, 3) + 1
    axis2 = mod(axis + 1, 3) + 1

    aabb(1, :) = origin(:)
    aabb(2, :) = origin(:)
    aabb(1, axis) = origin(axis) - 0.5d0*height
    aabb(2, axis) = origin(axis) + 0.5d0*height
    aabb(1, axis1) = origin(axis1) - radius
    aabb(2, axis1) = origin(axis1) + radius
    aabb(1, axis2) = origin(axis2) - radius
    aabb(2, axis2) = origin(axis2) + radius
  end function aabb_from_centered_profile

  pure function aabb_from_sphere(origin, radius) result(aabb)
    real(kind=dp), intent(in) :: origin(3), radius
    real(kind=dp) :: aabb(2, 3)

    aabb(1, :) = origin(:) - radius
    aabb(2, :) = origin(:) + radius
  end function aabb_from_sphere

  pure function aabb_from_quad(corners) result(aabb)
    real(kind=dp), intent(in) :: corners(3, 4)
    real(kind=dp) :: aabb(2, 3)

    integer(kind=ip) :: iaxis

    do iaxis = 1, 3
      aabb(1, iaxis) = minval(corners(iaxis, :))
      aabb(2, iaxis) = maxval(corners(iaxis, :))
    end do
  end function aabb_from_quad

end module m_collision_accelerator

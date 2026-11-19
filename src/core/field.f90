module m_field
    !! This module defines vector fields and their grid representations.
    implicit none

    private
    public t_VectorField
    public t_VectorFieldGrid
    public new_VectorFieldGrid

    type, abstract :: t_VectorField
        !! Abstract type for a vector field
        integer :: n_elements
            !! Number of elements in the vector field
    contains
        procedure(vectorField_at), deferred :: at
            !! Get the vector field value at a given position.
    end type

    interface
        function vectorField_at(self, position) result(ret)
            !! Get the vector field value at a given position.
            import t_VectorField
            class(t_VectorField), intent(in) :: self
                !! Instance of the vector field
            double precision, intent(in) :: position(3)
                !! Position in 3D space
            double precision :: ret(self%n_elements)
                !! Vector field value at the given position
        end function
    end interface

    type, extends(t_VectorField) :: t_VectorFieldGrid
        !! A vector field defined on a grid
        integer :: nx
            !! Number of grid points in the x direction
        integer :: ny
            !! Number of grid points in the y direction
        integer :: nz
            !! Number of grid points in the z direction
        double precision, allocatable :: values(:, :, :, :)
            !! Values of the vector field at the grid points
        double precision, allocatable :: offsets(:, :)
            !! Per-element interpolation offsets in x, y, z.
        logical :: nodal = .true.
            !! All offsets vanish, so every element shares one stencil.
    contains
        procedure :: at => vectorFieldGrid_at
            !! Get the vector field value at a given position
    end type

contains

    function new_VectorFieldGrid(n_elements, nx, ny, nz, values, offsets) result(obj)
        !! Create a new vector field grid object with specified properties.

        integer, intent(in) :: n_elements
            !! Number of elements in the vector field
        integer, intent(in) :: nx
            !! Number of grid cells in the x direction
        integer, intent(in) :: ny
            !! Number of grid cells in the y direction
        integer, intent(in) :: nz
            !! Number of grid cells in the z direction
        double precision, intent(in) :: values(n_elements, 0:nx, 0:ny, 0:nz)
            !! Values of the vector field at the grid points
        double precision, intent(in), optional :: offsets(3, n_elements)
            !! Per-element interpolation offsets in x, y, z.
        type(t_VectorFieldGrid) :: obj
            !! A new vector field grid object with specified properties

        obj%n_elements = n_elements
        obj%nx = nx
        obj%ny = ny
        obj%nz = nz
        allocate (obj%values(n_elements, 0:nx, 0:ny, 0:nz))
        obj%values(:, :, :, :) = values(:, :, :, :)
        allocate (obj%offsets(3, n_elements))
        if (present(offsets)) then
            obj%offsets(:, :) = offsets(:, :)
        else
            obj%offsets(:, :) = 0d0
        end if
        obj%nodal = all(abs(obj%offsets) <= epsilon(1d0))
    end function

    function vectorFieldGrid_at(self, position) result(ret)
        !! Get the vector field value at a given position using linear interpolation.
        !!
        !! Positions outside the grid are linearly extrapolated from the edge
        !! cell; an axis with a single node is constant. The nodal cell is
        !! located once; an element locates again only along its offset axes
        !! (for example one axis per component of a staggered field).

        class(t_VectorFieldGrid), intent(in) :: self
            !! Instance of the vector field grid
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision :: ret(self%n_elements)
            !! Vector field value at the given position

        integer :: nodes(3), node_cell(3), node_next(3), cell(3), next(3), i, axis
        double precision :: node_fraction(3), fraction(3)

        nodes = [self%nx, self%ny, self%nz]
        do axis = 1, 3
            call locate_axis(nodes(axis), position(axis), node_cell(axis), node_next(axis), node_fraction(axis))
        end do

        if (self%nodal) then
            do i = 1, self%n_elements
                ret(i) = interpolate_element(self, i, node_cell, node_next, node_fraction)
            end do
            return
        end if

        do i = 1, self%n_elements
            cell = node_cell
            next = node_next
            fraction = node_fraction
            do axis = 1, 3
                if (self%offsets(axis, i) /= 0d0) then
                    call locate_axis(nodes(axis), position(axis) - self%offsets(axis, i), &
                                     cell(axis), next(axis), fraction(axis))
                end if
            end do
            ret(i) = interpolate_element(self, i, cell, next, fraction)
        end do
    end function

    pure subroutine locate_axis(nodes, position, cell, next, fraction)
        !! Lower and upper node index and fractional position along one axis.
        !! The lower node is clamped to [0, nodes-1] so outer points extrapolate;
        !! a single-node axis (nodes = 0) is constant.
        integer, intent(in) :: nodes
        double precision, intent(in) :: position
        integer, intent(out) :: cell, next
        double precision, intent(out) :: fraction

        if (nodes == 0) then
            cell = 0
            next = 0
            fraction = 0d0
        else
            cell = max(0, min(int(position), nodes - 1))
            next = cell + 1
            fraction = position - cell
        end if
    end subroutine

    pure function interpolate_element(self, element_index, cell, next, fraction) result(ret)
        !! Trilinear interpolation of one element on a located cell.
        class(t_VectorFieldGrid), intent(in) :: self
        integer, intent(in) :: element_index
        integer, intent(in) :: cell(3), next(3)
        double precision, intent(in) :: fraction(3)
        double precision :: ret

        double precision :: rp(3), rp1(3)
        double precision :: u00, u01, u10, u11, u0, u1

        rp = fraction
        rp1 = 1 - rp

        u00 = rp(1)*self%values(element_index, next(1), cell(2), cell(3)) &
              + rp1(1)*self%values(element_index, cell(1), cell(2), cell(3))
        u01 = rp(1)*self%values(element_index, next(1), next(2), cell(3)) &
              + rp1(1)*self%values(element_index, cell(1), next(2), cell(3))
        u10 = rp(1)*self%values(element_index, next(1), cell(2), next(3)) &
              + rp1(1)*self%values(element_index, cell(1), cell(2), next(3))
        u11 = rp(1)*self%values(element_index, next(1), next(2), next(3)) &
              + rp1(1)*self%values(element_index, cell(1), next(2), next(3))

        u0 = rp(2)*u01 + rp1(2)*u00
        u1 = rp(2)*u11 + rp1(2)*u10

        ret = rp(3)*u1 + rp1(3)*u0
    end function

end module

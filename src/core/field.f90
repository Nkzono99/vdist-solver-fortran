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
    end function

    function vectorFieldGrid_at(self, position) result(ret)
        !! Get the vector field value at a given position using linear interpolation.

        class(t_VectorFieldGrid), intent(in) :: self
            !! Instance of the vector field grid
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision :: ret(self%n_elements)
            !! Vector field value at the given position

        double precision :: p(3)
        integer :: ip(3), ip1(3)
        double precision :: rp(3), rp1(3)
        double precision :: u00(self%n_elements), u01(self%n_elements)
        double precision :: u10(self%n_elements), u11(self%n_elements)
        double precision :: u0(self%n_elements), u1(self%n_elements)

        integer :: i

        if (any(abs(self%offsets) > epsilon(1d0))) then
            do i = 1, self%n_elements
                ret(i) = vectorFieldGrid_component_at(self, i, position - self%offsets(:, i))
            end do
            return
        end if

        ! Linear interpolation
        p(:) = position(:)
        do i = 1, 3
            ! ip(i) を [0 … n_dim–1] の範囲に制限する
            ! もし position(i) < 0    → 0 に丸め込む
            !    position(i) ≥ n_dim → n_dim–1 に丸め込む
            ip(i) = int(p(i))
            if (ip(i) < 0)         ip(i) = 0
            if (ip(i) > self%nx-1 .and. i == 1) ip(i) = self%nx - 1
            if (ip(i) > self%ny-1 .and. i == 2) ip(i) = self%ny - 1
            if (ip(i) > self%nz-1 .and. i == 3) ip(i) = self%nz - 1

            ip1(i) = ip(i) + 1
        end do

        rp(:) = p(:) - ip(:)
        rp1(:) = 1 - rp(:)

        u00 = rp(1)*self%values(:, ip1(1), ip(2), ip(3)) &
              + rp1(1)*self%values(:, ip(1), ip(2), ip(3))
        u01 = rp(1)*self%values(:, ip1(1), ip1(2), ip(3)) &
              + rp1(1)*self%values(:, ip(1), ip1(2), ip(3))
        u10 = rp(1)*self%values(:, ip1(1), ip(2), ip1(3)) &
              + rp1(1)*self%values(:, ip(1), ip(2), ip1(3))
        u11 = rp(1)*self%values(:, ip1(1), ip1(2), ip1(3)) &
              + rp1(1)*self%values(:, ip(1), ip1(2), ip1(3))

        u0 = rp(2)*u01 + rp1(2)*u00
        u1 = rp(2)*u11 + rp1(2)*u10

        ret(:) = rp(3)*u1 + rp1(3)*u0
    end function

    function vectorFieldGrid_component_at(self, element_index, position) result(ret)
        !! Get one vector-field element at a shifted position.

        class(t_VectorFieldGrid), intent(in) :: self
            !! Instance of the vector field grid
        integer, intent(in) :: element_index
            !! Element index to interpolate
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision :: ret
            !! Vector field element value at the given position

        double precision :: p(3)
        integer :: ip(3), ip1(3)
        double precision :: rp(3), rp1(3)
        double precision :: u00, u01, u10, u11
        double precision :: u0, u1

        integer :: i

        p(:) = position(:)
        do i = 1, 3
            ip(i) = int(p(i))
            if (ip(i) < 0) ip(i) = 0
            if (ip(i) > self%nx - 1 .and. i == 1) ip(i) = self%nx - 1
            if (ip(i) > self%ny - 1 .and. i == 2) ip(i) = self%ny - 1
            if (ip(i) > self%nz - 1 .and. i == 3) ip(i) = self%nz - 1

            ip1(i) = ip(i) + 1
        end do

        rp(:) = p(:) - ip(:)
        rp1(:) = 1 - rp(:)

        u00 = rp(1)*self%values(element_index, ip1(1), ip(2), ip(3)) &
              + rp1(1)*self%values(element_index, ip(1), ip(2), ip(3))
        u01 = rp(1)*self%values(element_index, ip1(1), ip1(2), ip(3)) &
              + rp1(1)*self%values(element_index, ip(1), ip1(2), ip(3))
        u10 = rp(1)*self%values(element_index, ip1(1), ip(2), ip1(3)) &
              + rp1(1)*self%values(element_index, ip(1), ip(2), ip1(3))
        u11 = rp(1)*self%values(element_index, ip1(1), ip1(2), ip1(3)) &
              + rp1(1)*self%values(element_index, ip(1), ip1(2), ip1(3))

        u0 = rp(2)*u01 + rp1(2)*u00
        u1 = rp(2)*u11 + rp1(2)*u10

        ret = rp(3)*u1 + rp1(3)*u0
    end function

end module

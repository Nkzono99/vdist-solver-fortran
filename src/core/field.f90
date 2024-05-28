module m_field
    implicit none

    type, abstract :: t_VectorField
        integer :: n_elements
    contains
        procedure(vectorField_at), deferred :: at
    end type

    interface
        function vectorField_at(self, position) result(ret)
            import t_VectorField
            class(t_VectorField), intent(in) :: self
            double precision, intent(in) :: position(3)
            double precision :: ret(self%n_elements)
        end function
    end interface

    type, extends(t_VectorField) :: t_VectorFieldGrid
        integer :: nx
        integer :: ny
        integer :: nz
        double precision, allocatable :: values(:, :, :, :)
    contains

        procedure :: at => vectorFieldGrid_at
    end type

    contains

    function new_VectorFieldGrid(n_elements, nx, ny, nz, values) result(obj)
        integer, intent(in) :: n_elements
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        integer, intent(in) :: nz
        double precision, intent(in) :: values(n_elements, nx+1, ny+1, nz+1)
        type(t_VectorFieldGrid) :: obj

        obj%n_elements = n_elements
        obj%nx = obj%nx
        obj%ny = obj%ny
        obj%nz = obj%nz
        allocate(obj%values(n_elements, nx+1, ny+1, nz+1))
        obj%values(:, :, :, :) = values(:, :, :, :)
    end function

    function vectorFieldGrid_at(self, position) result(ret)
        class(t_VectorFieldGrid), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision :: ret(self%n_elements)
    end function

    end module

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
        double precision, intent(in) :: values(n_elements, 0:nx, 0:ny, 0:nz)
        type(t_VectorFieldGrid) :: obj

        obj%n_elements = n_elements
        obj%nx = nx
        obj%ny = ny
        obj%nz = nz
        allocate(obj%values(n_elements, 0:nx, 0:ny, 0:nz))
        obj%values(:, :, :, :) = values(:, :, :, :)
    end function

    function vectorFieldGrid_at(self, position) result(ret)
        class(t_VectorFieldGrid), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision :: ret(self%n_elements)

        double precision :: p(3)
        integer :: ip(3), ip1(3)
        double precision :: rp(3), rp1(3)
        double precision :: u00(6), u01(6), u10(6), u11(6), u0(6), u1(6)

        ! Linear interpolation
        p(:) = position(:)
        ip(:) = int(p(:))
        ip1(:) = ip(:) + 1

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

end module

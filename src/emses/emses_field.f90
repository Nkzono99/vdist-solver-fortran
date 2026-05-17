module m_emses_field
    !! EMSES field interpolation helpers.
    !!
    !! Space-charge E and B are supplied as relocated node values. Accumulated
    !! charge E is supplied on the original staggered electric-field locations
    !! and is interpolated with the same half-cell offsets used by MPIEMSES3D.

    use m_field, only: t_VectorField, t_VectorFieldGrid, new_VectorFieldGrid

    implicit none

    private
    public t_EMSESFieldGrid
    public new_EMSESFieldGrid

    type, extends(t_VectorField) :: t_EMSESFieldGrid
        !! Combined EMSES E/B field.
        type(t_VectorFieldGrid) :: relocated_eb
            !! Relocated space-charge E and B.
        type(t_VectorFieldGrid) :: accumulated_e
            !! Staggered accumulated-charge E.
    contains
        procedure :: at => emsesFieldGrid_at
    end type

contains

    function new_EMSESFieldGrid(nx, ny, nz, relocated_eb_values, accumulated_e_values) result(obj)
        !! Create a combined EMSES field grid.

        integer, intent(in) :: nx
            !! Number of grid cells in the x direction
        integer, intent(in) :: ny
            !! Number of grid cells in the y direction
        integer, intent(in) :: nz
            !! Number of grid cells in the z direction
        double precision, intent(in) :: relocated_eb_values(6, 0:nx, 0:ny, 0:nz)
            !! Relocated E/B values.
        double precision, intent(in) :: accumulated_e_values(3, 0:nx, 0:ny, 0:nz)
            !! Accumulated-charge E values on staggered component grids.
        type(t_EMSESFieldGrid) :: obj
            !! Combined field grid.

        double precision :: offsets(3, 3)

        offsets(:, :) = 0d0
        offsets(1, 1) = 0.5d0
        offsets(2, 2) = 0.5d0
        offsets(3, 3) = 0.5d0

        obj%n_elements = 6
        obj%relocated_eb = new_VectorFieldGrid(6, nx, ny, nz, relocated_eb_values)
        obj%accumulated_e = new_VectorFieldGrid(3, nx, ny, nz, accumulated_e_values, offsets)
    end function

    function emsesFieldGrid_at(self, position) result(ret)
        !! Get the combined E/B value at a particle position.

        class(t_EMSESFieldGrid), intent(in) :: self
            !! Instance of the field grid
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision :: ret(self%n_elements)
            !! Combined E/B value

        double precision :: accumulated_e(3)

        ret(:) = self%relocated_eb%at(position)
        accumulated_e(:) = self%accumulated_e%at(position)
        ret(1:3) = ret(1:3) + accumulated_e(:)
    end function

end module

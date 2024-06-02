module m_dust_charge_simulator
    use m_particle
    use m_simulator
    use m_field

    implicit none

    private
    public t_DustParticle
    public new_DustParticle
    public t_DustChargeSimulator
    public new_DustChargeSimulator

    double precision, parameter :: pi = acos(-1.0d0)
        !! The mathematical constant pi

    type :: t_DustParticle
        type(t_Particle) :: particle

        double precision :: charge
        double precision :: mass
        double precision :: radius

    contains

        procedure :: potential => dustParticle_potential

    end type

    type :: t_DustChargeSimulator
        integer :: nx
            !! Number of grid cells in the x direction
        integer :: ny
            !! Number of grid cells in the y direction
        integer :: nz
            !! Number of grid cells in the z direction
        integer :: nspec
        double precision, allocatable :: temperatures(:)
        class(t_VectorField), allocatable :: currents(:)

        double precision :: jph0 = 0d0

    contains

        procedure :: update => dustChargeSimulator_update

    end type

contains

    function new_DustParticle(charge, mass, radius, position, velocity) result(obj)
        !! Create a new particle object with specified properties.
        double precision, intent(in) :: charge
            !! Charge
        double precision, intent(in) :: mass
            !! Mass
        double precision, intent(in) :: radius
            !! Radius
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision, intent(in) :: velocity(3)
            !! Velocity in 3D space
        type(t_DustParticle) :: obj
            !! A new dust particle object with specified properties

        obj%charge = charge
        obj%mass = mass
        obj%radius = radius
        obj%particle = new_Particle(charge/mass, position, velocity)
    end function

    function dustParticle_potential(self) result(ret)
        class(t_DustParticle), intent(in) :: self
        double precision :: ret

        ! It has been nomalized by epsilon_0 == 1
        ret = self%charge/(4*pi*self%radius)
    end function

    function new_DustChargeSimulator(nx, ny, nz, &
                                     nspec, &
                                     temperatures, &
                                     currents, &
                                     jph0) result(obj)

        integer, intent(in) :: nx
            !! Number of grid points in the x direction
        integer, intent(in) :: ny
            !! Number of grid points in the y direction
        integer, intent(in) :: nz
            !! Number of grid points in the z direction
        integer, intent(in) :: nspec
        double precision, intent(in) :: temperatures(nspec)
        class(t_VectorField), intent(in) :: currents(nspec)
        double precision, intent(in), optional :: jph0

        type(t_DustChargeSimulator) :: obj

        integer :: ispec

        obj%nx = nx
        obj%ny = ny
        obj%nz = nz
        obj%nspec = nspec

        obj%temperatures = temperatures(:)

        obj%currents = currents

        if (present(jph0)) then
            obj%jph0 = jph0
        end if
    end function

    function dustChargeSimulator_update(self, dust, dt) result(dust_new)
        class(t_DustChargeSimulator), intent(in) :: self
        type(t_DustParticle), intent(in) :: dust
        double precision, intent(in) :: dt
        type(t_DustParticle) :: dust_new

        integer :: ispec
        double precision :: net_current
        double precision :: dust_potential

        dust_potential = dust%potential()

        net_current = 0d0

        ispec = 1
        if (ispec <= self%nspec) then
            block
                double precision :: current(3), j0
                double precision :: temperature

                current = self%currents(ispec)%at(dust%particle%position)
                j0 = norm2(current)

                temperature = self%temperatures(ispec)

                net_current = net_current + calculate_electron_current(dust_potential, j0, temperature)
            end block
        end if

        ispec = 2
        if (ispec <= self%nspec) then
            block
                double precision :: current(3), j0
                double precision :: temperature

                current = self%currents(ispec)%at(dust%particle%position)
                j0 = norm2(current)

                temperature = self%temperatures(ispec)

                net_current = net_current + calculate_ion_current(dust_potential, j0, temperature)
            end block
        end if

        ispec = 3
        if (ispec <= self%nspec) then
            block
                double precision :: current(3), j0
                double precision :: temperature

                current = self%currents(ispec)%at(dust%particle%position)
                j0 = norm2(current)

                temperature = self%temperatures(ispec)

                net_current = net_current + calculate_photoelectron_current(dust_potential, self%jph0, j0, temperature)
            end block
        end if

        block
            double precision :: charge_new

            charge_new = dust%charge + net_current*dt

            dust_new = new_DustParticle(charge_new, dust%mass, &
                                        dust%radius, &
                                        dust%particle%position, &
                                        dust%particle%velocity)
        end block
    end function

    function calculate_electron_current(phid, je0, Te) result(ret)
        double precision, intent(in) :: phid
        double precision, intent(in) :: je0
        double precision, intent(in) :: Te
        double precision :: ret

        if (phid >= 0) then
            ret = je0*(1 + phid/Te)
        else
            ret = je0*(1 - phid/Te)
        end if
    end function

    function calculate_ion_current(phid, ji0, Ti) result(ret)
        double precision, intent(in) :: phid
        double precision, intent(in) :: ji0
        double precision, intent(in) :: Ti
        double precision :: ret

        if (phid >= 0) then
            ret = ji0*exp(-phid/Ti)
        else
            ret = ji0*exp(phid/Ti)
        end if
    end function

    function calculate_photoelectron_current(phid, jph0, jphc0, Tph) result(ret)
        double precision, intent(in) :: phid
        double precision, intent(in) :: jph0
        double precision, intent(in) :: jphc0
        double precision, intent(in) :: Tph
        double precision :: ret

        if (phid >= 0) then
            ret = jph0*exp(-phid/Tph)*(1 + phid/Tph) + jphc0*(1 + phid/Tph)
        else
            ret = jph0 + jphc0*(1 - phid/Tph)
        end if
    end function

end module

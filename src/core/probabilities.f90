module m_Probabilities
    !! This module defines different probability distributions.

    implicit none

    double precision :: pi = acos(-1.0d0)
        !! The mathematical constant pi

    type, abstract :: t_Probability
        !! Abstract type for a probability distribution
        double precision :: coefficient = 1d0
            !! Coefficient for the probability distribution
    contains
        procedure(probabilities_at), deferred :: at
            !! Get the probability at a given position and velocity
    end type

    interface
        function probabilities_at(self, position, velocity) result(ret)
            !! Get the probability at a given position and velocity
            import t_Probability
            class(t_Probability), intent(in) :: self
                !! Instance of the probability distribution
            double precision, intent(in) :: position(3)
                !! Position in 3D space
            double precision, intent(in) :: velocity(3)
                !! Velocity in 3D space
            double precision :: ret
                !! Probability at the given position and velocity
        end function
    end interface

    type, extends(t_Probability) :: t_ZeroProbability
        !! A probability distribution that always returns zero
    contains
        procedure :: at => zeroProbability_at
            !! Get the probability (always zero) at a given position and velocity
    end type

    type, extends(t_Probability) :: t_MaxwellianProbability
        !! A Maxwellian probability distribution
        double precision :: locs(3)
            !! Locations (means) of the distribution
        double precision :: scales(3)
            !! Scales (standard deviations) of the distribution
    contains
        procedure :: at => maxwellianProbability_at
            !! Get the probability at a given position and velocity
    end type

    type, extends(t_Probability) :: tp_Probability
        !! A probability distribution that references another distribution
        class(t_Probability), pointer :: ref
            !! Reference to another probability distribution
    contains
        procedure :: at => p_Probability_at
            !! Get the probability from the referenced distribution
    end type

    private
    public t_Probability
    public t_ZeroProbability
    public new_ZeroProbability
    public t_MaxwellianProbability
    public new_MaxwellianProbability
    public tp_Probability

contains

    function new_ZeroProbability() result(obj)
        !! Create a new zero probability distribution object.

        type(t_ZeroProbability) :: obj
            !! A new zero probability distribution object
    end function

    function zeroProbability_at(self, position, velocity) result(ret)
        !! Get the probability (always zero) at a given position and velocity.

        class(t_ZeroProbability), intent(in) :: self
            !! Instance of the zero probability distribution
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision, intent(in) :: velocity(3)
            !! Velocity in 3D space
        double precision :: ret
            !! Probability (always zero) at the given position and velocity

        ret = 0d0
    end function

    function new_MaxwellianProbability(locs, scales, coefficient) result(obj)
        !! Create a new Maxwellian probability distribution object with specified properties.

        double precision, intent(in) :: locs(3)
            !! Locations (means) of the distribution
        double precision, intent(in) :: scales(3)
            !! Scales (standard deviations) of the distribution
        double precision, optional, intent(in) :: coefficient
            !! Coefficient for the probability distribution
        type(t_MaxwellianProbability) :: obj
            !! A new Maxwellian probability distribution object

        obj%locs = locs
        obj%scales = scales

        if (present(coefficient)) then
            obj%coefficient = coefficient
        else
            obj%coefficient = 1d0
        end if
    end function

    function maxwell_pdf(x, mu, sigma) result(ret)
        !! Calculate the Maxwellian probability density function value.

        double precision, intent(in) :: x
            !! Value to evaluate the PDF at
        double precision, intent(in) :: mu
            !! Mean of the distribution
        double precision, intent(in) :: sigma
            !! Standard deviation of the distribution
        double precision :: ret
            !! Resulting PDF value

        ret = 1/sqrt(2*pi*sigma*sigma)*exp(-(x - mu)*(x - mu)/(2*sigma*sigma))
    end function

    function maxwellianProbability_at(self, position, velocity) result(ret)
        !! Get the probability at a given position and velocity using Maxwellian distribution.

        class(t_MaxwellianProbability), intent(in) :: self
            !! Instance of the Maxwellian probability distribution
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision, intent(in) :: velocity(3)
            !! Velocity in 3D space
        double precision :: ret
            !! Probability at the given position and velocity

        ret = self%coefficient
        block
            integer :: i
            do i = 1, 3
                ret = ret*maxwell_pdf(velocity(i), self%locs(i), self%scales(i))
            end do
        end block
    end function

    function p_Probability_at(self, position, velocity) result(ret)
        !! Get the probability from the referenced distribution at a given position and velocity.

        class(tp_Probability), intent(in) :: self
            !! Instance of the referencing probability distribution
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision, intent(in) :: velocity(3)
            !! Velocity in 3D space
        double precision :: ret
            !! Probability at the given position and velocity from the referenced distribution

        ret = self%ref%at(position, velocity)
    end function

end module

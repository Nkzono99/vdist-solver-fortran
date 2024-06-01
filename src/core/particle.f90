module m_particle
    !! This module defines a particle with properties such as position, velocity, and charge-to-mass ratio.

    implicit none

    private
    public t_Particle
    public new_Particle

    type t_Particle
        double precision :: t
            !! Time since simulation began
        double precision :: q_m
            !! Charge to mass ratio
        double precision :: position(3)
            !! Position in 3D space
        double precision :: velocity(3)
            !! Velocity in 3D space
    end type

contains

    function new_Particle(q_m, position, velocity) result(obj)
        !! Create a new particle object with specified properties.
        double precision, intent(in) :: q_m
            !! Charge to mass ratio
        double precision, intent(in) :: position(3)
            !! Position in 3D space
        double precision, intent(in) :: velocity(3)
            !! Velocity in 3D space
        type(t_Particle) :: obj
            !! A new particle object with specified properties

        obj%t = 0d0
        obj%q_m = q_m
        obj%position = position
        obj%velocity = velocity
    end function

end module

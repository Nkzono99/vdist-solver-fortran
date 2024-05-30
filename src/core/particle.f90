module m_particle
    implicit none

    private
    public t_Particle
    public new_Particle

    type t_Particle
        double precision :: t
        double precision :: position(3)
        double precision :: velocity(3)
    end type

contains

    function new_Particle(position, velocity) result(obj)
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        type(t_Particle) :: obj

        obj%t = 0d0
        obj%position = position
        obj%velocity = velocity
    end function

end module

module m_particle
    implicit none

    private
    public t_Particle
    public new_Particle

    type t_Particle
        double precision :: t
        double precision :: q_m
        double precision :: position(3)
        double precision :: velocity(3)
    end type

contains

    function new_Particle(q_m, position, velocity) result(obj)
        double precision, intent(in) :: q_m
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        type(t_Particle) :: obj

        obj%t = 0d0
        obj%q_m = q_m
        obj%position = position
        obj%velocity = velocity
    end function

end module

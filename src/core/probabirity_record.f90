module m_probabirity_record
    use m_particle
    implicit none

    type t_ProbabirityRecord
        logical :: is_valid = .false.
        double precision :: t
        double precision :: probabirity
        type(t_Particle) :: particle
    end type

end module
module m_Probabirities
    implicit none

    double precision :: pi = acos(-1.0d0)

    type, abstract :: t_Probabirity
        double precision :: coefficient = 1d0
    contains
        procedure(probabirities_at), deferred :: at
    end type

    interface
        function probabirities_at(self, position, velocity) result(ret)
            import t_Probabirity
            class(t_Probabirity), intent(in) :: self
            double precision, intent(in) :: position(3)
            double precision, intent(in) :: velocity(3)
            double precision :: ret
        end function
    end interface

    type, extends(t_Probabirity) :: t_ZeroProbabirity
    contains
        procedure :: at => zeroProbabirity_at
    end type

    type, extends(t_Probabirity) :: t_MaxwellianProbabirity
        double precision :: locs(3)
        double precision :: scales(3)
    contains
        procedure :: at => maxwellianProbabirity_at
    end type

    type, extends(t_Probabirity) :: tp_Probabirity
        class(t_Probabirity), allocatable :: ref
    contains
        procedure :: at => p_Probabirity_at
    end type

contains

    function new_ZeroProbabirity() result(obj)
        type(t_ZeroProbabirity) :: obj
    end function

    function zeroProbabirity_at(self, position, velocity) result(ret)
        class(t_ZeroProbabirity), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        double precision :: ret

        ret = 0d0
    end function

    function new_MaxwellianProbabirity(locs, scales, coefficient) result(obj)
        double precision, intent(in) :: locs(3)
        double precision, intent(in) :: scales(3)
        double precision, optional, intent(in) :: coefficient
        type(t_MaxwellianProbabirity) :: obj

        obj%locs = locs
        obj%scales = scales

        if (present(coefficient)) then
            obj%coefficient = coefficient
        else
            obj%coefficient = 1d0
        end if
    end function

    function maxwell_pdf(x, mu, sigma) result(ret)
        double precision, intent(in) :: x
        double precision, intent(in) :: mu
        double precision, intent(in) :: sigma
        double precision :: ret

        ret = 1/sqrt(2*pi*sigma*sigma)*exp((x - mu)*(x - mu)/2*sigma*sigma)
    end function

    function maxwellianProbabirity_at(self, position, velocity) result(ret)
        class(t_MaxwellianProbabirity), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        double precision :: ret

        ret = self%coefficient
        block
            integer :: i
            do i = 1, 3
                ret = ret*maxwell_pdf(velocity(i), self%locs(i), self%scales(i))
            end do
        end block
    end function

    function p_Probabirity_at(self, position, velocity) result(ret)
        class(tp_Probabirity), intent(in) :: self
        double precision, intent(in) :: position(3)
        double precision, intent(in) :: velocity(3)
        double precision :: ret

        ret = self%ref%at(position, velocity)
    end function

end module

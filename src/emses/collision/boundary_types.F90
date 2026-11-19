!> Central registry for collision boundary type names.
module m_collision_boundary_types
  use iso_fortran_env, only: error_unit
  use m_kinds, only: ip
  implicit none
  private

  public :: is_collision_boundary_entry_empty
  public :: is_collision_boundary_type_supported
  public :: is_collision_object_boundary_type
  public :: is_collision_surface_boundary_type
  public :: require_supported_collision_boundary_entry
  public :: require_supported_simple_collision_boundary_type

contains

  pure logical function is_collision_boundary_entry_empty(name)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("", "none")
      is_collision_boundary_entry_empty = .true.
    case default
      is_collision_boundary_entry_empty = .false.
    end select
  end function is_collision_boundary_entry_empty

  pure logical function is_collision_object_boundary_type(name)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("rectangle", "rectangle-corners", "sphere", "circlex", "circley", "circlez", &
          "cuboid", "cylinderx", "cylindery", "cylinderz", &
          "open-cylinderx", "open-cylindery", "open-cylinderz", &
          "diskx", "disky", "diskz", &
          "plane-with-circlex", "plane-with-circley", "plane-with-circlez")
      is_collision_object_boundary_type = .true.
    case default
      is_collision_object_boundary_type = .false.
    end select
  end function is_collision_object_boundary_type

  pure logical function is_collision_surface_boundary_type(name)
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ("flat-surface", "rectangle-hole", "dent-rectangle-hole", &
          "cylinder-hole", "hyperboloid-hole", "ellipsoid-hole", &
          "plane-with-circle-hole")
      is_collision_surface_boundary_type = .true.
    case default
      is_collision_surface_boundary_type = .false.
    end select
  end function is_collision_surface_boundary_type

  pure logical function is_collision_boundary_type_supported(name)
    character(len=*), intent(in) :: name

    is_collision_boundary_type_supported = is_collision_boundary_entry_empty(name) .or. &
                                           is_collision_object_boundary_type(name) .or. is_collision_surface_boundary_type(name)
  end function is_collision_boundary_type_supported

  subroutine require_supported_collision_boundary_entry(name, context, itype)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: context
    integer(kind=ip), intent(in), optional :: itype

    if (is_collision_boundary_type_supported(name)) return

    if (present(itype)) then
      write (error_unit, '(a,a,a,i0,a,a)') "Error.Input: unsupported collision boundary type '", &
        trim(name), "' at boundary_types(", itype, ") while ", trim(context)
    else
      write (error_unit, '(a,a,a,a)') "Error.Input: unsupported collision boundary type '", &
        trim(name), "' while ", trim(context)
    end if
    error stop 1
  end subroutine require_supported_collision_boundary_entry

  subroutine require_supported_simple_collision_boundary_type(name, context)
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: context

    if (is_collision_boundary_entry_empty(name)) return
    if (trim(name) == "complex") return
    if (is_collision_surface_boundary_type(name)) return

    if (is_collision_object_boundary_type(name)) then
      write (error_unit, '(a,a,a)') "Error.Input: object collision boundary type '", &
        trim(name), "' requires boundary_type = 'complex'"
    else
      write (error_unit, '(a,a,a,a)') "Error.Input: unsupported collision boundary_type '", &
        trim(name), "' while ", trim(context)
    end if
    error stop 1
  end subroutine require_supported_simple_collision_boundary_type

end module m_collision_boundary_types

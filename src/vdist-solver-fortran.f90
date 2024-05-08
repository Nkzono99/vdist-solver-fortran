module vdist_solver_fortran
  implicit none
  private

  public :: say_hello
contains
  subroutine say_hello
    print *, "Hello, vdist-solver-fortran!"
  end subroutine say_hello
end module vdist_solver_fortran

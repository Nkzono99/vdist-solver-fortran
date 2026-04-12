module m_kinds
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none

  private
  public :: ip, lp, sp, dp

  !> Alias for 32-bit integer kind (int32), named ip (integer precision)
  integer, parameter :: ip = int32
  !> Alias for 64-bit integer kind (int64), named lp (long precision)
  integer, parameter :: lp = int64

  !> Alias for single-precision real kind (real32), named sp
  integer, parameter :: sp = real32
  !> Alias for double-precision real kind (real64), named dp
  integer, parameter :: dp = real64
end module

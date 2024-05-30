module allcom
    implicit none

    integer, parameter :: max_nspec = 10
    integer, parameter :: nboundary_types = 10
    integer, parameter :: max_nepl = 15

    ! /esorem/
    integer :: emflag = 0

    ! /plasma/
    double precision :: wp(max_nspec)
    double precision :: wc
    double precision :: phixy
    double precision :: phiz

    ! /tmgrid/
    integer :: nx, ny, nz

    ! /system/
    integer :: nspec
    integer :: npbnd(3, max_nspec)

    ! /intp/
    double precision :: qm(max_nspec)
    double precision :: path(max_nspec)
    double precision :: peth(max_nspec)
    double precision :: vdri(max_nspec)
    double precision :: vdthz(max_nspec)
    double precision :: vdthxy(max_nspec)

    ! /ptcond/
    real(kind=8)    :: zssurf = -9999.0d0
    real(kind=8)    :: xlrechole(2), ylrechole(2), zlrechole(2)
    real(kind=8)    :: xurechole(2), yurechole(2), zurechole(2)

    ! Boundary type
    ! "flat-surface"
    ! "rectangle-hole"
    ! "cylinder-hole"
    ! "hyperboloid-hole"
    ! "rectangle"
    ! "circle"
    ! "cuboid"
    ! "disk"
    character(len=30) :: boundary_type = "flat-surface"
    character(len=30) :: boundary_types(nboundary_types) = "none"
    double precision :: cylinder_origin(nboundary_types, 3) = 0.0d0
    double precision :: cylinder_radius(nboundary_types) = 0.0d0
    double precision :: cylinder_height(nboundary_types) = 0.0d0
    !> Minimum radius relative to the entrance radius of the hyperbolic hole (r_min == rcurv*r_max).
    !> If rcurv == 1.0d0, then r_min == r_max.
    !> Else if rcurv == 0.5d0 then r_min == 0.5d0 * r_max
    double precision :: rcurv = 1.0d0
    double precision :: rectangle_shape(nboundary_types, 6) = 0.0d0
    double precision :: sphere_origin(nboundary_types, 3) = 0.0d0
    double precision :: sphere_radius(nboundary_types) = 0.0d0
    double precision :: circle_origin(nboundary_types, 3) = 0.0d0
    double precision :: circle_radius(nboundary_types) = 0.0d0
    double precision :: cuboid_shape(nboundary_types, 6) = 0.0d0
    double precision :: disk_origin(nboundary_types, 3) = 0.0d0
    double precision :: disk_height(nboundary_types) = 0.0d0
    double precision :: disk_radius(nboundary_types) = 0.0d0
    double precision :: disk_inner_radius(nboundary_types) = 0.0d0

    ! /emissn/
    integer :: nflag_emit(max_nspec)
    integer :: nepl(max_nspec)
    double precision :: curf
    integer :: nemd(max_nepl)
    double precision :: curfs(max_nepl)
    double precision :: xmine(max_nepl), xmaxe(max_nepl)
    double precision :: ymine(max_nepl), ymaxe(max_nepl)
    double precision :: zmine(max_nepl), zmaxe(max_nepl)
    double precision :: thetaz(max_nepl), thetaxy(max_nepl)

end module

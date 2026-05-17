module m_emses_simulator_builder
    !! Construction helpers for EMSES simulators.
    !!
    !! Encapsulates the logic to build `t_ESSimulator` from an input namelist
    !! plus grid arrays. Separated from the C-facing API (`m_emses_solver`) so
    !! that builder logic can be tested and evolve independently of the ctypes
    !! bridge.

    use, intrinsic :: iso_c_binding

    use finbound, only: t_Boundary, t_BoundaryList, new_BoundaryList, &
                        t_PlaneXYZ, &
                        new_PlaneX, new_PlaneY, new_PlaneZ, &
                        t_RectangleXYZ, &
                        new_RectangleX, new_RectangleY, new_RectangleZ

    use m_vdsolverf_core
    use m_allcom
    use m_namelist, only: read_namelist
    use m_emses_boundaries, only: create_simple_collision_boundaries
    use m_photoelectron_raycast, only: new_PhotoelectronRaycastProbability, &
                                       t_PhotoelectronRaycastProbability

    use m_maxwell_flux_erf, only: solve_density_from_flux_erf

    implicit none

    private
    public create_simulator
    public destroy_simulator

contains

    subroutine destroy_simulator(simulator)
        !! Release all heap-allocated state owned by the simulator: the
        !! boundary list as well as each probability function (including any
        !! occlusion boundary list held by the raycast photoelectron prob).
        !!
        !! `select type` is used instead of a polymorphic destroy() method
        !! so that concrete probability classes that do not own heap state
        !! (ZeroProbability, MaxwellianProbability) don't need to declare a
        !! no-op override.
        type(t_ESSimulator), intent(inout) :: simulator
        integer :: i

        call simulator%boundaries%destroy

        if (allocated(simulator%probability_functions)) then
            do i = 1, size(simulator%probability_functions)
                if (associated(simulator%probability_functions(i)%ref)) then
                    select type (p => simulator%probability_functions(i)%ref)
                    type is (t_PhotoelectronRaycastProbability)
                        call p%blocking_boundaries%destroy
                    end select
                    deallocate (simulator%probability_functions(i)%ref)
                end if
            end do
            deallocate (simulator%probability_functions)
        end if
    end subroutine

    function create_simulator(inppath, length, lx, ly, lz, ebvalues, ispec, max_probability_types) result(simulator)
        !! Create and initialize a new ES simulator object.

        character(1, c_char), intent(in) :: inppath(*)
            !! Path to the input file
        integer(c_int), value, intent(in) :: length
            !! Length of the input path
        integer(c_int), value, intent(in) :: lx
            !! Number of grid cells in the x direction
        integer(c_int), value, intent(in) :: ly
            !! Number of grid cells in the y direction
        integer(c_int), value, intent(in) :: lz
            !! Number of grid cells in the z direction
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
            !! Electric and magnetic field values
        integer(c_int), value, intent(in) :: ispec
            !! Species index
        integer(c_int), value, intent(in) :: max_probability_types
            !! Maximum number of probability types
        type(t_ESSimulator) :: simulator
            !! Simulator object

        type(t_VectorFieldGrid), target :: eb

        type(t_BoundaryList) :: boundaries

        type(tp_Probability), allocatable :: probability_functions(:)

        integer :: n_probability_functions

        n_probability_functions = 0

        allocate (probability_functions(max_probability_types))

        block
            character(length) :: s
            integer :: i

            do i = 1, length
                s(i:i) = inppath(i)
            end do
            call read_namelist(s)
        end block

        block
            integer :: isdoms(2, 3)
            integer :: boundary_conditions(3)

            isdoms = reshape([[0, lx], [0, ly], [0, lz]], [2, 3])

            call register_inner_boundary_probability(isdoms, ispec, &
                                                     probability_functions, &
                                                     n_probability_functions)

            boundaries = create_simple_collision_boundaries(isdoms, tag=n_probability_functions)

            ! boundaries = new_BoundaryList()

            eb = new_VectorFieldGrid(6, lx, ly, lz, ebvalues(1:6, 1:lx + 1, 1:ly + 1, 1:lz + 1))

            call add_probability_boundaries(boundaries, ispec, n_probability_functions, probability_functions)

            boundary_conditions(:) = npbnd(:, ispec)
            simulator = new_ESSimulator(lx, ly, lz, &
                                        boundary_conditions, &
                                        eb, &
                                        boundaries, &
                                        probability_functions)
        end block
    end function

    subroutine register_inner_boundary_probability(isdoms, ispec, &
                                                   probability_functions, &
                                                   n_probability_functions)
        !! Register the probability function used when a backtraced particle
        !! collides with an internal surface/object.
        !!
        !! Default: `t_ZeroProbability` (absorbing).
        !! When `use_raycast .and. nflag_emit(ispec) == 2`, a
        !! `t_PhotoelectronRaycastProbability` is used instead. A separate
        !! collision-boundary list (built from the same input geometry) is
        !! handed to it as the occlusion test target so the main simulator
        !! boundary list can be destroyed independently.
        integer, intent(in) :: isdoms(2, 3)
        integer, intent(in) :: ispec
        type(tp_Probability), intent(inout) :: probability_functions(:)
        integer, intent(inout) :: n_probability_functions

        type(t_BoundaryList) :: blocking_boundaries
        double precision :: sun_direction(3)
        double precision :: vmean(3), vthermal(3)

        n_probability_functions = n_probability_functions + 1

        if (use_raycast .and. nflag_emit(ispec) == 2) then
            blocking_boundaries = create_simple_collision_boundaries(isdoms)

            vmean = vdri_vector(ispec)
            vthermal = vth_vector(ispec)
            sun_direction = resolve_sun_direction(ispec)

            allocate (probability_functions(n_probability_functions)%ref, &
                      source=new_PhotoelectronRaycastProbability( &
                      vmean, vthermal, sun_direction, blocking_boundaries))
        else
            allocate (probability_functions(n_probability_functions)%ref, &
                      source=new_ZeroProbability())
        end if
    end subroutine

    function resolve_sun_direction(ispec) result(ret)
        !! Direction from the collision point toward the sun.
        !!
        !! Uses `ray_zenith_angle_deg(ispec)` and `ray_azimuth_angle_deg(ispec)`
        !! when those are set (sentinel: 9999d0 means "unset"). Otherwise falls
        !! back to `vdthz(ispec)` / `vdthxy(ispec)`. The returned vector is the
        !! negation of the vdri-style rotation of +z, normalized — i.e. the
        !! opposite of the drift vector per the user-requested convention.
        integer, intent(in) :: ispec
        double precision :: ret(3)

        double precision :: zenith_deg, azimuth_deg
        double precision :: drift_unit(3)
        double precision, parameter :: SENTINEL = 9000d0

        if (ray_zenith_angle_deg(ispec) < SENTINEL) then
            zenith_deg = ray_zenith_angle_deg(ispec)
        else
            zenith_deg = vdthz(ispec)
        end if

        if (ray_azimuth_angle_deg(ispec) < SENTINEL) then
            azimuth_deg = ray_azimuth_angle_deg(ispec)
        else
            azimuth_deg = vdthxy(ispec)
        end if

        drift_unit = [0d0, 0d0, 1d0]
        drift_unit = rot3d_y(drift_unit, -zenith_deg*DEG2RAD)
        drift_unit = rot3d_z(drift_unit, azimuth_deg*DEG2RAD)

        if (norm2(drift_unit) > 0d0) then
            ret = -drift_unit/norm2(drift_unit)
        else
            ret = [0d0, 0d0, -1d0]
        end if
    end function

    subroutine add_probability_boundaries(boundaries, &
                                          ispec, &
                                          n_probability_functions, &
                                          probability_functions)
        !! Add probability boundaries to the simulator.

        type(t_BoundaryList), intent(inout) :: boundaries
            !! Boundary list
        integer, intent(in) :: ispec
            !! Species index
        integer, intent(inout) :: n_probability_functions
            !! Number of probability functions
        type(tp_Probability), intent(inout) :: probability_functions(:)
            !! Array of probability functions

        call add_external_boundaries
        call add_emission_surface(priority=1)

    contains

        subroutine add_external_boundaries
            !! Add external boundaries to the simulator.

            double precision :: vmean(3), vthermal(3)
                !! Mean velocity and thermal velocity
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZ), pointer :: pplane

            integer :: tag_zero
                !! Tag for zero probability boundary
            integer :: tag_vdist
                !! Tag for maxwellian velocity distribution boundary

            vmean(:) = vdri_vector(ispec)
            vthermal(:) = vth_vector(ispec)

            allocate (probability_functions(n_probability_functions + 1)%ref, &
                      source=new_ZeroProbability())
            n_probability_functions = n_probability_functions + 1
            tag_zero = n_probability_functions

            if (nflag_emit(ispec) == 0) then
                allocate (probability_functions(n_probability_functions + 1)%ref, &
                          source=new_MaxwellianProbability(vmean, vthermal))
                n_probability_functions = n_probability_functions + 1
                tag_vdist = n_probability_functions
            else
                tag_vdist = tag_zero
            end if

            if (npbnd(1, ispec) == 2) then ! X-Boundary
                ! X lower boundary
                allocate (pplane)
                pplane = new_PlaneX(0d0)
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)

                ! X higher boundary
                allocate (pplane)
                pplane = new_PlaneX(dble(nx))
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(2, ispec) == 2) then ! Y-Boundary
                ! Y lower boundary
                allocate (pplane)
                pplane = new_PlaneY(0d0)
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)

                ! Y higher boundary
                allocate (pplane)
                pplane = new_PlaneY(dble(ny))
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(3, ispec) == 2) then ! Z-Boundary
                ! Z lower boundary
                allocate (pplane)
                pplane = new_PlaneZ(0d0)
                pbound => pplane
                if (zssurf < 0d0) then
                    pbound%material%tag = tag_vdist
                else
                    pbound%material%tag = tag_zero
                end if
                call boundaries%add_boundary(pbound)

                ! Z higher boundary
                allocate (pplane)
                pplane = new_PlaneZ(dble(nz))
                pbound => pplane
                pbound%material%tag = tag_vdist
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(1, ispec) == 3) then ! X-Boundary
                ! X lower boundary
                allocate (pplane)
                pplane = new_PlaneX(0d0)
                pbound => pplane
                pbound%material%tag = tag_zero
                call boundaries%add_boundary(pbound)

                ! X higher boundary
                allocate (pplane)
                pplane = new_PlaneX(dble(nx))
                pbound => pplane
                pbound%material%tag = tag_zero
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(2, ispec) == 3) then ! Y-Boundary
                ! Y lower boundary
                allocate (pplane)
                pplane = new_PlaneY(0d0)
                pbound => pplane
                pbound%material%tag = tag_zero
                call boundaries%add_boundary(pbound)

                ! Y higher boundary
                allocate (pplane)
                pplane = new_PlaneY(dble(ny))
                pbound => pplane
                pbound%material%tag = tag_zero
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(3, ispec) == 3) then ! Z-Boundary
                ! Z lower boundary
                allocate (pplane)
                pplane = new_PlaneZ(0d0)
                pbound => pplane
                if (zssurf < 0d0) then
                    pbound%material%tag = tag_zero
                else
                    pbound%material%tag = tag_zero
                end if
                call boundaries%add_boundary(pbound)

                ! Z higher boundary
                allocate (pplane)
                pplane = new_PlaneZ(dble(nz))
                pbound => pplane
                pbound%material%tag = tag_zero
                call boundaries%add_boundary(pbound)
            end if
        end subroutine

        subroutine add_emission_surface(priority)
            integer, intent(in) :: priority

            integer :: iepl_start, iepl_end
            integer :: iepl
            integer :: tag_vdist

            if (nepl(ispec) == 0) then
                return
            end if

            if (ispec == 1) then
                iepl_start = 1
            else
                iepl_start = sum(nepl(1:ispec - 1)) + 1
            end if

            iepl_end = sum(nepl(1:ispec))

            do iepl = iepl_start, iepl_end
                block
                    double precision :: vmean(3)
                    double precision :: vthermal(3)
                    double precision :: density, flux

                    vmean = emission_vdri_vector(ispec, iepl)
                    vthermal = emission_vth_vector(ispec, iepl)

                    flux = curf(ispec)
                    if (curfs(iepl) >= 0d0) then
                        flux = curfs(ispec)
                    end if

                    density = solve_density_from_flux_erf( &
                              flux, &
                              vmean(abs(nemd(iepl))), &
                              vthermal(abs(nemd(iepl))) &
                              )

                    allocate (probability_functions(n_probability_functions + 1)%ref, &
                              source=new_MaxwellianProbability(vmean, vthermal, density))
                    n_probability_functions = n_probability_functions + 1
                    tag_vdist = n_probability_functions
                end block

                block
                    class(t_RectangleXYZ), pointer :: prect
                    class(t_Boundary), pointer :: pbound
                    double precision :: origin(3), wx, wy, wz

                    origin(:) = [xmine(iepl), ymine(iepl), zmine(iepl)]
                    wx = xmaxe(iepl) - xmine(iepl)
                    wy = ymaxe(iepl) - ymine(iepl)
                    wz = zmaxe(iepl) - zmine(iepl)

                    if (abs(nemd(iepl)) == 1) then
                        allocate (prect, source=new_RectangleX(origin, wy, wz))
                    else if (abs(nemd(iepl)) == 2) then
                        allocate (prect, source=new_RectangleY(origin, wz, wx))
                    else if (abs(nemd(iepl)) == 3) then
                        allocate (prect, source=new_RectangleZ(origin, wx, wy))
                    end if
                    prect%priority = priority

                    pbound => prect
                    pbound%material%tag = tag_vdist
                    call boundaries%add_boundary(pbound)
                end block
            end do
        end subroutine

    end subroutine

end module

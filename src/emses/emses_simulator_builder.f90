module m_emses_simulator_builder
    !! Construction helpers for EMSES simulators.
    !!
    !! Encapsulates the logic to build `t_ESSimulator` and `t_DustChargeSimulator`
    !! from an input namelist plus grid arrays. Separated from the C-facing API
    !! (`m_emses_solver`) so that builder logic can be tested and evolve
    !! independently of the ctypes bridge.

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
    use m_photoelectron_raycast, only: new_PhotoelectronRaycastProbability

    use m_maxwell_flux_erf, only: solve_density_from_flux_erf

    implicit none

    private
    public create_simulator
    public create_dust_charge_simulator

contains

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
            sun_direction = -vmean

            allocate (probability_functions(n_probability_functions)%ref, &
                      source=new_PhotoelectronRaycastProbability( &
                      vmean, vthermal, sun_direction, blocking_boundaries))
        else
            allocate (probability_functions(n_probability_functions)%ref, &
                      source=new_ZeroProbability())
        end if
    end subroutine

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

    function create_dust_charge_simulator(inppath, length, lx, ly, lz, nspecies, current_values, jph0) result(simulator)
        !! Create and initialize a new dust charge simulator object.

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
        integer(c_int), value, intent(in) :: nspecies
            !! Number of species
        real(c_double), intent(in) :: current_values(3*nspecies, lx + 1, ly + 1, lz + 1)
            !! Electric and magnetic field values
        real(c_double), intent(in) :: jph0
            !! Photoelectrons current
        type(t_DustChargeSimulator) :: simulator
            !! Simulator object

        block
            character(length) :: s
            integer :: i
            do i = 1, length
                s(i:i) = inppath(i)
            end do
            call read_namelist(s)
        end block

        block
            double precision, allocatable :: temperatures(:)
            type(t_VectorFieldGrid) :: currents

            allocate (temperatures(nspecies))
            temperatures = path(1:nspecies)*path(1:nspecies)/abs(qm(1:nspecies)) ! [eV in EMSES-U]

            currents = new_VectorFieldGrid(3*nspecies, lx, ly, lz, current_values(:, :, :, :))

            simulator = new_DustChargeSimulator(lx, ly, lz, nspecies, temperatures, currents, jph0)
        end block
    end function

end module

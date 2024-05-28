module m_emses_solver
    use m_simulator
    use m_particle
    use m_boundary_list
    use emses_boundaries
    use allcom
    use m_namelist
    use iso_c_binding
    use m_field
    use m_vector
    implicit none

contains

    subroutine get_probabirities(inppath, &
                                 nx, ny, nz, &
                                 efvalues, &
                                 bfvalues, &
                                 ispec, &
                                 npcls, &
                                 positions, &
                                 velocities, &
                                 dt, &
                                 max_step, &
                                 max_probabirity_types, &
                                 return_probabirities, &
                                 return_positions, &
                                 return_velocities &
                                 ) bind(c)
        character(1, c_char), intent(in) :: inppath(:)
        integer(c_int), intent(in) :: nx
        integer(c_int), intent(in) :: ny
        integer(c_int), intent(in) :: nz
        real(c_double), intent(in) :: efvalues(1, nx + 1, ny + 1, nz + 1)
        real(c_double), intent(in) :: bfvalues(1, nx + 1, ny + 1, nz + 1)
        integer(c_int), intent(in) :: ispec
        integer(c_int), intent(in) :: npcls
        real(c_double), intent(in) :: positions(3, npcls)
        real(c_double), intent(in) :: velocities(3, npcls)
        real(c_double), intent(in) :: dt
        integer(c_int), intent(in) :: max_step
        integer(c_int), optional, intent(in) :: max_probabirity_types
        real(c_double), intent(out) :: return_probabirities(npcls)
        real(c_double), intent(out) :: return_positions(3, npcls)
        real(c_double), intent(out) :: return_velocities(3, npcls)

        type(t_ESSimulator) :: simulator
        type(t_VectorFieldGrid), target :: ef, bf

        type(t_BoundaryList) :: boundaries
        type(tp_Probabirity), allocatable :: probabirity_functions(:)
        integer :: n_probabirity_functions = 0
        integer :: ipcl

        call read_namelist(inppath)

        if (present(max_probabirity_types)) then
            allocate (probabirity_functions(max_probabirity_types))
        else
            allocate (probabirity_functions(100))
        end if

        block
            integer :: isdoms(2, 3)
            integer :: boundary_conditions(3)
            class(t_VectorField), pointer :: pef, pbf

            isdoms = reshape([[0, nx], [0, ny], [0, nz]], [2, 3])
            boundaries = create_simple_collision_boundaries(isdoms, tag=0)

            probabirity_functions(0)%ref = new_ZeroProbabirity()
            n_probabirity_functions = n_probabirity_functions + 1

            ef = new_VectorFieldGrid(3, nx, ny, nz, efvalues(:, :, :, :))
            bf = new_VectorFieldGrid(3, nx, ny, nz, bfvalues(:, :, :, :))
            pef => ef
            pbf => bf

            call add_probabirity_boundaries(boundaries, ispec, n_probabirity_functions, probabirity_functions)

            simulator = new_ESSimulator(nx, ny, nz, &
                                        qm(ispec), &
                                        boundary_conditions, &
                                        pef, &
                                        pbf, &
                                        boundaries, &
                                        probabirity_functions)
        end block

        do ipcl = 1, npcls
            block
                type(t_ProbabirityRecord) :: record
                type(t_Particle) :: particle

                particle = new_Particle(positions(:, ipcl), velocities(:, ipcl))
                record = simulator%calculate_probabirity(particle, dt, max_step)

                if (record%is_valid) then
                    return_probabirities(ipcl) = record%probabirity
                    return_positions(:, ipcl) = record%particle%position(:)
                    return_velocities(:, ipcl) = record%particle%velocity(:)
                else
                    return_probabirities(ipcl) = -1.0d0
                end if
            end block
        end do
    end subroutine

    subroutine add_probabirity_boundaries(boundaries, &
                                          ispec, &
                                          n_probabirity_functions, &
                                          probabirity_functions)
        type(t_BoundaryList), intent(inout) :: boundaries
        integer, intent(in) :: ispec
        integer, intent(inout) :: n_probabirity_functions
        type(tp_Probabirity), intent(inout) :: probabirity_functions(:)

        call add_external_boundaries

    contains

        subroutine add_external_boundaries
            double precision :: vmean(3), vthermal(3)
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZ), pointer :: pplane

            vmean(:) = vdri_vector(ispec)
            vthermal(:) = vth_vector(ispec)

            probabirity_functions(n_probabirity_functions + 1)%ref = &
                new_MaxwellianProbabirity(vmean, vthermal)
            n_probabirity_functions = n_probabirity_functions + 1

            if (npbnd(1, ispec) == 2) then ! X-Boundary
                ! X lower boundary
                allocate (pplane)
                pplane = new_PlaneX(0d0)
                pbound => pplane
                pbound%material%tag = n_probabirity_functions
                call boundaries%add_boundary(pbound)

                ! X higher boundary
                allocate (pplane)
                pplane = new_PlaneX(dble(nx))
                pbound => pplane
                pbound%material%tag = n_probabirity_functions
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(2, ispec) == 2) then ! Y-Boundary
                ! Y lower boundary
                allocate (pplane)
                pplane = new_PlaneY(0d0)
                pbound => pplane
                pbound%material%tag = n_probabirity_functions
                call boundaries%add_boundary(pbound)

                ! Y higher boundary
                allocate (pplane)
                pplane = new_PlaneY(dble(ny))
                pbound => pplane
                pbound%material%tag = n_probabirity_functions
                call boundaries%add_boundary(pbound)
            end if

            if (npbnd(3, ispec) == 2) then ! Z-Boundary
                ! Z lower boundary
                if (zssurf < 0d0) then
                    allocate (pplane)
                    pplane = new_PlaneZ(0d0)
                    pbound => pplane
                    pbound%material%tag = n_probabirity_functions
                    call boundaries%add_boundary(pbound)
                end if

                ! Z higher boundary
                allocate (pplane)
                pplane = new_PlaneZ(dble(nz))
                pbound => pplane
                pbound%material%tag = n_probabirity_functions
                call boundaries%add_boundary(pbound)
            end if
        end subroutine

    end subroutine

    function vdri_vector(ispec) result(ret)
        integer, intent(in) :: ispec
        double precision :: ret(3)

        ret(:) = [0d0, 0d0, vdri(ispec)]
        ret(:) = rot3d_y(ret, vdthz(ispec))
        ret(:) = rot3d_z(ret, vdthxy(ispec))
    end function

    function vth_vector(ispec) result(ret)
        integer, intent(in) :: ispec
        double precision :: ret(3)

        ret(:) = [peth(ispec), peth(ispec), path(ispec)]
        ret(:) = rot3d_y(ret, phiz)
        ret(:) = rot3d_z(ret, phixy)
    end function

end module

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
                                 length, &
                                 lx, ly, lz, &
                                 ebvalues, &
                                 ispec, &
                                 npcls, &
                                 positions, &
                                 velocities, &
                                 dt_multiplier, &
                                 max_step, &
                                 max_probabirity_types, &
                                 return_probabirities, &
                                 return_positions, &
                                 return_velocities &
                                 ) bind(c)
        character(1, c_char), intent(in) :: inppath(*)
        integer(c_int), value, intent(in) :: length
        integer(c_int), value, intent(in) :: lx
        integer(c_int), value, intent(in) :: ly
        integer(c_int), value, intent(in) :: lz
        real(c_double), intent(in) :: ebvalues(6, lx + 1, ly + 1, lz + 1)
        integer(c_int), value, intent(in) :: ispec
        integer(c_int), value, intent(in) :: npcls
        real(c_double), intent(in) :: positions(3, npcls)
        real(c_double), intent(in) :: velocities(3, npcls)
        real(c_double), value, intent(in) :: dt_multiplier
        integer(c_int), value, intent(in) :: max_step
        integer(c_int), value, intent(in) :: max_probabirity_types
        real(c_double), intent(out) :: return_probabirities(npcls)
        real(c_double), intent(out) :: return_positions(3, npcls)
        real(c_double), intent(out) :: return_velocities(3, npcls)

        type(t_ESSimulator) :: simulator
        type(t_VectorFieldGrid), target :: eb

        type(t_BoundaryList) :: boundaries
        type(tp_Probabirity), allocatable :: probabirity_functions(:)
        integer :: n_probabirity_functions = 0
        integer :: ipcl

        block
            character(length) :: s
            integer :: i
            do i = 1, length
                s(i:i) = inppath(i)
            end do
            call read_namelist(s)
        end block

        allocate (probabirity_functions(max_probabirity_types))

        block
            integer :: isdoms(2, 3)
            integer :: boundary_conditions(3)
            class(t_VectorField), pointer :: peb

            allocate (probabirity_functions(n_probabirity_functions + 1)%ref, source=new_ZeroProbabirity())
            n_probabirity_functions = n_probabirity_functions + 1

            isdoms = reshape([[0, lx], [0, ly], [0, lz]], [2, 3])
            boundaries = create_simple_collision_boundaries(isdoms, tag=n_probabirity_functions)

            eb = new_VectorFieldGrid(6, lx, ly, lz, &
                                     ebvalues(1:6, 1:lx + 1, 1:ly + 1, 1:lz + 1))
            peb => eb

            call add_probabirity_boundaries(boundaries, ispec, n_probabirity_functions, probabirity_functions)

            boundary_conditions(:) = npbnd(:, ispec)
            simulator = new_ESSimulator(lx, ly, lz, &
                                        qm(ispec), &
                                        boundary_conditions, &
                                        peb, &
                                        boundaries, &
                                        probabirity_functions)
        end block

        do ipcl = 1, npcls
            ! print *, ipcl, '/', npcls
            block
                type(t_ProbabirityRecord) :: record
                type(t_Particle) :: particle

                particle = new_Particle(positions(:, ipcl), velocities(:, ipcl))
                record = simulator%calculate_probabirity(particle, dt_multiplier*dt, max_step)

                if (record%is_valid) then
                    return_probabirities(ipcl) = record%probabirity
                    return_positions(:, ipcl) = record%particle%position(:)
                    return_velocities(:, ipcl) = record%particle%velocity(:)
                else
                    return_probabirities(ipcl) = -1.0d0
                end if
            end block
        end do

        call boundaries%destroy
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
        call add_external_boundaries

    contains

        subroutine add_external_boundaries
            double precision :: vmean(3), vthermal(3)
            class(t_Boundary), pointer :: pbound
            type(t_PlaneXYZ), pointer :: pplane

            integer :: tag_zero
            integer :: tag_vdist

            vmean(:) = vdri_vector(ispec)
            vthermal(:) = vth_vector(ispec)

            allocate (probabirity_functions(n_probabirity_functions + 1)%ref, &
                      source=new_ZeroProbabirity())
            n_probabirity_functions = n_probabirity_functions + 1
            tag_zero = n_probabirity_functions

            allocate (probabirity_functions(n_probabirity_functions + 1)%ref, &
                      source=new_MaxwellianProbabirity(vmean, vthermal))
            n_probabirity_functions = n_probabirity_functions + 1
            tag_vdist = n_probabirity_functions

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

program test_public_api
    !! Compile-time regression guard for the public module surface.
    !!
    !! Each `use ... only:` clause fails to compile if the corresponding
    !! symbol is removed from the module's public list, so successfully
    !! reaching runtime here is itself the positive test result.

    use m_emses_solver, only: es_get_probabilities => get_probabilities, &
                              es_get_backtraces => get_backtraces
    use m_emses_simulator_builder, only: b_create_simulator => create_simulator, &
                                         b_create_dust_charge_simulator => create_dust_charge_simulator
    use m_vdsolverf, only: u_get_backtraces => get_backtraces, &
                           u_get_probabilities => get_probabilities

    implicit none

    print *, "test_public_api: all expected public symbols are importable."
end program

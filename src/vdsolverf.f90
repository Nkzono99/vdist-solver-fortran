module m_vdsolverf
    !! Umbrella module re-exporting the public API of vdist-solver-fortran.
    use m_emses_solver, only: estimate_velocity_range_map, get_probabilities, &
                              get_probabilities_octree, get_backtraces
end module m_vdsolverf

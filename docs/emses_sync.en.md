# Synchronising with EMSES

> Lang: [日本語](emses_sync.md) | **English**

This solver replays MPIEMSES3D output, so its geometry, particle sources and
equations of motion must match EMSES. EMSES code is not shared as a library;
vdist follows it in two ways:

- **Geometry is copied.** The internal-boundary builders are byte-identical copies
  of files in EMSES `src/physics/collision/`.
- **Sources and motion are written here.** The data layouts differ, so the EMSES
  formulas are re-implemented in vdist code and tests check them against the EMSES
  definitions.

The last EMSES commit reviewed is recorded only in `SYNCED_COMMIT` in
[`scripts/check-emses-sync.sh`](../scripts/check-emses-sync.sh).

## Copied files

| vdist | MPIEMSES3D |
|---|---|
| `src/emses/collision/objects.F90` | `src/physics/collision/objects.F90` |
| `src/emses/collision/surfaces.F90` | `src/physics/collision/surfaces.F90` |
| `src/emses/collision/boundary_types.F90` | `src/physics/collision/boundary_types.F90` |
| `src/emses/collision/boundary_helpers.F90` | `src/physics/collision/boundary_helpers.F90` |
| `src/emses/collision/collision_accelerator.F90` | `src/physics/collision/collision_accelerator.F90` |

Do not edit these files here; fix them in EMSES and copy them again. Two files are
vdist-specific:

- `src/emses/collision/collision_boundary_input_config.F90` fills the same type as
  EMSES `src/main/config/collision_boundary_input_config.F90` from the namelist
  variables in `m_allcom`.
- `src/emses/collision/emses_boundaries.F90` assigns the probability tag.

Keep these EMSES files within 132 columns so they compile with the gfortran default.

## Re-implemented behaviour

| Behaviour | vdist | EMSES reference |
|---|---|---|
| Outer-boundary drift and field-frame thermal spread | `vdri_vector`, `vth_vector`, `thermal_basis` (`src/emses/allcom.f90`) | `initialization_simulation_setup.F90` (`calculate_drift_velocities`), `domain_boundary_source.F90`, `m_thermal_frame.F90` |
| Surface-source mean and widths, patch current and incidence | `surface_source_parameters`, `emission_vdri_vector`, `patch_flux` | `surface_emission_source.F90`, `direct_surface_emitter.F90` |
| Raycast photoelectron flux, normal, tangents and sun direction | `src/emses/photoelectron_raycast.f90`, `resolve_sun_direction` | `raycast_surface_emitter.F90`, `raycast_sampler.F90` |
| One-sided flux moment | `src/utils/maxwell_flux_erf.f90` | `m_velocity_distribution.F90` (`velocity_flux_rate`) |
| Boris update and grid interpolation | `src/core/simulator.f90`, `src/core/field.f90` | `particle_solver.F90`, `m_trilinear_interp.F90` |
| Namelist keys | `src/emses/namelist.f90`, `vdsolverf/emses/tmpolary_input.py` | `namelist_reader.F90` |

The formulas are summarised in [Physics](physics.en.md#supported-source-distributions).
EMSES settings that vdist cannot reproduce raise `ValueError` instead of being
ignored. Python (`vdsolverf/emses/source_support.py`) checks the original EMSES
namelist; Fortran (`check_emses_input`) checks the temporary vdist namelist.

## Procedure

```sh
scripts/check-emses-sync.sh /path/to/MPIEMSES3D origin/main
```

1. If a mirrored file shows `DIFFERS`, copy the EMSES file and run `fpm test`. For
   new shapes, also update `check_geometry` (`src/emses/emses_input_validation.f90`)
   and the namelist keys.
2. Read each commit listed under "EMSES commits since ..." and port the changes that
   affect the behaviour in the table. Add a test for every behaviour change.
3. After the review, set `SYNCED_COMMIT` to the EMSES commit you compared against.

On KUDPC, build and test on compute nodes.

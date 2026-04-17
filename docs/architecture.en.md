# Architecture

> Lang: [日本語](architecture.md) | **English**

`vdist-solver-fortran` is a Fortran numerics core exposed to Python through
a `ctypes` bridge. This document maps the repository layout to the modules
that actually run at query time.

## Top-level layout

```
src/                  Fortran sources (built by fpm)
  core/               Generic physics and solver infrastructure
  emses/              EMSES-specific boundary building and C API
  utils/              Small shared helpers
  vdsolverf.f90       Umbrella module (m_vdsolverf) re-exporting the C API
vdsolverf/            Python package
  core/               Dataclasses (Particle, DustParticle, PhaseGrid)
  emses/              ctypes wrapper, temporary-input builder, geotype helpers
fpm.toml              fpm build configuration (shared library, test auto-discovery)
Makefile              Wraps `fpm install` + platform-specific shared-lib linking
setup.py              Python build; shells out to `make`
test/                 fpm test programs + m_test_helpers module
docs/                 This directory
```

## Fortran module graph

```
m_vdsolverf                           (src/vdsolverf.f90)
  └── m_emses_solver                  (src/emses/emses_solver.f90)      -- C API (bind(c))
        └── m_emses_simulator_builder (src/emses/emses_simulator_builder.f90)
              ├── m_allcom            (src/emses/allcom.f90)            -- namelist globals
              ├── m_namelist          (src/emses/namelist.f90)          -- file IO
              ├── m_emses_boundaries  (src/emses/collision/*.F90)       -- geometry build
              ├── m_photoelectron_raycast (src/emses/photoelectron_raycast.f90)
              └── m_vdsolverf_core    (src/core/vdsolverf_core.f90)     -- aggregate core re-export
                    ├── m_particle
                    ├── m_field
                    ├── m_probabilities
                    ├── m_simulator
                    ├── m_dust_charge_simulator
                    └── m_solver
```

`m_emses_solver` is the only module that exposes C symbols. The builder
module handles simulator construction (including the raycast probability
wiring) and cleanup via `destroy_simulator`. Everything reachable from the
C API lives through the umbrella, so Python consumers see exactly three
entry points: `get_backtraces`, `get_probabilities`, `get_backtrace_dust`.

## Python package

```
vdsolverf.core         Particle, DustParticle, PhaseGrid dataclasses
vdsolverf.emses.wrapper
  _load_dll(...)       Resolves the platform-specific shared library
  get_backtrace(...)   Single-particle convenience wrapper
  get_backtraces(...)  Multi-particle ctypes call
  get_probabilities(...)
  get_dust_backtrace(...)
  create_relocated_ebvalues / create_relocated_current_values
                       Assembles EB/current field arrays from emout
vdsolverf.emses.tmpolary_input
  TempolaryInput       Context manager writing a minimal plasma-vdsolverf.inp
                       from emout data, then cleaning it up on exit.
  TMP_INP_KEYS         Whitelist of namelist keys forwarded to Fortran.
vdsolverf.emses.geotype
  Converts geotype primitives into boundary-type / boundary-shape tuples.
```

## Crossing the boundary

1. Python gathers EB fields (and currents for dust mode) into contiguous
   `numpy` arrays matching the Fortran shapes declared on the `bind(c)`
   subroutines — see
   [`.claude/rules/fortran-python-interop.md`](../.claude/rules/fortran-python-interop.md)
   for the argument-alignment rules.
2. `TempolaryInput` writes the filtered namelist to
   `data.directory / plasma-vdsolverf.inp`.
3. `_load_dll` resolves `libvdist-solver-fortran.so` / `.dylib` / `.dll`.
4. The ctypes call enters one of the three `bind(c)` entry points, which
   builds the simulator, runs the backtrace loop, and calls
   `destroy_simulator` on exit.
5. The temporary namelist file is deleted when `TempolaryInput.__exit__`
   fires.

## Extension points

- **New probability function:** create a `type, extends(t_Probability)` in
  a new module under `src/core/` or `src/emses/`, allocate it in
  `register_inner_boundary_probability` (for inner boundaries) or inside
  `add_probability_boundaries` (for outer planes), and teach
  `destroy_simulator` about any heap-owned state via an extra `select type`
  branch.
- **New namelist key:** append to the appropriate group in
  `src/emses/namelist.f90`, add the key to `TMP_INP_KEYS` in
  `vdsolverf/emses/tmpolary_input.py`, and document it in
  [Namelist reference](namelist.en.md).
- **New C entry point:** add a `bind(c)` subroutine to
  `src/emses/emses_solver.f90`, mirror its signature in `wrapper.py`
  (`argtypes`, `restype`, call site), and run
  [`sync-wrapper-interface`](../.claude/skills/sync-wrapper-interface/SKILL.md)
  to audit the alignment.

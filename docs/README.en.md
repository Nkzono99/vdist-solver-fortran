# Documentation

> Lang: [日本語](README.md) | **English**

Reference material for `vdist-solver-fortran`. For a short project overview
and install instructions see the [top-level README](../README.en.md).

## Table of contents

| Document | What's inside |
|---|---|
| [Usage](usage.en.md) | Python-side recipes: single backtrace, multi-particle backtrace, phase-space probability solver |
| [Velocity range autorange](autorange.en.md) | Per-cell range estimation with `estimate_velocity_range_map`, validation, diagnostics, and performance notes |
| [Namelist reference](namelist.en.md) | Supported EMSES `plasma.inp` groups (`/ptcond/`, `/emissn/`), with raycast-photoelectron parameters |
| [Physics](physics.en.md) | Liouville theorem, Maxwellian surface emission, raycast photoelectron model |
| [Architecture](architecture.en.md) | Repository layout, Fortran/Python boundary, module responsibilities |
| [Development](development.en.md) | Build, test, and extension guide (`fpm`, `.venv`, adding probabilities, CI checklist) |

## Audience

- **Researchers** running EMSES post-processing jobs — start with
  [Usage](usage.en.md), [Velocity range autorange](autorange.en.md), and
  the [Namelist reference](namelist.en.md).
- **Contributors** editing Fortran or Python — skim
  [Architecture](architecture.en.md) then read [Development](development.en.md).
- **Readers checking the physics** — [Physics](physics.en.md) derives what
  each probability function computes.

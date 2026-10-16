# Changelog

All notable changes to this project are documented in this file.

## [Unreleased]

### Added

- Added field-bundle support for velocity-range estimation and validation, plus
  Python-free `mode='velocity_range'` batches calling the same autorange core.
  Added `write_range_batch` / `read_range_batch_results`, C-order range arrays,
  status/settings/timing summaries and resume checks over all input contents.
  MPI range jobs compute only on rank 0 with OpenMP.
- Added public `write_batch_results` for prepared dust-current batches, including
  atomic completion markers and int64 diagnostic counts.
- Added the optional serial/MPI `vdsolverf-batch` executable, dynamic chunk
  dispatch, OpenMP target processing, atomic completion markers and resumable
  npy outputs, plus Python `write_batch` / `read_batch_results` helpers.
- Added Fortran npy v1/v2 C-order little-endian float64/int32/int64 I/O,
  standalone field-bundle loading, and serial/library/resume regressions.
- Added Fortran unit-density dust OML flux moments with terminal-box trapezoid
  integration, sampled variation estimates, bounded per-thread scratch, and the
  `get_dust_currents` Python API / `DustCurrentResult` model.
- Added `leaf_is_terminal` octree diagnostics without changing existing C APIs.
- Added portable field bundle writers/readers (npy v1.0, little-endian float64,
  C order) for reusing relocated fields in octree and dust-current evaluations.
- Added analytic Maxwellian, native ctypes, OpenMP, quadrature, and bundle
  regression coverage and bilingual dust-current documentation.

### Fixed

- Fixed emission-patch density normalisation: use the erf helper's argument
  order, width `sqrt(2)*sigma`, and drift projected along signed `nemd`.
  Emission-surface probabilities change scale; at zero drift with the same
  previously positive normal spread (e.g. +/-Z patches), they increase by
  `sqrt(2*pi)` (about 2.5066) relative to the old swapped-argument call.
  Normal/tangential Gaussian spreads now use `path`/`peth` on every patch axis.
  This also restores X/Y sources whose former signed thermal width gave zero density.
- Prevented installed compiler-specific modules under `include/` from shadowing
  the current fpm build's modules when changing compiler flags or toolchains.

## [1.7.1] - 2026-10-12

### Fixed

- Fixed Maxwellian and raycast photoelectron source distributions in oblique
  magnetic fields by evaluating velocity deviations in the principal thermal
  frame. Isotropic distributions remain rotation invariant, and anisotropic
  distributions retain their full covariance. Autorange support uses Cartesian
  marginal standard deviations (#34).
- Fixed boundary collision times to interpolate the actual signed step
  duration, including adaptive time steps (#35).
- Fixed autorange deposits across periodic boundaries by interpolating the
  unwrapped trajectory before mapping each sample into the periodic domain.
  Multiple periodic axes and multiple crossings per step are supported (#36).
- Initialized unused probability pointers to prevent simulator cleanup crashes.
- Guarded optional surface-coverage arguments before reading their values,
  avoiding compiler-dependent crashes during boundary construction.

### Added

- Added native Fortran regressions and Python shared-library integration tests
  for source rotations, collision times, and periodic deposits. CI runs all
  Python tests against the library built from the source distribution.

## [1.7.0] - 2026-10-06

### Added

- Added independent `use_electric_field` and `use_magnetic_field` options to
  backtrace, probability, octree, autorange, validation, and MPI/srun APIs.
  Both default to `True`; disabled fields are not loaded from EMSES outputs.
  Electric-field opt-out includes space-charge and accumulated-charge fields,
  and magnetic-field opt-out includes the background field.
- Added regression coverage for signed gyromotion, magnetic speed and parallel
  velocity conservation, nonuniform-field inverse steps, periodic crossings,
  and field-option forwarding through Python and MPI/srun wrappers.

### Fixed

- Fixed the Boris magnetic rotation denominator to use the squared vector norm,
  preserving speed in oblique magnetic fields.
- Fixed inverse Boris updates to recover position before evaluating the field
  and undoing the velocity update, matching MPIEMSES3D's forward update order.
  Forward tracing uses the updated velocity for the drift, and field sampling
  wraps periodic crossings into the physical domain.

## [1.6.0] - 2026-07-09

### Added

- Added Fortran-backed `estimate_velocity_range_map` for EMSES per-cell
  velocity-range estimation, including progress logging, sparse accumulation,
  validation support, cell-level access, persistence, and usage documentation.
- Added Fortran-backed `get_probabilities_octree` for adaptive sparse velocity
  space probability evaluation at selected spatial points.
- Added `VelocityOctreeResult` helpers for per-spatial sample access,
  per-spatial leaf access, compact leaf sample slicing, and result diagnostics.
- Added octree and autorange documentation covering algorithms, KUDPC execution
  examples, status handling, visualization, and tuning.

### Changed

- Reduced autorange accumulator memory by using sparse per-thread accumulation
  instead of dense per-thread full-grid buffers.
- Refined octree diagnostics so compact leaf sample spans point into the
  compact returned sample arrays and actual sample/leaf counts are included in
  metadata.
- Expanded Python particle and phase-grid helper coverage.

### Fixed

- Fixed EMSES probability edge cases so invalid probability records retain final
  state information and no-collision solver results are initialized.
- Fixed wrapper helper contracts and added regression tests for octree input
  bounds, output compaction, and capacity/status handling.
- Fixed emission source handling to use per-surface emission current.

## [1.5.0] - 2026-05-18

### Added

- Added raycast photoelectron probability support, including namelist exposure
  and Python wrapper documentation.
- Added accumulated-charge electric-field separation using `phibk` /
  `phibksp*`, with space-charge field relocation aligned to MPIEMSES3D.
- Added optional particle-parallel MPI backend and `srun` launcher wrappers.
- Added PyPI source-distribution packaging with Fortran sources, build
  profiles, and GitHub Actions for CI and trusted publishing.
- Added repository-local agent skills and bilingual documentation structure.

### Changed

- Reduced Python electric-field assembly temporaries for lower memory pressure.
- Split Fortran tests into focused fpm test programs and refactored EMSES
  simulator construction into dedicated modules.
- Renamed internal Fortran modules toward consistent `m_*` naming.
- Updated README installation guidance to prefer `pip install
  vdist-solver-fortran`.

### Removed

- Removed dust-simulator-related code now maintained in a separate repository.

### Fixed

- Fixed temporary-input and dust-charging issues before dust code removal.
- Fixed photoelectron density initialization and `curfs` handling.

[1.7.1]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.7.0...v1.7.1
[1.7.0]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.6.0...v1.7.0
[1.6.0]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.5.0...v1.6.0
[1.5.0]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.4.2...v1.5.0

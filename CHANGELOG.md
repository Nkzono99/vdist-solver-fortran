# Changelog

All notable changes to this project are documented in this file.

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

[1.6.0]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.5.0...v1.6.0
[1.5.0]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.4.2...v1.5.0

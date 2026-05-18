# Changelog

All notable changes to this project are documented in this file.

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

[1.5.0]: https://github.com/Nkzono99/vdist-solver-fortran/compare/v1.4.2...v1.5.0

# vdist-solver-fortran

> Lang: [日本語](README.md) | **English**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14018863.svg)](https://doi.org/10.5281/zenodo.14018863)
[![CI](https://github.com/Nkzono99/vdist-solver-fortran/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Nkzono99/vdist-solver-fortran/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/vdist-solver-fortran)](https://pypi.org/project/vdist-solver-fortran/)

Velocity distribution solver for Python, implemented in Fortran.

The core computes phase-space probability densities by combining
backward-in-time particle tracing over an EMSES simulation with
boundary-tagged source distributions (Maxwellian, raycast photoelectron,
absorbing). The shared library is driven from Python through a `ctypes`
wrapper.

## Requirements

- `gfortran`
- `make`
- `fpm`
- Python 3.7+ (development uses 3.12 in `.venv/`)

## Install

PyPI installation is recommended.  During pip builds, `make install` builds
the Fortran shared library and bundles it into the Python package.

> [!Note]
> The build should also succeed on macOS but is not CI-tested.

```bash
python -m pip install -U pip setuptools wheel
python -m pip install vdist-solver-fortran
```

You can also install the development version directly from GitHub:

```bash
python -m pip install "git+https://github.com/Nkzono99/vdist-solver-fortran.git"
```

Pip builds use `INSTALL_PROFILE=auto` by default.  Override it with
`INSTALL_PROFILE=generic` or `INSTALL_PROFILE=camphor` when needed.

```bash
INSTALL_PROFILE=generic python -m pip install vdist-solver-fortran
```

## Quick start

```python
import emout
from vdsolverf.core import Particle
from vdsolverf.emses import get_backtrace

data = emout.Emout("EMSES-simulation-directory")

ts, probability, positions, velocities = get_backtrace(
    directory=data.directory,
    ispec=0,                             # 0 electron, 1 ion, 2 photoelectron
    istep=-1,
    particle=Particle([32, 32, 400], [0, 0, -10]),
    dt=data.inp.dt,
    max_step=300_000,
    output_interval=1,
    use_adaptive_dt=False,
    use_electric_field=True,             # False disables the electric field
    use_magnetic_field=True,             # False disables B, including background B
)
```

These keyword arguments also apply to probability evaluation, velocity-range
estimation, and the MPI wrappers. Both default to `True`. Output files for a
disabled field are not required.
Backward tracing reverses MPIEMSES3D's ordinary Boris update order. See
[Physics](docs/physics.en.md#boris-updates-and-backward-tracing) for magnetic
rotation and the conditions for reversing a step.

Photoelectron evaluation requires `use_raycast = .true.` in the EMSES
namelist `/emissn/` — see the
[Namelist reference](docs/namelist.en.md#raycast-photoelectron).

## Documentation

| Document | Summary |
|---|---|
| [Usage](docs/usage.en.md) | Python recipes: single and multi-particle backtrace, phase-space probability |
| [Velocity-space octree probability solver](docs/octree_probabilities.en.md) | Adaptive/sparse velocity-space search, algorithm, execution workflow, and sample/leaf output |
| [Velocity range autorange](docs/autorange.en.md) | Per-cell velocity range estimation, validation, and diagnostics |
| [Namelist reference](docs/namelist.en.md) | Supported `plasma.inp` groups and parameters |
| [Physics](docs/physics.en.md) | Liouville theorem, Maxwellian emission, raycast photoelectron |
| [Architecture](docs/architecture.en.md) | Fortran / Python layout, module graph, extension points |
| [Development](docs/development.en.md) | Build, test, contribute |

## Example notebook

- [Phase probability distribution solver & multiple backtraces](https://nbviewer.org/github/Nkzono99/examples/blob/main/examples/vdist-solver-fortran/example.ipynb)

## License

Apache License 2.0. See [LICENSE](LICENSE).

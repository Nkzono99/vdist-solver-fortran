# vdist-solver-fortran

> Lang: [日本語](README.md) | **English**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14018863.svg)](https://doi.org/10.5281/zenodo.14018863)

Velocity distribution solver for Python, implemented in Fortran.

The core computes phase-space probability densities by combining
backward-in-time particle tracing over an EMSES simulation with
boundary-tagged source distributions (Maxwellian, raycast photoelectron,
absorbing). The shared library is driven from Python through a `ctypes`
wrapper.

## Requirements

- `gfortran`
- Python 3.7+ (development uses 3.12 in `.venv/`)

## Install

Installation scripts are currently only guaranteed to work on Linux and
Windows.

> [!Note]
> The build should also succeed on macOS but is not CI-tested.

```bash
pip install git+https://github.com/Nkzono99/vdist-solver-fortran.git
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
)
```

Photoelectron evaluation requires `use_raycast = .true.` in the EMSES
namelist `/emissn/` — see the
[Namelist reference](docs/namelist.en.md#raycast-photoelectron).

## Documentation

| Document | Summary |
|---|---|
| [Usage](docs/usage.en.md) | Python recipes: single and multi-particle backtrace, phase-space probability, dust charging |
| [Namelist reference](docs/namelist.en.md) | Supported `plasma.inp` groups and parameters |
| [Physics](docs/physics.en.md) | Liouville theorem, Maxwellian emission, raycast photoelectron |
| [Architecture](docs/architecture.en.md) | Fortran / Python layout, module graph, extension points |
| [Development](docs/development.en.md) | Build, test, contribute |

## Example notebook

- [Phase probability distribution solver & multiple backtraces](https://nbviewer.org/github/Nkzono99/examples/blob/main/examples/vdist-solver-fortran/example.ipynb)

## License

Apache License 2.0. See [LICENSE](LICENSE).

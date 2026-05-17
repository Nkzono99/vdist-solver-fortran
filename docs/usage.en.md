# Usage

> Lang: [日本語](usage.md) | **English**

All examples assume you have an EMSES run on disk (e.g. `./my-run/`) that
`emout` can load, and that `vdist-solver-fortran` is installed in the active
Python environment.

> [!Note]
> `ispec` is a zero-based species index in the Python API.
> Conventional mapping: `0` = electron, `1` = ion, `2` = photoelectron.
> Photoelectron evaluation requires `use_raycast = .true.` in `/emissn/` —
> see [Namelist reference](namelist.en.md#raycast-photoelectron).

## Single backtrace

Trace one particle backward in time from a phase-space point and record the
probability of reaching a probability-tagged boundary.

```python
import emout
import matplotlib.pyplot as plt
from vdsolverf.core import Particle
from vdsolverf.emses import get_backtrace

data = emout.Emout("my-run")

particle = Particle(position=[32, 32, 400], velocity=[0, 0, -10])

ts, probability, positions, velocities = get_backtrace(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particle=particle,
    dt=data.inp.dt,
    max_step=300_000,
    output_interval=1,
    use_adaptive_dt=False,
)

plt.plot(positions[:, 0], positions[:, 2])
plt.gcf().savefig("backtrace.png")
```

## Multi-particle backtrace

Trace many particles in parallel (OpenMP) and overlay their paths weighted
by probability.

```python
import emout
import matplotlib.pyplot as plt
import numpy as np
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_backtraces

data = emout.Emout("my-run")

NVX, NVZ = 50, 50
phase_grid = PhaseGrid(
    x=32, y=32, z=130,
    vx=(-100, 100, NVX),
    vy=0,
    vz=(-400, -360, NVZ),
)

particles = phase_grid.create_particles()

ts, probabilities, positions, velocities, last_indexes = get_backtraces(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=10_000,
    output_interval=1,
    use_adaptive_dt=False,
    n_threads=4,
)

maxp = np.nanmax(probabilities)
for p, pos, li in zip(probabilities, positions, last_indexes):
    if np.isnan(p):
        continue
    alpha = min(1.0, p / maxp)
    plt.scatter(pos[:li, 0], pos[:li, 2], s=0.1, color="black", alpha=alpha)

plt.gcf().savefig("backtraces.png")
```

`last_indexes[i]` is the number of meaningful samples for particle `i`;
beyond that, positions are zero-padded.

## Phase-space probability solver

Skip trajectory storage and only return the probability and the last-step
phase-space point for each particle.

```python
import emout
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_probabilities

data = emout.Emout("my-run")

phase_grid = PhaseGrid(
    x=32, y=32, z=130,
    vx=(-100, 100, 50),
    vy=0,
    vz=(-400, -360, 50),
)

particles = phase_grid.create_particles()

probabilities, ret_particles = get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    use_adaptive_dt=False,
    n_threads=4,
)
```

## Selecting a specific shared-library path

The wrapper auto-detects OS and loads the bundled shared library, but you
can override:

```python
from vdsolverf.emses import get_backtrace

get_backtrace(..., system="linux", library_path="/custom/path/libvdist-solver-fortran.so")
```

Accepted values for `system`: `"auto"` (default), `"linux"`, `"darwin"`,
`"windows"`.

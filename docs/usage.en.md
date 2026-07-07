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

## Estimate per-cell velocity ranges before solving probabilities

You can forward-trace deterministic support particles from EMSES open-boundary
and emission-surface source distributions in Fortran, producing per-cell
velocity ranges. The returned `VelocityRangeMap` stores separate `vx/vy/vz`
ranges for each spatial cell, and `create_particles()` flattens those adaptive
ranges into a particle list accepted by `get_probabilities`.

```python
from vdsolverf.emses import estimate_velocity_range_map, get_probabilities

range_map = estimate_velocity_range_map(
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=0.25,
    max_step=30_000,
    use_adaptive_dt=True,
    coverage_sigma=4.0,
    safety_factor=1.25,
    collect_moments=False,
)

particles, index = range_map.create_particles(velocity_bins=(16, 8, 16))

probabilities, ret_particles = get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    use_adaptive_dt=False,
)

probability_grid = index.reshape(probabilities)
```

`range_map.count` is the number of envelope support-point hits. It is a
diagnostic sampling count, not a physical density. Cells with `count == 0`
are skipped by `create_particles()` by default. Pass `collect_moments=True`
when you also need the `mean_v` / `cov_v` diagnostics. See
[Per-cell velocity range autorange](autorange.en.md) for parameter choices
and validation.

## MPI particle parallelism

The existing `vdsolverf.emses.get_*` functions remain the OpenMP/threaded
entry points.  Use the optional MPI backend explicitly when you want
particle-parallel execution.

```python
from vdsolverf.emses.mpi import get_probabilities

probabilities, ret_particles = get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    n_threads=2,        # threads per rank
)
```

Use this form when the Python script itself is launched with MPI, for example
`srun -n 8 python script.py`.  `mpi4py` is an optional dependency; install it
only in MPI environments with `pip install "vdist-solver-fortran[mpi]"`.

To launch Slurm from an ordinary Python process, use the launcher wrapper:

```python
from vdsolverf.emses.mpi import srun_get_probabilities

probabilities, ret_particles = srun_get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    ntasks=8,
    n_threads=2,
    cpus_per_task=2,
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

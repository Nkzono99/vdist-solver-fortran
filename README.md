# vdist-solver-fortran

Velocity distribution solver for python written in Fortran.

## Requirements

- gfortran

## Install

Installation scripts are currently only guaranteed to work on Linux.

> [!Note]
> This script should work on Windows and Mac.
> 
> However, it has not been tested and is not guaranteed to work.

```
pip install git+https://github.com/Nkzono99/vdist-solver-fortran.git
```

## Example code

-  [Phase Probability Distribution Solver & Muliple Backtraces](https://nbviewer.org/github/Nkzono99/examples/blob/main/examples/vdist-solver-fortran/example.ipynb)

## Usage

### Backtrace

```python
import matplotlib.pyplot as plt
from vdsolver.core import Particle
from vdsolverf.emses import get_backtrace

position = [32, 32, 400]
velocity = [0, 0, -10]

particle = Particle(position, velocity)

ts, probabirity, positions, velocities = get_backtrace(
    directory="EMSES-simulation-directory",
    ispec=0,
    istep=-1,
    particle=particle,
    dt=0.002,
    max_step=300000,
    use_adaptive_dt=False,
)

print("Number of trace-points:", len(positions))

plt.plot(positions[:, 0], positions[:, 2])
plt.gcf().savefig("backtrace.png")
```

### Multiple Backtrace

```python
import matplotlib.pyplot as plt
import numpy as np
from vdsolver.core import Particle, PhaseGrid
from vdsolverf.emses import get_backtraces

NVX = 50
NVZ = 50
phase_grid = PhaseGrid(
    x=32,
    y=32,
    z=130,
    vx=(-100, 100, NVX),
    vy=0,
    vz=(-400, -360, NVZ),
)

phases = phase_grid.create_grid()
particles = []
for phase in phases.reshape(-1, phases.shape[-1]):
    pos = phase[:3].copy()
    vel = phase[3:].copy()
    pcl = Particle(pos, vel)
    particles.append(pcl)

ts, probabirities, positions, velocities = get_backtraces(
    directory="EMSES-simulation-directory",
    ispec=0,
    istep=-1,
    particles=particles,
    dt=0.002,
    max_step=10000,
    use_adaptive_dt=False,
    n_threads=4,
)

plt.figure(figsize=(15, 15))

data.phisp[-1, :, int(data.inp.ny//2), :].val_si.plot(use_si=False)

maxp = np.array(probabirities).max()
for probabirity, positions in zip(probabirities, positions_list):
    if np.isnan(probabirity):
        continue
    alpha = min(1.0, probabirity / maxp)

    plt.scatter(positions[:, 0], positions[:, 2], s=0.1, color='black', alpha=alpha)

plt.gcf().savefig("backtraces.png")
```

### Phase Probability Distribution Solver

```python
from vdsolver.core import Particle
from vdsolverf.emses import get_probabilities

particles = [
    Particle([16, 16, 400], [0, 0, -20]),
    Particle([16, 16, 400], [0, 0, -30]),
    ]

probabilities, ret_particles = get_probabilities(
        directory="EMSES-simulation-directory",
        ispec=0, # 0: electron, 1: ion, 2: photoelectron(not supported yet)
        istep=-1,
        particles=particles,
        dt=0.002,
        max_step=30000,
        adaptive_dt=False,
        n_threads=4,
    )

print(probabilities)
print(ret_particles)
```

```python
from vdsolver.core import Particle
from vdsolverf.emses import get_probabilities

NVX = 50
NVZ = 50
phase_grid = PhaseGrid(
    x=32,
    y=32,
    z=130,
    vx=(-100, 100, NVX),
    vy=0,
    vz=(-400, -360, NVZ),
)

phases = phase_grid.create_grid()
particles = []
for phase in phases.reshape(-1, phases.shape[-1]):
    pos = phase[:3].copy()
    vel = phase[3:].copy()
    pcl = Particle(pos, vel)
    particles.append(pcl)

probabilities, ret_particles = get_probabilities(
        directory="EMSES-simulation-directory",
        ispec=0, # 0: electron, 1: ion, 2: photoelectron(not supported yet)
        istep=-1,
        particles=particles,
        dt=0.002,
        max_step=30000,
        adaptive_dt=False,
        n_threads=4,
    )

phases = phases.reshape(NVZ, NVX, 6)
VX = phases[:, :, 3]
VZ = phases[:, :, 5]
probs = probs.reshape(NVZ, -1)
plt.pcolormesh(VX, VZ, probs, shading="auto")
plt.colorbar()
plt.show()
```

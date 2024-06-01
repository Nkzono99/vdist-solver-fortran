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

## Usage

### Backtrace

```python
import matplotlib.pyplot as plt
from vdsolver.core import Particle
from vdsolverf.emses import get_backtrace

position = [32, 32, 400]
velocity = [0, 0, -10]
max_step = 100000

particle = Particle(position, velocity)

ts, positions, velocities = get_backtrace(
    directory="EMSES-simulation-directory",
    ispec=0,
    istep=-1,
    particle=particle,
    dt=0.002,
    max_step=max_step,
    use_adaptive_dt=False,
)

print("Number of trace-points:", len(positions))

plt.plot(positions[:, 0], positions[:, 2])
plt.gcf().savefig("backtrace.png")
```

### Phase Probabirity Distribution Solver

```python
from vdsolver.core import Particle
from vdsolverf.emses import get_probabirities

particles = [
    Particle([16, 16, 400], [0, 0, -20]),
    Particle([16, 16, 400], [0, 0, -30]),
    ]

probabirities, ret_particles = get_probabirities(
        directory="EMSES-simulation-directory",
        ispec=0, # 0: electron, 1: ion, 2: photoelectron(not supported yet)
        istep=-1,
        particles=particles,
        dt=0.002,
        max_step=100000,
        adaptive_dt=False,
    )

print(probabirities)
print(ret_particles)
```

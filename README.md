# vdist-solver-fortran
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14018863.svg)](https://doi.org/10.5281/zenodo.14018863)

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

ts, probability, positions, velocities = get_backtrace(
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

ts, probabilities, positions, velocities = get_backtraces(
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

maxp = np.array(probabilities).max()
for probability, positions in zip(probabilities, positions_list):
    if np.isnan(probability):
        continue
    alpha = min(1.0, probability / maxp)

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

### NOTE: Placement of internal boundary objects (collision detection)

There are various methods to set internal boundaries for EMSES parameters, but this solver supports the following options. Please add these settings to the plasma.inp file.

First, specify the inner boundary type in the boundary_type or boundary_types(itype) field, and then configure the associated parameters.

#### Namelist Parameters
```
&ptcond
  boundary_type = 'none'|
                  'flat-surface'|
                  'rectangle-hole'|'cylinder-hole'|'hyperboloid-hole'|'ellipsoid-hole'|
                  'rectangle[xyz]'|'circle[x/y/z]'|'cuboid'|'disk[x/y/z]
                  'complex'

  ! Use if boundary_type is '****-surface' or '****-hole'.
  zssurf = Surface Height [grid]

  ! Use if boundary_type is '****-hole'.
  [x/y/z][l/u]pc = Hole [lower/upper] limit grid（[x/y/z] coordinate) [grid]

  ! Use if boundary_type is 'complex'
  boundary_types(ntypes) = <boundary_type>|

  ! Use if boundary_types(itype) is 'rectangle'
  rectangle_shape(ntypes, 6) = Rectangle location (xmin, xmax, ymin, ymax, zmin, zmax)

  ! Use if boundary_types is 'circle[x/y/z]'
  circle_origin(ntypes, 3) = Circle center coordinates
  circle_radius(ntypes) = Circle radius

  ! Use if boundary_types(itype) is 'cuboid'
  cuboid_shape(ntypes, 6) = Cuboid location (xmin, xmax, ymin, ymax, zmin, zmax)

  ! Use if boundary_types(itype) is 'disk[x/y/z]
  disk_origin(ntypes, 3) = Disk center bottom coordinates
  disk_height(ntypes) = Disk height (= thickness)
  disk_radius(ntypes) = Disk outer radius
  disk_inner_radius(ntypes) = Disk inner radius

  ! Rotation angle of all boundaries [deg].
  boundary_rotation_deg(3) = 0d0, 0d0, 0d0
&
```

### NOTE: Setting the particle emission surface

EMSES offers various settings for emission surfaces, and this solver supports the following options. Please add these settings to the plasma.inp file.

```
&emissn
    curf(nspec) = Current density corresponding to particle species (has not been supported)

    curfs(nepl) = Current density corresponding to the emitting surface (this has priority over curf) (has not been supported)

    nflag_emit(nspec) = For non-zero numbers, the emission surface setting works

    nepl(nspec) = Number of emission surfaces of each particle species

    nemd(nepl) = Normal direction of emission surface (-: - direction, +: + direction, 1: x direction, 2: y direction, 3: z direction)

    xmine(nepl), xmaxe(nepl), ymine(nepl), ymaxe(nepl), zmine(nepl), zmaxe(nepl)
        = Range of emission surfaces

    thetaz, thetaxy = Angle with magnetic field related to thermal velocity of emitted particles [deg]
&
```

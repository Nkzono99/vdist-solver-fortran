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
import emout
import matplotlib.pyplot as plt
from vdsolverf.core import Particle
from vdsolverf.emses import get_backtrace

data = emout.Emout("EMSES-simulation-directory")

position = [32, 32, 400]
velocity = [0, 0, -10]

ispec = 0 # 0: electron, 1: ion, 2: photoelectron
particle = Particle(position, velocity)

ts, probability, positions, velocities = get_backtrace(
    directory=data.directory,
    ispec=ispec,
    istep=-1,
    particle=particle,
    dt=data.inp.dt,
    max_step=300000,
    use_adaptive_dt=False,
)

print("Number of trace-points:", len(positions))

plt.plot(positions[:, 0], positions[:, 2])
plt.gcf().savefig("backtrace.png")
```

### Multiple Backtrace

```python
import emout
import matplotlib.pyplot as plt
import numpy as np
from vdsolverf.core import Particle, PhaseGrid
from vdsolverf.emses import get_backtraces

data = emout.Emout("EMSES-simulation-directory")

NVX = 50
NVZ = 50

ispec = 0 # 0: electron, 1: ion, 2: photoelectron
phase_grid = PhaseGrid(
    x=32,
    y=32,
    z=130,
    vx=(-100, 100, NVX),
    vy=0,
    vz=(-400, -360, NVZ),
)

phases = phase_grid.create_grid()
particles = phase_grid.create_particles()

ts, probabilities, positions_list, velocities_list, last_indexes = get_backtraces(
    directory=data.directory,
    ispec=ispec,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
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
import emout
from vdsolverf.core import Particle
from vdsolverf.emses import get_probabilities

data = emout.Emout("EMSES-simulation-directory")

ispec = 0 # 0: electron, 1: ion, 2: photoelectron
particles = [
    Particle([16, 16, 400], [0, 0, -20]),
    Particle([16, 16, 400], [0, 0, -30]),
    ]

probabilities, ret_particles = get_probabilities(
        directory=data.directory,
        ispec=ispec,
        istep=-1,
        particles=particles,
        dt=data.inp.dt,
        max_step=30000,
        adaptive_dt=False,
        n_threads=4,
    )

print(probabilities)
print(ret_particles)
```

```python
import emout
from vdsolver.core import Particle
from vdsolverf.emses import get_probabilities

data = emout.Emout("EMSES-simulation-directory")

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
particles = phase_grid.create_particles()

probabilities, ret_particles = get_probabilities(
        directory=data.directory,
        ispec=0, # 0: electron, 1: ion, 2: photoelectron(not supported yet)
        istep=-1,
        particles=particles,
        dt=data.inp.dt,
        max_step=30000,
        use_adaptive_dt=False,
        n_threads=4,
    )

phases = phases.reshape(NVZ, NVX, 6)
VX = phases[:, :, 3]
VZ = phases[:, :, 5]
probabilities = probabilities.reshape(NVZ, -1)
plt.pcolormesh(VX, VZ, probabilities, shading="auto")
plt.colorbar()
plt.show()
```

### NOTE: Placement of internal boundary objects (collision detection)

There are various methods to set internal boundaries for EMSES parameters, but this solver supports the following options. Please add these settings to the plasma.inp file.

First, specify the inner boundary type in the boundary_type or boundary_types(itype) field, and then configure the associated parameters.
Of course, it also supports the placement of objects via geotype.

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

  ! Geotype parameters
  npc = Number of geotype objects.
    geotype(npc) = Type of the object (0-1: cuboid, 2: cylinder, 3: sphere).

  ! For cuboid:
  xlpc(npc), xupc(npc), ylpc(npc), yupc(npc), zlpc(npc), zupc(npc) 
    = Boundaries of the cuboid in each dimension (x, y, z).

  ! For cylinder:
  bdyalign(npc) = Axis alignment of the cylinder (1: x-axis, 2: y-axis, 3: z-axis).
  bdyedge(1:2, npc) = Edge coordinates (lower and upper bounds).
  bdyradius(npc) = Radius of the cylinder.
  bdycoord(1:2, npc) = Center coordinates along the axis.

  ! For sphere:
  bdyradius(npc) = Radius of the sphere.
  bdycoord(1:3, npc) = Center coordinates of the sphere.
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

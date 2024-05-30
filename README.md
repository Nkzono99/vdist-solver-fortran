# vdist-solver-fortran

## Install

### Fortran Library

Installation scripts are currently only available for linux (using gfortran).

```
git clone https://github.com/Nkzono99/vdist-solver-fortran.git
source source_me.sh
vdist-solver-fortran.git
cd vdist-solver-fortran
make
```

### Python Wrapper
```
pip install git+https://github.com/Nkzono99/vdist-solver-fortran.git
```


## Usage

```python
from vdsolver.core import Particle
from vdsolverf import get_probabirities

particles = [
    Particle([16, 16, 400], [0, 0, -20]),
    Particle([16, 16, 400], [0, 0, -30]),
    ]

probabirities, ret_particles = get_probabirities(
        directory="EMSES-simulation-directory",
        ispec=0, # 0: electron, 1: ion, 2: photoelectron(not supported yet)
        istep=-1,
        particles=particles,
        dt_multiplier=1,
        max_step=100000,
    )

print(probabirities)
print(ret_particles)
```

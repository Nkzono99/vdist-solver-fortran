# Per-Cell Velocity Range Autorange

> Lang: [日本語](autorange.md) | **English**

`estimate_velocity_range_map` creates deterministic support particles from
EMSES source distributions, forward-traces them in Fortran, and estimates
separate `vx/vy/vz` ranges for each spatial cell. Instead of choosing one
global velocity box by hand, you can build a `VelocityRangeMap` and then
create the particle list passed to `get_probabilities`.

## Basic Workflow

```python
import emout
from vdsolverf.emses import estimate_velocity_range_map, get_probabilities

data = emout.Emout("my-run")

range_map = estimate_velocity_range_map(
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=0.25,
    max_step=30_000,
    use_adaptive_dt=True,
    coverage_sigma=4.0,
    safety_factor=1.25,
    n_threads=4,
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
    n_threads=4,
)

probability_grid = index.reshape(probabilities)
```

`probability_grid` has shape `(nz, ny, nx, nvz, nvy, nvx)`. Cells where no
range was estimated are skipped by `create_particles()` by default.

## Validation Pass

The support particles are not a strict nonlinear envelope. When you need a
more conservative range, run a coarse `get_probabilities` pass and expand only
cells whose velocity-box edge still contains significant probability.

```python
from vdsolverf.emses import (
    estimate_velocity_range_map,
    validate_and_expand_velocity_range_map,
)

range_map = estimate_velocity_range_map(
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=0.25,
    max_step=30_000,
    use_adaptive_dt=True,
    coverage_mode="relative_density",
    eps_rel=1e-6,
    safety_factor=1.25,
    n_threads=4,
)

range_map = validate_and_expand_velocity_range_map(
    range_map,
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=data.inp.dt,
    max_step=30_000,
    coarse_bins=(8, 4, 8),
    edge_threshold=1e-3,
    expand_factor=1.5,
    max_iter=2,
    n_threads=4,
)
```

Cells with `status == 3` were expanded by validation. Validation runs extra
`get_probabilities` calls, so keep `coarse_bins` much smaller than the final
production velocity grid.

## Main Parameters

| Parameter | Guideline |
|---|---|
| `dt` | Forward-trace step size. With `use_adaptive_dt=True`, start around `0.25` to `0.5` to limit per-step movement. |
| `max_step` | Maximum number of forward steps. Increase it when reflected particles need time to return. |
| `use_adaptive_dt` | Prefer `True` for range estimation to reduce skipped cells. |
| `coverage_sigma` | Source Maxwellian support radius. Use `3.0` for exploration, `4.0` as a default, and `5.0` for conservative runs. |
| `coverage_mode` / `eps_rel` | With `coverage_mode="relative_density"`, the code uses `sqrt(-2 log eps_rel)`. `eps_rel=1e-6` is about `5.26 sigma`. |
| `safety_factor` | Expands deposited min/max ranges about the center. Start with `1.25`; use `1.5` or more if tails are missed. |
| `source_samples_per_cell` | Tangential sampling count on each source surface cell. `1` is fastest; increase it when source-surface variation is strong. |
| `velocity_bins` | Velocity grid size per valid cell in `create_particles()`. The order is `(nvx, nvy, nvz)`. |
| `n_threads` | Number of OpenMP threads in Fortran. Defaults to `OMP_NUM_THREADS`, or `1` if unset. |
| `collect_moments` | Computes `mean_v` / `cov_v` when `True`. It uses more memory, so the default is `False`. |

## Return Values and Diagnostics

`estimate_velocity_range_map` returns a `VelocityRangeMap`. Most arrays have
shape `(nz, ny, nx)`.

| Attribute | Meaning |
|---|---|
| `vx_min`, `vx_max`, `vy_min`, `vy_max`, `vz_min`, `vz_max` | Per-cell velocity ranges. Cells with no hit contain `nan`. |
| `count` | Number of support-particle hits. This is a sampling diagnostic, not a physical density. |
| `weight_sum` | Sum of support-particle weights. Currently intended as a diagnostic. |
| `mean_v`, `cov_v` | Velocity-moment diagnostics, valid only when `collect_moments=True`. |
| `status` | `0`: OK, `1`: LOW_COUNT, `2`: FALLBACK/no hit, `3`: EDGE_EXPANDED. |
| `confidence` | Simple score, `min(count / minimum_count, 1)`. |
| `metadata` | Estimation settings such as `coverage_sigma` and `safety_factor`. |

`VelocityRangeMap.create_particles()` returns `(particles, index)`.
`particles` is a flattened particle list accepted by `get_probabilities`.
`index.reshape(probabilities)` maps flattened probabilities back to
`(nz, ny, nx, nvz, nvy, nvx)`.

## Accessing Velocity Ranges and 6D Distributions

`VelocityRangeMap` does not store the 6D probability distribution itself. It
stores the per-cell velocity ranges. The saved velocity-axis information is
the min/max pair for each spatial cell.

```python
range_map.vx_min      # shape: (nz, ny, nx)
range_map.vx_max
range_map.vy_min
range_map.vy_max
range_map.vz_min
range_map.vz_max
range_map.valid_mask  # shape: (nz, ny, nx)
```

Recover the velocity axes for one cell from the `velocity_bins` used for
visualization or particle creation.

```python
iz, iy, ix = 10, 20, 30
nvx, nvy, nvz = 16, 8, 8

vx = np.linspace(range_map.vx_min[iz, iy, ix], range_map.vx_max[iz, iy, ix], nvx)
vy = np.linspace(range_map.vy_min[iz, iy, ix], range_map.vy_max[iz, iy, ix], nvy)
vz = np.linspace(range_map.vz_min[iz, iy, ix], range_map.vz_max[iz, iy, ix], nvz)
```

Use the `index` returned by `create_particles()` to view `get_probabilities`
results as a 6D array.

```python
velocity_bins = (16, 8, 8)  # (nvx, nvy, nvz)
particles, index = range_map.create_particles(velocity_bins=velocity_bins)

probabilities, _ = get_probabilities(..., particles=particles)
prob_grid = index.reshape(probabilities)

print(prob_grid.shape)
# (nz, ny, nx, nvz, nvy, nvx)

cell_probability = prob_grid[iz, iy, ix]
# shape: (nvz, nvy, nvx)
```

`range_map.save()` stores the range arrays such as `vx_min/vx_max` plus
diagnostics. It does not store `prob_grid` or every `vx/vy/vz` grid point.
After loading, pass the same `velocity_bins` to reconstruct the same axes.

## Cell-Level Access

Use `range_map[iz, iy, ix]` to inspect or sample a single spatial cell. It
returns a `VelocityRangeCell` view.

```python
cell = range_map[iz, iy, ix]

print(cell.position)     # cell center [x, y, z]
print(cell.vmin)         # [vx_min, vy_min, vz_min]
print(cell.vmax)         # [vx_max, vy_max, vz_max]
print(cell.count)
print(cell.status)

vx, vy, vz = cell.velocity_axes((16, 8, 8))
particles, index = cell.create_particles(velocity_bins=(16, 8, 8))
probabilities, _ = get_probabilities(..., particles=particles)

cell_probability = index.reshape(probabilities)
# shape: (nvz, nvy, nvx)
```

For invalid cells, `cell.valid == False`, and `cell.create_particles()` returns
an empty particle list by default.

## Saving and Loading

A `VelocityRangeMap` returned by `estimate_velocity_range_map` remembers the
source `data.directory`, `ispec`, and `istep`. Therefore, `save()` without
arguments writes to the default file name:

```python
path = range_map.save()
print(path)
# data.directory / "vdsolverf-velocity-range-map-ispec0-istep-1.npz"
```

Load from the same default path with:

```python
from vdsolverf.core import VelocityRangeMap

range_map = VelocityRangeMap.load(
    directory=data.directory,
    ispec=0,
    istep=-1,
)
```

Use `range_map.save("path/to/range-map.npz")` for an explicit path, and
`VelocityRangeMap.load("path/to/range-map.npz")` to load from one.

## Reflection and Scope

Forward tracing uses the actual EM fields and boundaries, so reflection caused
by the potential distribution or by reflecting boundary conditions is included
if those trajectories pass through the cells. Typical misses come from too
small `max_step`, too large `dt`, or support particles too sparse to represent
a reflected velocity lobe.

The current source-envelope implementation covers EMSES open boundaries and
emission surfaces. Raycast photoelectrons and secondary electrons generated at
inner boundaries are not yet automatic source-envelope inputs. For those
components, combine this feature with manual velocity ranges or validation.

## Performance Notes

The estimator runs in Fortran and parallelizes over source patches with OpenMP.
Deposits write to thread-local accumulators and merge at the end; it does not
use `atomic` deposits. This avoids atomic contention, but memory usage scales
roughly with `number of cells * number of threads`. `collect_moments=True`
adds more thread-local moment arrays, so start with `False`.

The included microbenchmark can be built manually:

```bash
mpiifort -O3 -qopenmp -Iinclude benchmarks/bench_emses_autorange.f90 \
  -Llib -lvdist-solver-fortran -Wl,-rpath,"$PWD/lib" \
  -o /tmp/bench_emses_autorange

/tmp/bench_emses_autorange 4 32 32 32 96 5
```

Arguments are `n_threads lx ly lz max_step reps`. On KUDPC login nodes, run the
benchmark through `tssrun` or inside a batch job, not directly on the login node.

# Per-Cell Velocity Range Autorange

> Lang: [日本語](autorange.md) | **English**

`estimate_velocity_range_map` creates deterministic support particles from
EMSES source distributions, forward-traces them in Fortran, and estimates
separate `vx/vy/vz` ranges for each spatial cell. Instead of choosing one
global velocity box by hand, you can build a `VelocityRangeMap` and then
create the particle list passed to `get_probabilities`.

External-source support widths use Cartesian marginal standard deviations of
the rotated thermal covariance. Each trajectory segment is interpolated to
its unwrapped endpoint, then each deposited point is mapped into the periodic
domain. This applies to several periodic axes and multiple crossings within
one step.

## Field bundles and standalone tier-2 batches

Both `estimate_velocity_range_map(..., field_bundle=...)` and
`validate_and_expand_velocity_range_map(..., field_bundle=...)` accept a bundle
directory or `FieldBundle`, just like
[`get_dust_currents`](dust_currents.en.md). They reuse the relocated fields and
trimmed namelist without opening EMSES output files. The requested `directory`,
`istep`, species and `use_electric_field` / `use_magnetic_field` must match the
bundle. With the same settings and thread count, bundled and ordinary inputs
produce the same ranges. Validation also uses the bundle for its coarse
probability evaluations.

```python
from vdsolverf.io import write_field_bundle, write_range_batch, read_range_batch_results
from vdsolverf.emses import estimate_velocity_range_map, validate_and_expand_velocity_range_map

bundle = write_field_bundle("my-run", -1, [0, 1], "my-run/field_bundle")
range_map = estimate_velocity_range_map(
    "my-run", 0, -1, dt=0.25, max_step=30_000,
    n_threads=8, field_bundle=bundle,
)
range_map = validate_and_expand_velocity_range_map(
    range_map, "my-run", 0, -1, dt=0.1, max_step=30_000,
    n_threads=8, field_bundle=bundle,
)

# Prepare once in Python; the compute job below needs only Fortran.
batch_nml = write_range_batch(
    "range-batch", bundle, 0, dt=0.25, max_step=30_000,
    use_adaptive_dt=True, coverage_mode="relative_density", eps_rel=1e-6,
    safety_factor=1.25, source_samples_per_cell=1, minimum_count=4,
    accumulator_cache_size=20_000, n_threads=8,
)
```

```bash
# Build/run on compute nodes at KUDPC; -t always uses H:MM:SS.
tssrun -p gr20001b -t 0:40:00 --rsc p=1:t=8:c=8 make batch
tssrun -p gr20001b -t 0:40:00 --rsc p=1:t=8:c=8 \
  build/batch-serial/bin/vdsolverf-batch range-batch/batch.nml
```

```python
range_map = read_range_batch_results("range-batch")
range_map.save()  # Original directory, species and istep are restored.
# Feed each selected cell's six bounds into write_batch for current integration.
```

`&vdsolverf_batch mode='velocity_range'` requires only `batch.nml`; it never
reads `targets.npy` or `psi.npy`. `mode='dust_currents'` remains the default.
Range keys are `bundle_dir`, zero-based `ispec`, `dt`, `max_step`,
`use_adaptive_dt`, `coverage_sigma`, `safety_factor`, `source_samples_per_cell`,
`minimum_count`, `accumulator_cache_size`, `n_threads`, `output_dir`, and optional
`max_probability_types`. The writer resolves `coverage_mode` / `eps_rel` to
`coverage_sigma`. Defaults match the Python range API (`dt=0.25`,
`max_step=1000`); set the trace duration explicitly for production.
Omit `n_threads` to use the execution environment's `OMP_NUM_THREADS`.
Field flags are inherited from the bundle. The mode runs one OpenMP job;
an MPI executable computes on rank 0 only and the other ranks exit cleanly.

Under `output_dir/velocity_range/` (default `out/velocity_range/`), the eight
`vx_min/vx_max/vy_min/vy_max/vz_min/vz_max/weight_sum/confidence.npy` files are
float64, and `count.npy` / `status.npy` are int32. All have C-order shape
`(nz, ny, nx)`. `summary.json` in the same directory contains status counts,
settings, wall time and an input fingerprint, and is atomically written last as
the completion marker. Moments are omitted; the reader supplies NaN views.

Resume is enabled by default (also `--resume`): an unchanged completed job
leaves all output files untouched. The executable streams a CRC64 fingerprint
over `batch.nml`, `bundle.json`, the trimmed namelist and the selected species'
field npy, and includes the effective thread count. This reads the entire field
file even when skipping. Any input-byte change or missing output recomputes the
map; `--no-resume` always recomputes it. Keep input files immutable during a run.
See [dust-current batches](dust_currents.en.md#standalone-batch) for builds,
MPI submission and the subsequent current-integration stage.

## Reproducibility

Support particles and trajectories are deterministic. Each worker owns its
cache; global min/max, counts, weights and velocity moments are merged under
cell locks. Cache flush order can change between runs, including runs with
the same thread count. The optional `mean_v` and `cov_v` therefore have
floating-point summation differences; they are not bit-identical outputs.
On the 3 x 2 x 4 nonuniform-field fixture, repeated 1/2/8-thread ifx runs
differed by less than `1e-15` in these normalized diagnostics, while bounds,
counts, weights, status and confidence matched bit for bit. For this fixture,
compare moments with `rtol=1e-12, atol=1e-12`; longer or differently scaled
traces require an appropriate absolute tolerance. Use one thread when exact
diagnostic reproducibility is required.

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

Note: `estimate_velocity_range_map` defaults to `use_adaptive_dt=False`.
Commenting out that argument uses the ordinary time-step interpretation. When
you want to reduce skipped cells, explicitly pass `use_adaptive_dt=True`. In
that mode, `dt` is an approximate grid-distance limit per step, not
`data.inp.dt` from the EMSES run. For example, `dt=0.004` with
`max_step=30000` reaches only about 120 grid cells from the source in a
straight path.

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
| `dt` | Forward-trace step size. With `use_adaptive_dt=True`, this is a grid-distance limit per step, not the EMSES `data.inp.dt`. Start around `0.25` to `0.5`. |
| `max_step` | Maximum number of forward steps. Increase it when reflected particles need time to return. |
| `use_adaptive_dt` | Defaults to `False`. Pass `True` explicitly when you want to reduce skipped cells. |
| `coverage_sigma` | Source Maxwellian support radius. Use `3.0` for exploration, `4.0` as a default, and `5.0` for conservative runs. |
| `coverage_mode` / `eps_rel` | With `coverage_mode="relative_density"`, the code uses `sqrt(-2 log eps_rel)`. `eps_rel=1e-6` is about `5.26 sigma`. |
| `safety_factor` | Expands deposited min/max ranges about the center. Start with `1.25`; use `1.5` or more if tails are missed. |
| `source_samples_per_cell` | Tangential sampling count on each source surface cell. `1` is fastest; increase it when source-surface variation is strong. |
| `velocity_bins` | Velocity grid size per valid cell in `create_particles()`. The order is `(nvx, nvy, nvz)`. |
| `n_threads` | Number of OpenMP threads in Fortran. Defaults to `OMP_NUM_THREADS`, or `1` if unset. |
| `collect_moments` | Computes `mean_v` / `cov_v` when `True`. It uses more memory, so the default is `False`. When `False`, the wrapper does not allocate dense moment arrays and `range_map.mean_v` / `cov_v` are NaN views. |
| `show_progress` | Shows a progress bar in the same style as `get_probabilities` when `True`. Set `False` to keep batch logs quiet. |
| `accumulator_cache_size` | Number of hit cells cached by each OpenMP thread before flushing. The default is `20000`. Larger values reduce flush synchronization but use more memory. |

## Memory and Parallelism

`estimate_velocity_range_map` avoids deposit-time atomics by giving each thread
a bounded sparse cache. Once a cache reaches `accumulator_cache_size`, it is
reduced into the global range map using striped cell locks.

Memory use is roughly:

```text
global dense range arrays
+ EB field
+ n_threads * accumulator_cache_size * cache_entry_size
```

The implementation no longer allocates dense thread-local arrays with shape
`all spatial cells * n_threads`. For large grids with 112 threads, start with
the default `accumulator_cache_size=20000` and increase it only if flush
synchronization dominates.

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
Deposits use bounded per-thread sparse caches and flush under cell locks, rather
than per-deposit atomics. Memory scales with the dense map plus
`n_threads * accumulator_cache_size` cache entries. `collect_moments=True`
adds global moment arrays and cache moment storage, so start with `False`.

The included microbenchmark can be built manually:

```bash
mpiifort -O3 -qopenmp -Iinclude benchmarks/bench_emses_autorange.f90 \
  -Llib -lvdist-solver-fortran -Wl,-rpath,"$PWD/lib" \
  -o /tmp/bench_emses_autorange

/tmp/bench_emses_autorange 4 32 32 32 96 5
```

Arguments are `n_threads lx ly lz max_step reps`. On KUDPC login nodes, run the
benchmark through `tssrun` or inside a batch job, not directly on the login node.

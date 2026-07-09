# Octree Velocity-Space Probability Solver

> Lang: [日本語](octree_probabilities.md) | **English**

`get_probabilities_octree` adaptively subdivides velocity space with an octree
for each requested spatial point and evaluates `get_probabilities`-equivalent
probabilities in Fortran. Use it when a dense uniform 6D array is too expensive
and the probability structure is localized or multi-lobed in velocity space.

## Basic Example

```python
import emout
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_probabilities_octree

data = emout.Emout("my-run")

phase_grid = PhaseGrid(
    x=(120, 180, 16),
    y=64,
    z=(300, 420, 16),
    vx=(-8.0e6, 8.0e6, 2),
    vy=(-4.0e6, 4.0e6, 2),
    vz=(-8.0e6, 8.0e6, 2),
)

result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    phase_grid=phase_grid,
    dt=data.inp.dt,
    max_step=30_000,
    use_adaptive_dt=False,
    scout_bins=(11, 9, 11),
    max_depth=5,
    max_samples_per_cell=80_000,
    max_leaves_per_cell=4096,
    n_threads=32,
)
```

`PhaseGrid` `x/y/z` defines the spatial points. `vx/vy/vz` define the root
velocity box for every spatial point. The velocity bin counts in `PhaseGrid` are
not used by the octree, only the velocity limits are used.

Three-element `PhaseGrid` tuples mean `(start, end, count)`. For example,
`x=(120, 180, 16)` creates 16 points between 120 and 180 with `np.linspace`;
it does not mean a stride of 16.

The `max_depth=5`, `max_samples_per_cell=80_000`, and
`max_leaves_per_cell=4096` values above are a budget-limited practical example,
not enough capacity for a full worst-case depth-5 tree. Increase capacities or
split the velocity box manually when many points return status `2` or `4`.

## Arbitrary Points And Per-Point Bounds

You can pass explicit spatial points and per-point velocity bounds instead of a
`PhaseGrid`.

```python
positions = [
    [120.5, 64.5, 320.5],
    [121.5, 64.5, 320.5],
]

velocity_bounds = [
    [-8e6, 8e6, -4e6, 4e6, -8e6, 8e6],
    [-6e6, 6e6, -3e6, 3e6, -9e6, 7e6],
]

result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    position=positions,
    velocity_bounds=velocity_bounds,
    dt=data.inp.dt,
    max_step=30_000,
)
```

`velocity_bounds` accepts `(6,)`, `(3, 2)`, `(nspatial, 6)`, or
`(nspatial, 3, 2)`. A single `(6,)` or `(3, 2)` bound is broadcast to all
spatial points.

## Result

The return value is `VelocityOctreeResult`. It stores compact sample and octree
box arrays, not a dense 6D grid.

| Attribute | Shape | Meaning |
|---|---:|---|
| `spatial_points` | `(nspatial, 3)` | Evaluated spatial points |
| `velocities` | `(nsample, 3)` | Velocity samples |
| `probabilities` | `(nsample,)` | Probability at each sample; misses are `nan` |
| `spatial_index` | `(nsample,)` | Spatial point index for each sample |
| `leaf_bounds` | `(nleaf, 6)` | Velocity bounds of evaluated octree boxes, including split intermediate boxes |
| `leaf_value_min/max` | `(nleaf,)` | Min/max sampled probability in each box |
| `leaf_depth` | `(nleaf,)` | Octree depth |
| `leaf_sample_start/count` | `(nleaf,)` | Sample span for each box in the compact `velocities/probabilities` arrays |
| `status` | `(nspatial,)` | Per-spatial status |
| `sample_count`, `leaf_count` | `(nspatial,)` | Per-spatial output counts |
| `metadata` | `dict` | Runtime parameters, `actual_sample_count`, and `actual_leaf_count` |

Status values:

| Value | Meaning |
|---:|---|
| `0` | OK |
| `1` | No positive probability signal was found in the root box |
| `2` | `max_samples_per_cell` was reached |
| `3` | The root box was expanded at least once. Edge signal may or may not remain |
| `4` | `max_leaves_per_cell` was reached |
| `5` | Invalid velocity bounds |

Nonzero status points can still contain scout or partial samples. Check status
before visualization or downstream reductions.

```python
import numpy as np

bad = np.flatnonzero(result.status != 0)
if bad.size:
    print("non-OK octree points:", bad[:10], "status:", result.status[bad[:10]])
```

## Visualization

Select samples for one spatial point with `spatial_index`.

```python
import numpy as np

i = 0
if result.status[i] != 0:
    print("warning: non-OK octree status", result.status[i])

v = result.velocities_for_spatial(i)
p = result.probabilities_for_spatial(i)

valid = np.isfinite(p)
ax.scatter(v[valid, 0], v[valid, 2], c=p[valid], s=2)
```

For octree box projections, use `leaf_spatial_index` and `leaf_bounds`.
The `leaf_*` name is kept for output compatibility; these arrays include
evaluated intermediate boxes as well as terminal leaves.

```python
leaf_mask = result.leaf_mask(i)
boxes = result.leaf_bounds[leaf_mask]
value_range = result.leaf_value_max[leaf_mask] - result.leaf_value_min[leaf_mask]
```

To select samples that belong to one leaf, slice the compact sample arrays.

```python
leaf_index = np.flatnonzero(leaf_mask)[0]
sample_slice = result.leaf_sample_slice(leaf_index)
leaf_velocities = result.velocities[sample_slice]
leaf_probabilities = result.probabilities[sample_slice]
```

## Parameters

| Argument | Guideline |
|---|---|
| `scout_bins` | Root-box scout samples. Increase this for narrow separated lobes. |
| `max_depth` | Maximum octree subdivision depth. Higher values increase resolution and cost rapidly. |
| `refine_threshold_rel` | Split a box when `pmax - pmin` relative to `pmax` exceeds this value. |
| `edge_threshold_rel` | Threshold for detecting significant probability on the root-box boundary. |
| `expand_factor` / `max_expansions` | Expand root boxes that still have edge signal. |
| `max_samples_per_cell` | Per-spatial sample capacity. Increase when status `2` appears. |
| `max_leaves_per_cell` | Per-spatial octree-box capacity. Increase when status `4` appears. |
| `n_threads` | Fortran OpenMP thread count. Parallelism is over spatial points. |

If every box is split to depth `D`, the worst-case number of evaluated boxes is
`(8 ** (D + 1) - 1) / 7`. Non-root boxes are sampled with `3x3x3` points, so
set output capacities with margin.

## Notes

- This API evaluates probabilities directly, so it does not propagate the
  envelope approximation errors of `estimate_velocity_range_map`.
- If the initial `scout_bins` completely misses a narrow lobe, that lobe will
  not trigger refinement. For multiple narrow Maxwellian lobes, increase
  `scout_bins` or split the velocity box manually and evaluate multiple calls.
- `leaf_bounds` includes evaluated intermediate boxes as well as terminal boxes.
- The result is a sparse sample/box representation. To save it, use `np.savez`
  with arrays such as `result.velocities`, `result.probabilities`, and
  `result.leaf_bounds`.

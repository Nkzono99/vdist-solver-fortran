# Fortran Octree Probability Evaluation Design

## Purpose

Add a high-speed adaptive velocity-space evaluator for EMSES probability
queries. The feature targets case-study workloads where users evaluate
`get_probabilities` on roughly 1k to 10k spatial cells and want compressed
diagnostics, not stored 6D probability arrays.

The new API must preserve existing public APIs:

- `get_backtrace`
- `get_backtraces`
- `get_probabilities`
- `estimate_velocity_range_map`

The new functionality is added as `get_probabilities_octree`.

## Core Direction

The heavy work is implemented in Fortran. Python provides argument validation,
ctypes binding, result objects, and optional reducer orchestration.

The evaluator accepts PhaseGrid-like 6D bounds, but does not build a dense 6D
grid. It treats the spatial axes as a list of spatial sample points and builds
an independent 3D velocity octree for each spatial point.

```text
spatial point 0 -> velocity octree
spatial point 1 -> velocity octree
...
```

This keeps the adaptive logic aligned with the actual target: a velocity
distribution per spatial cell or sample point.

## Public API

Python entry point:

```python
result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    phase_grid=PhaseGrid(
        x=(100, 300, 50),
        y=100.0,
        z=(200, 800, 100),
        vx=(-8e6, 8e6, 2),
        vy=(-8e6, 8e6, 2),
        vz=(-8e6, 8e6, 2),
    ),
    scout_bins=(9, 9, 9),
    max_depth=6,
    max_samples_per_cell=50_000,
    refine_threshold_rel=1e-4,
    edge_threshold_rel=1e-4,
    expand_factor=1.5,
    max_expansions=3,
    dt=data.inp.dt,
    max_step=30_000,
    use_adaptive_dt=False,
    n_threads=112,
)
```

Alternative for one spatial point:

```python
result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    position=[x, y, z],
    velocity_bounds=((-8e6, 8e6), (-8e6, 8e6), (-8e6, 8e6)),
    ...
)
```

The two forms are mutually exclusive. `phase_grid` is the multi-point form,
and `position` plus `velocity_bounds` is the single-point form.

## Result Model

The result is ragged, not dense 6D:

```python
VelocityOctreeResult
    spatial_points          # shape: (nspatial, 3)
    velocities              # shape: (nsample, 3)
    probabilities           # shape: (nsample,)
    spatial_index           # shape: (nsample,)
    leaf_bounds             # shape: (nleaf, 6), vx0,vx1,vy0,vy1,vz0,vz1
    leaf_value_min          # shape: (nleaf,)
    leaf_value_max          # shape: (nleaf,)
    leaf_sample_start       # shape: (nleaf,)
    leaf_sample_count       # shape: (nleaf,)
    status                  # shape: (nspatial,)
    metadata
```

The MVP exposes ragged samples and leaf diagnostics. Follow-up helpers can add
cell-wise reduction and interpolation:

```python
result.reduce_by_cell(reducer)
result.to_regular_grid(spatial_index, velocity_bins=(50, 50, 50))
```

The MVP does not save dense 6D data. If the caller needs 50x50x50 data for a
cell, it is reconstructed or evaluated for that cell only.

## Fortran Entry Point

Add a C ABI entry point in `src/emses/emses_solver.f90`:

```fortran
subroutine get_probabilities_octree(...) bind(c)
```

The implementation lives in a new module:

```text
src/emses/emses_octree_probabilities.f90
```

The Fortran entry point receives:

- EMSES input path and relocated E/B arrays, matching existing wrappers.
- Spatial points as `(3, nspatial)`.
- Initial velocity bounds as `(6, nspatial)` or a single shared `(6)`.
- Octree controls.
- Output buffers sized by Python from conservative capacity estimates.
- Counters returning actual sample and leaf counts.

The first implementation uses capacity-limited flat arrays instead of pointer
graphs. This is faster, easier to bind to Python, and more predictable on HPC.

## Octree Algorithm

For each spatial point:

1. Create one root velocity box from the initial bounds.
2. Evaluate a deterministic scout grid inside the root box.
3. If probability is significant on any root boundary, expand the root bounds
   and repeat up to `max_expansions`.
4. Push the root box into a breadth-first queue.
5. Sample each queued box. The root uses `scout_bins`; child boxes currently use
   `3x3x3` samples.
6. Split a box when:
   - `depth < max_depth`
   - sample and leaf budgets remain
   - `pmax > 0`
   - `(pmax - pmin) / pmax >= refine_threshold_rel`
7. Stop when no box requests refinement, `max_samples_per_cell` is reached, or
   `max_leaves_per_cell` is reached.

The splitter is 3D octree: every refined velocity box is split into 8 children.

## Multiple Maxwellian Lobes

The algorithm must not assume a single connected distribution. It supports
multiple lobes through:

- A scout grid over the entire velocity box.
- Independent queued boxes.
- Edge expansion when any lobe reaches a velocity boundary.
- Refinement based on local box scores, not only global center values.

Limitations are explicit: a lobe narrower than the scout spacing can be missed.
Users control this with `scout_bins`, `max_depth`, and
`max_samples_per_cell`. A deterministic low-discrepancy exploration mode can be
added after the MVP if narrow disconnected lobes are still missed.

## Performance Design

The Fortran implementation should evaluate particles in batches, not one leaf
at a time through Python. The solver construction and EMSES field loading happen
once per call.

Parallelism:

- OpenMP parallelizes over spatial points first.
- If `nspatial` is small, OpenMP parallelizes over leaf/sample batches.
- Each thread uses local append buffers for samples and leaves.
- Thread buffers are reduced into global arrays in chunks to avoid per-sample
  atomics.

Memory:

- Output capacity is controlled by:
  - `max_samples_per_cell * nspatial`
  - `max_leaves_per_cell * nspatial`
- If capacity is exceeded, the spatial point is marked
  `OUTPUT_CAPACITY_EXCEEDED` and the current best leaves are returned.
- Returned arrays are sliced to actual counts in Python.

The first implementation should prioritize fast probability evaluation and
predictable memory over a sophisticated tree data structure.

## Reducer Integration

Reducer streaming is a follow-up layer. The octree evaluator first returns
ragged samples and leaves for selected cells. Then a second API can stream
reduced output:

```python
reduced = compute_reduced_probabilities_octree(
    ...,
    reducer=reducer,
    reducer_output_shape=(k,),
    max_batch_bytes="2GB",
)
```

This separation keeps the Fortran octree evaluator testable and avoids baking
Python callback behavior into the Fortran hot path.

## Error Handling and Status

Per-spatial-point status values:

- `0`: OK
- `1`: NO_SIGNAL
- `2`: SAMPLE_LIMIT_REACHED
- `3`: ROOT_EXPANDED
- `4`: OUTPUT_CAPACITY_EXCEEDED
- `5`: INVALID_BOUNDS

Fortran returns status arrays and actual counts. Python raises only for invalid
global arguments or ABI failures. Per-cell numerical issues are represented in
the result status.

## Testing

Python tests:

- API argument validation for mutually exclusive `phase_grid` and
  `position`/`velocity_bounds`.
- ctypes signature test with a fake DLL.
- result slicing from overallocated Fortran buffers.
- multi-lobe toy result handling.

Fortran tests:

- Octree refinement on a synthetic Gaussian probability function.
- Two separated Gaussian lobes are both retained in evaluated octree boxes.
- Edge expansion triggers when a lobe touches the root velocity boundary.
- Sample limit returns `STATUS_SAMPLE_LIMIT_REACHED`.
- Threaded result matches serial result for deterministic inputs.

Integration smoke:

- Small EMSES run with one spatial point and low sample budget.
- Real-data short run through `tssrun` to verify progress, memory, and
  OpenMP behavior.

## Documentation

Add:

```text
docs/octree_probabilities.md
docs/octree_probabilities.en.md
```

The docs must state that the method is adaptive and does not guarantee finding
features narrower than the scout grid. They should recommend starting with
small `nspatial`, inspecting leaf diagnostics, and then scaling to 1k to 10k
cells.

## Implementation Order

1. Add Python result dataclasses and public API skeleton.
2. Add fake-DLL Python tests for ABI and argument behavior.
3. Add Fortran synthetic octree module tests.
4. Implement Fortran flat-array octree evaluator for synthetic probability
   functions.
5. Wire the evaluator to EMSES `t_Solver%calculate_probability`.
6. Add ctypes wrapper and public export.
7. Add docs and small smoke examples.
8. Run Python tests, focused Fortran tests, full Fortran tests, build shared
   library, and real-data smoke through `tssrun`.

## Non-Goals

- Do not replace `get_probabilities`.
- Do not store dense 6D arrays.
- Do not implement full 6D octree over spatial and velocity axes in the MVP.
- Do not implement Python callbacks inside Fortran.
- Do not solve cell-to-cell autorange transport in this feature.

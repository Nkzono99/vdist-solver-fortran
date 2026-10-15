# Dust charging current moments

> Lang: [日本語](dust_currents.md) | **English**

`vdsolverf.emses.get_dust_currents` evaluates the OML velocity integral inside
Fortran, using one adaptive velocity octree per spatial point for every potential
in `psi_grid`. It returns moments, without exporting the velocity samples.

## Definition and units

All inputs and outputs use **EMSES normalised units**, as in
[`get_probabilities_octree`](octree_probabilities.en.md). `ispec` is zero-based in
Python and converted to one-based at the C boundary. For the selected species,
`qm` is read from the namelist:

```text
F(psi) = integral p(v) |v| max(0, 1 - 2 qm psi / |v|^2) d^3v
j_s(psi) = q_s n_s F(psi)
I_s(psi) = pi a^2 j_s(psi)
```

`p(v)` is the solver's density-one boundary-source PDF, transported to the
requested point by backtracing. Density, species charge, geometric area and SI
conversion are the caller's responsibility. `psi` is the dust potential relative
to the local plasma potential, in the normalisation compatible with `qm` and
velocity. Velocities are in the stationary dust frame; no dust-velocity shift is
applied. The calculation includes neither trapped populations disconnected from
a source boundary nor feedback from the dust onto the plasma. OML assumes a dust
radius much smaller than the local Debye length and negligible background
potential variation across its collection region. Near-surface collection and
emission require their own physical model.

For an isotropic, non-drifting Maxwellian, `path = peth = vth_1d = sqrt(kT/m)`
(the one-dimensional standard deviation, **not** `sqrt(2 kT/m)`):

```text
chi = qm psi / vth_1d^2
<v> = sqrt(8/pi) vth_1d
F = <v> exp(-chi)    for chi > 0
F = <v> (1 - chi)    for chi <= 0
```

## Usage and results

```python
from vdsolverf.emses import get_dust_currents
from vdsolverf.io import write_field_bundle

bundle = write_field_bundle("my-run", -1, [0, 1], "my-run/field_bundle")
result = get_dust_currents(
    "my-run", ispec=0, istep=-1,
    positions=[[12.5, 8.5, 4.5], [13.5, 8.5, 4.5]],
    velocity_bounds=[-5., 5., -5., 5., -5., 5.],
    psi_grid=[-1., 0., 1.], dt=0.1, max_step=30_000,
    n_threads=8, field_bundle=bundle,
)
# Current density in matching normalised units:
j = species_charge * species_density * result.flux
print(result.status, result.sample_count, result.terminal_box_count)
```

`velocity_bounds` accepts `(6,)`, `(3, 2)`, `(nspatial, 6)` or `(nspatial, 3, 2)`;
a single box broadcasts to all points. `positions` accepts `(3,)` or
`(nspatial, 3)`. `psi_grid` is a nonempty finite one-dimensional array; its
input order is preserved. `DustCurrentResult` is exported from `vdsolverf.core`.

| Attribute | Python shape | Meaning |
|---|---|---|
| `positions` | `(nspatial, 3)` | Requested positions |
| `psi_grid` | `(npsi,)` | Requested potentials |
| `flux`, `flux_error` | `(nspatial, npsi)` | Density-one flux moment and sampled range estimate |
| `status` | `(nspatial,)` | Octree status |
| `sample_count` | `(nspatial,)` | Stored backtrace samples, including intermediate boxes |
| `terminal_box_count` | `(nspatial,)` | Integrated boxes with no children |
| `metadata` | `dict` | Input parameters, units and integration rule |

Field flags, platform/library overrides and temporary namelist suffix follow
the existing octree API. Octree defaults are identical: `scout_bins=(9,9,9)`,
`max_depth=6`, `max_samples_per_cell=50000`, `max_leaves_per_cell=4096`,
`refine_threshold_rel=edge_threshold_rel=1e-4`, `expand_factor=1.5`,
`max_expansions=3`. `max_probability_types=100`, `use_adaptive_dt=False`;
`n_threads=None` follows `OMP_NUM_THREADS`, defaulting to one.

## Integration and error estimate

The root scout grid includes the box faces (`linspace`); child boxes use a
`3 x 3 x 3` endpoint lattice, with vx varying fastest. Each **terminal box's own
samples** receive tensor-product trapezoid weights: along each axis the endpoint
weights are `1/2`, interior weights are `1`, multiplied by the grid spacing.
Intermediate parents are excluded. The integration domain is the **final,
possibly expanded root**. Missing/invalid backtraces contribute zero.

For attraction the integrand has an integrable `1/|v|` singularity at zero.
The isolated zero-velocity lattice node is assigned zero, including when its
backtrace cannot reach a boundary. This choice avoids infinities, but needs a
resolution/convergence check for strong attraction. It is not an analytic
singular-cell correction.

For each potential the reported estimate is

```text
flux_error = sum_terminal_boxes (max sampled integrand - min sampled integrand) * box_volume
```

It is finite and nonnegative for ordinary finite inputs. This is a sampled
variation estimate, **not a statistical confidence interval or a rigorous error
bound**: it excludes unsampled narrow lobes, the origin singularity and tails
outside the final root. It does not estimate backtrace/time-step error.

Before splitting, the current evaluator reserves all eight children plus the
samples for an entire sibling batch. Within a finite budget, fully sampled
terminal boxes are prioritised by the largest difference (over psi) between
their own-lattice and corner-only trapezoid moments. This schedules refinement
using a quadrature error indicator; the probability split criterion is unchanged.
The output integral remains the own-lattice trapezoid rule. If a budget prevents a split, the
sampled parent stays terminal, so terminal boxes still partition the root.
Status `2` or `4` therefore means a coarser complete-domain moment, requiring
convergence checks. If the sample capacity cannot even hold the root lattice,
no moment is evaluated: status `2`, zero flux/error and zero counts. This budget
handling differs deliberately from the existing sample-returning octree, which
retains its original partial-output behaviour.

| Status | Meaning |
|---:|---|
| 0 | Completed requested refinement |
| 1 | No positive source-connected probability |
| 2 | Sample capacity prevents refinement (or root sampling) |
| 3 | Root expanded; significant edge signal may remain |
| 4 | Box capacity prevents refinement |
| 5 | Invalid velocity bounds; zero output |

Capacity statuses supersede the expansion status. Box capacity takes precedence
if both budgets stop splits. The root is not re-scored after the final allowed
expansion, matching the existing octree. Inspect status and convergence before
using a table in DRIFT.

## Terminal flags on octree results

`get_probabilities_octree(...).leaf_is_terminal` is a boolean `(nleaf,)` array.
Use `result.leaf_mask(i) & result.leaf_is_terminal` to select terminal boxes for
one position. Existing `leaf_*` fields retain their meaning and include
intermediate boxes. The original C symbol is unchanged; a new
`get_probabilities_octree_with_terminals` symbol adds flags from the traversal.
On completed trees their volumes sum to the final root volume. Original octree
capacity statuses `2/4` can leave an incomplete tree and do not guarantee this
partition. For saved results or older libraries without the new symbol, flags
are derived from parent/child containment of returned boxes; unsampled children
cannot be inferred on truncated trees (`metadata['terminal_flags']='returned_tree'`
for the older-library fallback).

Both APIs accept `field_bundle=` as a bundle directory, `bundle.json` path or
`FieldBundle` object. Directory, step, species and field flags must match the
bundle; mismatches raise `ValueError`. A bundle call never opens EMSES outputs
or rewrites its namelist.

## Field bundle format (0.1.0)

`vdsolverf.io.write_field_bundle(directory, istep, ispecs, output_dir, *,
use_electric_field=True, use_magnetic_field=True)` returns a `FieldBundle`.
`read_field_bundle(path)` reads metadata only; `bundle.get_ebvalues(ispec)` (or
`bundle[ispec]`) lazily loads a read-only `np.load(..., mmap_mode='r')` array.

| File | Contents |
|---|---|
| `plasma-vdsolverf.inp` | Exactly the namelist produced by `TempolaryInput`, including geometry conversion and indexed keys |
| `ebvalues_ispec{N}.npy` | Zero-based species `N`; relocated fields, shape `(nz+1, ny+1, nx+1, 9)` |
| `bundle.json` | `schema="vdsolverf.field_bundle"`, `version="0.1.0"`, absolute source `directory`, `istep`, `nx/ny/nz`, `ispecs`, both field flags, `created_at` (UTC ISO 8601), `vdsolverf_version` |

Arrays use **npy v1.0, little-endian float64 (`descr='<f8'`), C order**
(`fortran_order=False`). The first six bytes are `0x93 NUMPY`, followed by version
bytes `1,0`, a two-byte little-endian unsigned header length, and a padded ASCII
Python dictionary terminated by newline. Its keys are `descr`, `fortran_order`,
`shape`. Binary values follow in C order, field component fastest, then x, y, z.
There are no pickle objects or compression. The Fortran view of this data is
`ebvalues(9, nx+1, ny+1, nz+1)`. Read and validate the header length; do not assume
a fixed byte offset. The standalone Fortran reader accepts both v1.0 and v2.0.

## Standalone batch

`vdsolverf-batch [--resume|--no-resume] batch.nml` loads a field bundle and
constructs a simulator once per rank. Execution needs neither Python nor
mpi4py. MPI is enabled only in the dedicated build; the default `make`,
`make install` and pip shared-library build require no MPI. The dependency-free
Fortran npy reader/writer accepts v1.0/v2.0, C order and little-endian
`<f8` / `<i4` / `<i8`; other dtypes, big endian and Fortran order are rejected.

```bash
make batch                       # gfortran, OpenMP serial executable
make batch-mpi                   # prefer mpif90, fall back to mpiifx
# Other wrappers: make batch-mpi MPI_FC=mpif90 MPI_FLAGS='-O3 -fopenmp'
build/batch-serial/bin/vdsolverf-batch batch/batch.nml
```

The MPI executable is `build/batch-mpi/bin/vdsolverf-batch`. Separate prefixes
retain both binaries. Override `BATCH_FLAGS` / `MPI_FLAGS` and `BATCH_PREFIX` /
`MPI_PREFIX` as needed. Build and run on compute nodes at KUDPC.

Here `mpif90` uses GCC 8 with Intel MPI; this is the validated configuration.
Intel Fortran (`BATCH_FC=ifx` / `MPI_FC=mpiifx`) uses
`-O3 -qopenmp -fp-model precise`. The former nonzero-field nondeterminism was
traced to shared compiler-generated descriptors in ifx 2023.2, including an
allocatable probability-wrapper copy on every boundary collision. The solver
now evaluates the immutable probability object directly, copies polymorphic
particle components explicitly, and serializes solver construction/destruction.
The nonuniform-field regression evaluates 19 positions twice with eight
threads and compares every output bit, including octree samples and final
particle states, against a serial evaluation:

```bash
fpm test --flag "-fopenmp" --target test_thread_determinism
fpm test --compiler ifx --flag "-O3 -qopenmp -fp-model precise" --target test_thread_determinism
```

Use the same compiler and floating-point flags when requiring bit identity.
Cross-compiler or fast-math builds need numerical tolerances. Autorange's
optional mean/covariance diagnostics additionally use locked, unordered sums;
see [autorange reproducibility](autorange.en.md#reproducibility).

```python
import numpy as np
from vdsolverf.io import write_field_bundle, write_batch, read_batch_results

bundle = write_field_bundle("my-run", -1, [0], "field_bundle")
targets = np.array([
    [12.5, 8.5, 4.5, -5., 5., -5., 5., -5., 5.],
    [13.5, 8.5, 4.5, -5., 5., -5., 5., -5., 5.],
])
batch_nml = write_batch("batch", bundle, 0, targets, [-1., 0., 1.],
                        dt=0.1, max_step=30_000, chunk_size=16)
# After running vdsolverf-batch batch/batch.nml in an external job:
result = read_batch_results("batch")  # DustCurrentResult in original target order
```

`write_batch` accepts a bundle path or `FieldBundle` and returns the batch.nml
Path. Bundle paths passed to the helper resolve from the caller's cwd and are
stored relative to batch.nml. `output_dir` is relative to batch.nml, or absolute.
Omitting `dt` uses the bundle's `/tmgrid/ dt`. Omit `n_threads` (or pass `None`)
to use the OpenMP setting at execution time.

```text
batch/
  batch.nml
  targets.npy    # float64 C order (n,9): x,y,z,vxmin,vxmax,vymin,vymax,vzmin,vzmax
  psi.npy        # float64 C order (npsi,), EMSES potentials, order preserved
  out/
```

Both arrays must be nonempty. Positions and psi must be finite; invalid
velocity boxes return native status `5`. `ispec` is zero-based, as in Python.
All inputs and outputs use the EMSES units and unit-density flux defined above.

The namelist contains one `&vdsolverf_batch` group:

| Key | Default / meaning |
|---|---|
| `bundle_dir`, `ispec`, `dt` | Required in handwritten input; bundle path may be relative to batch.nml |
| `max_step`, `use_adaptive_dt`, `max_probability_types` | `30000`, `.false.`, `100` |
| `scout_bins`, `max_depth` | `9,9,9`, `6` |
| `max_samples_per_cell`, `max_leaves_per_cell` | `50000`, `4096` |
| `refine_threshold_rel`, `edge_threshold_rel` | Both `1e-4` |
| `expand_factor`, `max_expansions` | `1.5`, `3` |
| `chunk_size`, `output_dir` | `16`, `'out'` |
| `n_threads` | OpenMP runtime default, normally `OMP_NUM_THREADS` (1 without OpenMP) |
| `use_electric_field`, `use_magnetic_field` | Both `.true.`; must match bundle flags exactly |

Field-flag mismatches stop execution. Generate a separate bundle with the
desired flags to change which fields are used.

### Chunks and resume

Chunk start indices are **zero-based**. Only the final chunk may be shorter.

| File (example start=0) | dtype / C-order shape |
|---|---|
| `chunk_000000000_flux.npy` | `<f8` `(count,npsi)` |
| `chunk_000000000_flux_error.npy` | `<f8` `(count,npsi)` |
| `chunk_000000000_status.npy` | `<i4` `(count,)` |
| `chunk_000000000_sample_count.npy` | `<i8` `(count,)` |
| `chunk_000000000_terminal_box_count.npy` | `<i8` `(count,)` |
| `chunk_000000000.done` | Empty marker, atomically renamed last after closing all five arrays |

Resume is enabled by default (`--resume` is also accepted). A chunk is skipped
only when its marker and all five arrays exist. Without a marker, remaining
`.npy` / `.tmp` files are replaced by recomputation. `--no-resume` removes each
marker before recomputing all chunks. Keep the bundle, targets, psi, solver
settings and chunk_size immutable while resuming. After changing inputs, use
a fresh output directory or `--no-resume`; input freshness is not hashed.
Concurrent jobs writing the same output directory are unsupported.

After all ranks finish, rank 0 atomically writes `summary.json` with
`status_counts` over all output targets, `n_targets`, `npsi`,
`wall_time_seconds` (loading, computation and waiting), `ranks`, `threads`,
and `chunks_computed` / `chunks_skipped` for this invocation.
`read_batch_results` validates every marker, dtype and shape and rejects
partial results. The summary appears in `result.metadata['summary']`; collection
also works if all chunks finished just before a job died without a summary.
Batch diagnostic counts retain int64; existing ctypes counts remain int32.

### MPI and KUDPC

With two or more ranks, rank 0 only dispatches work. Each worker requests its
next chunk after finishing the previous one; rank 0 replies with a chunk index
or a termination instruction. Keeping the dispatcher out of the octree prevents
costly targets from delaying work distribution. The one-rank MPI and serial
builds compute all chunks on rank 0. Inside a worker, the existing evaluator
uses OpenMP to schedule targets dynamically.

```bash
# Compile on a compute node, launching the build controller once.
tssrun -p <partition> -t 0:40:00 --rsc p=1:t=4:c=4 make batch-mpi

# Add build/batch-mpi/bin to PATH before submission.
tssrun -p <partition> -t 2:00:00 --rsc p=16:t=8:c=8 srun vdsolverf-batch batch/batch.nml
```

The KUDPC smoke test on 2026-10-15 succeeded with this form, from the repository
root. GCC 8.5 + Intel MPI, four ranks and four threads per rank, produced
bit-identical arrays to serial for 19 targets / 10 chunks. `--resume` skipped
all ten chunks. No PMI environment changes or hydra fallback were needed.

```bash
tssrun -p gr20001b -t 0:40:00 --rsc p=4:t=4:c=4 \
  srun build/batch-mpi/bin/vdsolverf-batch --no-resume build/batch-smoke/mpi/batch.nml
```

Generate the fixture with `.venv/bin/python test/batch_mpi_smoke.py prepare
build/batch-smoke`; after running serial, use the script's `compare
build/batch-smoke` command to verify MPI output.

Set `n_threads` no higher than the allocated threads/cores. Every rank keeps
its own fields, including the dispatcher. Smaller chunks improve balancing
but increase message and file counts; 16 is a starting point. Threads beyond
the number of targets in a chunk cannot be used effectively. Scratch scales
with worker threads and octree budgets, rather than the total target count.
The existing octree also copies its solver/simulator for each active target;
include those field copies when estimating peak memory per worker.

## Performance and verification

OpenMP distributes spatial points dynamically; a velocity sweep is reused for
all potentials. Sample scratch scales as `n_threads * max_samples_per_cell`
(and box scratch as `n_threads * max_leaves_per_cell`), not total samples over
all points. Outputs scale as `nspatial * npsi`. Root-expansion scout evaluations
are not included in `sample_count`, matching the existing octree.

Use bundles to amortise EMSES loading and relocation across repeated calls.
Choose bounds using physical/thermal estimates or autorange, and check budget
statuses, scout spacing, tails, time step and refinement independently. A
one-percent Maxwellian regression with the default octree settings does not
guarantee one-percent accuracy for arbitrary VDFs or velocity bounds.

On KUDPC run backtrace work on a compute node, matching threads to resources:

```bash
tssrun -p gr20001b -t 0:40:00 --rsc p=1:t=8:c=8 .venv/bin/python make_currents.py
```

## C ABI

`evaluate_dust_currents` is a void `bind(c)` entry point. Argument order:

```text
inppath, length, lx, ly, lz, ebvalues, ispec, nspatial,
spatial_points, velocity_bounds, npsi, psi_grid, dt, max_step,
use_adaptive_dt, max_probability_types, scout_nvx, scout_nvy, scout_nvz,
max_depth, max_samples_per_cell, max_leaves_per_cell,
refine_threshold_rel, edge_threshold_rel, expand_factor, max_expansions,
flux, flux_error, status, sample_count, terminal_box_count, n_threads
```

All scalar integers are `integer(c_int), value`, all scalar reals are
`real(c_double), value`, except `n_threads`, which is an optional integer **by
reference** (ctypes passes `byref(c_int(...))`). The path is `character(c_char)`
bytes plus byte length. Real arrays are `real(c_double)`; integer output arrays
are `integer(c_int)`. Fortran shapes are `(3,nspatial)`, `(6,nspatial)`,
`psi_grid(npsi)`, moments `(npsi,nspatial)` and diagnostics `(nspatial)`.
Python uses C-contiguous reversed shapes, float64/int32, and converts `ispec+1`.
`get_probabilities_octree_with_terminals` has the original octree argument list
plus a final `integer(c_int)` output flag array (0/1), capacity
`nspatial * max_leaves_per_cell`. Existing C entry points retain their ABI.

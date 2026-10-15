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
a fixed byte offset. A standalone Fortran npy reader and MPI batch executable
are deferred to the later integration phase.

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

# Physics

> Lang: [日本語](physics.md) | **English**

This solver computes phase-space probability densities by combining
backward-in-time particle tracing with boundary-tagged source
distributions. The two central ideas are Liouville's theorem and the set
of `t_Probability` implementations in
[`src/core/probabilities.f90`](../src/core/probabilities.f90) and
[`src/emses/photoelectron_raycast.f90`](../src/emses/photoelectron_raycast.f90).

## Liouville's theorem

For a collisionless system in a force field that is a smooth function of
position and time, the one-particle distribution function
$f(\mathbf{x}, \mathbf{v}, t)$ is constant along particle trajectories:

$$
\frac{df}{dt} = \partial_t f + \mathbf{v} \cdot \nabla_{\mathbf{x}} f + \mathbf{a} \cdot \nabla_{\mathbf{v}} f = 0.
$$

Given an observation point $(\mathbf{x}_0, \mathbf{v}_0)$ at time $t_0$, we
trace the orbit backward in time. If the orbit reaches a boundary $\Sigma$
at $(\mathbf{x}_\Sigma, \mathbf{v}_\Sigma)$ with a known source distribution
$f_\Sigma$, then

$$
f(\mathbf{x}_0, \mathbf{v}_0, t_0) = f_\Sigma(\mathbf{x}_\Sigma, \mathbf{v}_\Sigma).
$$

`solver_calculate_probability` in
[`src/core/solver.f90`](../src/core/solver.f90) implements exactly this
lookup: it integrates the equations of motion backward, and on boundary
collision it dispatches to the boundary's tagged probability function.

## Boris updates and backward tracing

The magnetic rotation uses `s = 2*t/(1 + sum(t*t))`. The denominator is a
single scalar containing the squared norm of `t`. In the absence of E, this
preserves speed and the velocity component parallel to B for any field direction.

MPIEMSES3D's ordinary Boris step updates velocity using the field at the current
position, then drifts with the updated velocity. Backward tracing reverses this
order: first recover the previous position using the current velocity, then undo
the velocity update using the field at that recovered position. Field sampling
wraps the position into the periodic domain after a periodic boundary crossing.

A step and its inverse recover position and velocity to roundoff when they use
the same static field and the same time step, without boundary collisions. This
does not guarantee reversal of an entire trajectory with adaptive steps or
collision handling.

Collision position, velocity, and time use the same segment fraction `r`.
The time is `t_before + r * dt_step`, including the sign and any adaptive
change of the step duration.

## Supported source distributions

### Zero (absorbing)

`t_ZeroProbability%at` returns 0 unconditionally — used for walls that
trap backtraced particles without contributing any source.

### Shifted Maxwellian

`t_MaxwellianProbability%at` returns the 3D Maxwell–Boltzmann density

$$
f_M(\mathbf{v}) = C \prod_{i=1}^{3} \frac{1}{\sqrt{2\pi}\sigma_i}
  \exp\!\left(-\frac{u_i^2}{2 \sigma_i^2}\right),
$$

where $\boldsymbol\mu$ = `locs` is the drift velocity, $\boldsymbol\sigma$
= `scales` gives principal-axis thermal spreads, and $C$ = `coefficient` is
a multiplicative weight. Evaluation uses $\mathbf u = R^T(\mathbf v-\boldsymbol\mu)$.
The optional `basis` ($R$) defaults to the identity matrix.

EMSES external sources and raycast photoelectrons use `locs = vdri_vector(ispec)`
and `scales = [peth(ispec), peth(ispec), path(ispec)]`. `thermal_basis()` applies
the same `phiz`, `phixy` rotation as the drift to the principal axes. The
covariance $R\,\mathrm{diag}(\sigma_i^2)R^T$ retains correlations in oblique
anisotropic distributions and leaves isotropic distributions invariant.
`vth_vector(ispec)` returns the square roots of its diagonal entries, the
Cartesian marginal standard deviations used for autorange support. Axis-aligned
emission patches use `path` along the signed `nemd` normal and `peth` along
both tangents (standard deviations remain positive).

### Emission-patch normalisation

For an explicit emission patch, let $u=\mathrm{sign}(nemd)\,\mu_{|nemd|}$
be the drift toward the plasma along the emission direction and let
$\sigma=path$. The one-sided moment of the unit-density Gaussian is

$$
I(u,\sigma)=\frac{u}{2}\left[1+\mathrm{erf}\left(\frac{u}{\sqrt{2}\sigma}\right)\right]
+\frac{\sigma}{\sqrt{2\pi}}\exp\left(-\frac{u^2}{2\sigma^2}\right).
$$

The coefficient is `C = curf/I`, with nonnegative `curfs` overriding `curf`.
Thus the emitted one-sided moment recovers that input. MPIEMSES3D samples
the normal emitted velocities from the flux-weighted half-Gaussian $v_n f(v_n)$,
with Gaussian standard deviation `path`; tangents use `peth`. The transported
phase-space distribution here is $C f_M$, without another velocity factor.
`curf/curfs` are current-density magnitudes in EMSES units: MPIEMSES3D divides
them by the macro-particle charge magnitude to obtain particle injection counts.
Here $C$ is a charge-equivalent density (number density for unit physical charge
magnitude); macro-particle charge/weight must not be applied to this continuous
distribution. For another physical charge magnitude, divide by it to obtain
number density. This source is not a density-one PDF.

The erf utility uses width $\sqrt{2}\sigma$, rather than $\sigma$.
Correcting the former swapped arguments and width convention gives
$C=curf\sqrt{2\pi}/\sigma$ at zero drift. For the same previously positive
normal spread (e.g. +/-Z patches), probabilities rise by
$\sqrt{2\pi}\simeq2.5066$. X/Y patches additionally correct the former
signed/permuted thermal spreads; formerly zero sources can become nonzero,
so a finite common scale factor does not describe those cases.
See [dust-current units](dust_currents.en.md#definition-and-units).

### Raycast photoelectron

Activated by `use_raycast = .true.` with `nflag_emit(ispec) == 2`.
Interpretation: a backtraced particle hitting an internal surface may have
been a photoelectron freshly emitted from that surface. Two physical
conditions gate the probability.

**1. Outward emission half-space.** Photoelectrons are emitted from the
surface pointing away from the material. Approximating the outward normal
with the sun direction $\hat{\mathbf n}$ (unit vector from surface toward
the sun), the particle's velocity must satisfy $\mathbf{v} \cdot
\hat{\mathbf n} > 0$. Otherwise the density is 0.

**2. Unoccluded illumination.** The photoelectron only exists if sunlight
reached the surface. A ray from the collision point in direction
$\hat{\mathbf n}$ is cast into a blocking-boundary list (the same internal
surfaces/objects that could cast shadows). If the ray hits anything with
$t > 0$ the surface is shaded and the density is 0.

When both conditions pass, the probability is

$$
f_{\rm PE}(\mathbf{v}) = 2 \cdot C \prod_{i=1}^{3}
  \frac{1}{\sqrt{2\pi}\sigma_i}
  \exp\!\left(-\frac{u_i^2}{2 \sigma_i^2}\right),
$$

with $\boldsymbol\mu$ = `vdri_vector(ispec)` and $\boldsymbol\sigma$ =
`[peth(ispec), peth(ispec), path(ispec)]`; $\mathbf u$ uses the same inverse
rotation into principal axes as above. The factor of 2 is the exact half-space normalization
when $\mu_\parallel = 0$; for shifted distributions it is a first-order
approximation that remains accurate as long as $|\boldsymbol\mu \cdot
\hat{\mathbf n}|$ is not dominated by the thermal spread.

The sun direction is resolved by `resolve_sun_direction(ispec)` in
[`src/emses/emses_simulator_builder.f90`](../src/emses/emses_simulator_builder.f90):

1. If `ray_zenith_angle_deg(ispec) < 9000d0`, use it; otherwise fall back
   to `vdthz(ispec)`. Same for azimuth / `vdthxy`.
2. Start from $[0, 0, 1]$, rotate by $-\zeta$ around the $y$-axis, then
   $\phi$ around $z$.
3. Return the negation of the resulting unit vector (pointing from the
   emission surface toward the sun, per the user-requested convention).

## Coordinate and unit conventions

- All positions are in EMSES grid units ($[0, nx] \times [0, ny] \times
  [0, nz]$).
- Velocities are in EMSES natural units; the magnitude convention matches
  `vdri`, `path`, `peth` from `/intp/`.
- Boundaries codes (`npbnd`): `0` periodic, `1` reflective, `2` open
  (Maxwellian source candidate), `3` absorbing.
- Species index `ispec` is **1-based in Fortran** and **0-based in the
  Python API** — the wrapper adjusts automatically.

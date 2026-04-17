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

## Supported source distributions

### Zero (absorbing)

`t_ZeroProbability%at` returns 0 unconditionally — used for walls that
trap backtraced particles without contributing any source.

### Shifted Maxwellian

`t_MaxwellianProbability%at` returns the 3D Maxwell–Boltzmann density

$$
f_M(\mathbf{v}) = C \prod_{i=1}^{3} \frac{1}{\sqrt{2\pi}\sigma_i}
  \exp\!\left(-\frac{(v_i - \mu_i)^2}{2 \sigma_i^2}\right),
$$

where $\boldsymbol\mu$ = `locs` is the drift velocity, $\boldsymbol\sigma$
= `scales` is the thermal spread, and $C$ = `coefficient` is a multiplicative
weight. The builder populates these from EMSES parameters via
`vdri_vector(ispec)` and `vth_vector(ispec)` in
[`src/emses/allcom.f90`](../src/emses/allcom.f90).

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
  \exp\!\left(-\frac{(v_i - \mu_i)^2}{2 \sigma_i^2}\right),
$$

with $\boldsymbol\mu$ = `vdri_vector(ispec)` and $\boldsymbol\sigma$ =
`vth_vector(ispec)`. The factor of 2 is the exact half-space normalization
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

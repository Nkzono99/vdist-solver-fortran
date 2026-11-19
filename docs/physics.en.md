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

The EMSES sources below follow the formulas MPIEMSES3D uses to generate particles.
[Synchronising with EMSES](emses_sync.en.md) lists the corresponding EMSES sources
and the update procedure.

**Outer-boundary injection** (`npbnd = 2`, `nflag_emit = 0`) uses
`locs = vdri_vector(ispec)`, `scales = [peth, peth, path]` and `basis = thermal_basis()`.
As in the EMSES initial loading and outer-boundary injection (4.23.0 and later,
`m_thermal_frame.F90`), the thermal spread is defined in the field frame with covariance
$R\,\mathrm{diag}(peth^2, peth^2, path^2)R^T$; oblique anisotropic distributions keep their
Cartesian correlations and isotropic ones are rotation invariant. The drift equals EMSES
`vdtx/vdty/vdtz`:

$$
\boldsymbol\mu = (vdx, vdy, vdz) + R\,(spe\cos(speth),\ spe\sin(speth),\ spa)
+ vdri\,(\sin\vartheta\cos\varphi,\ \sin\vartheta\sin\varphi,\ \cos\vartheta)
$$

with $\vartheta$ = `vdthz`, $\varphi$ = `vdthxy`; all angles are in degrees.
$R$ = `thermal_basis()` holds the EMSES background-field frame (`t11`..`t33`) as
columns; its third column equals the background field packed by the Python wrapper,
$\hat{\mathbf B}_0 = (\sin\phi_z\cos\phi_{xy},\ \sin\phi_z\sin\phi_{xy},\ \cos\phi_z)$.
`vth_vector(ispec)` returns the square roots of the covariance diagonal, the Cartesian
marginal standard deviations that EMSES also uses as boundary-table widths; autorange uses
them as support widths.

### Emission-patch normalisation

Surface sources (emission patches and raycast photoelectrons) use, in local
coordinates (outward normal, tangent 1, tangent 2), the mean `[spa + abvdem, spe, spe]`
and the standard deviations `[path, peth, peth]`. EMSES draws the normal velocity from
a flux-weighted one-sided Gaussian table with mean `spa + abvdem` and both tangents
independently from one Gaussian table with mean `spe`; `speth` is not used
(`surface_emission_source.F90`, `direct_surface_emitter.F90`). The tangents keep
their sign for either patch orientation (sign of `nemd`).

With $u = spa + abvdem$ and $\sigma = path$, the one-sided moment of the
unit-density Gaussian is

$$
I(u,\sigma)=\frac{u}{2}\,\mathrm{erfc}\left(-\frac{u}{\sqrt{2}\sigma}\right)
+\frac{\sigma}{\sqrt{2\pi}}\exp\left(-\frac{u^2}{2\sigma^2}\right).
$$

$\mathrm{erfc}$ replaces $1+\mathrm{erf}$ so that a drift toward the surface
($u \ll -\sigma$) does not cancel. The coefficient is `C = J/I` and the source
returns $C f_M$ only in the outward half-space $\mathbf v\cdot\hat{\mathbf n} > 0$,
so its outward first moment recovers $J$.

$J$ is chosen in the EMSES order: `curfs(iepl)` if positive, otherwise `curf(ispec)`
if positive. When neither is positive, EMSES falls back to `flpf`, `flpfs` or an
`npin`-based density, which this solver does not read, so the input is rejected.
For photoelectron patches (`nflag_emit = 2` with `pe_ray_cast = .false.`) and
`decrease_current_according_to_incident_angle` (default `.true.`), $J$ is multiplied by
$\mathrm{sign}(nemd)\,\hat d_{|nemd|}$ and set to zero when negative, as in EMSES;
$\hat{\mathbf d}$ is the incidence direction from the scalars `thetaz` and `thetaxy`
(degrees). EMSES shifts patch-emitted velocities by `vdri`, which moves the
half-space; species with patches and `vdri ≠ 0` are rejected.

`curf/curfs` are current densities in EMSES units: MPIEMSES3D divides them by the
macro-particle charge magnitude to obtain particle injection counts. Here $C$ is a
charge-equivalent density (number density for unit physical charge magnitude);
macro-particle charge/weight must not be applied to this continuous distribution.
For another physical charge magnitude, divide by it to obtain number density. This
source is not a density-one PDF. See
[dust-current units](dust_currents.en.md#definition-and-units).

### Raycast photoelectron

Activated by `use_raycast = .true.` (`pe_ray_cast` in EMSES inputs) with
`nflag_emit(ispec) == 2`. EMSES casts rays from a parallel camera along the sun
direction; at the first emitting surface hit (`material%tag /= -1`) it emits one
particle from the surface source above, using the hit normal $\hat{\mathbf n}$ and
the tangent basis of `raycast_sampler.F90`. Rays per projected area are fixed, so the
emitted flux per surface area is $|curf|\max(\hat{\mathbf n}\cdot\hat{\mathbf s}, 0)$,
with $\hat{\mathbf s}$ the unit vector from the surface toward the sun.

This solver casts a ray back along the camera direction from a point just sunward
of the collision point to recover the hit face and its normal, and returns

$$
f_{\rm PE}(\mathbf v) = \frac{|curf|\,\max(\hat{\mathbf n}\cdot\hat{\mathbf s},0)}{I}
\prod_{i=1}^{3}\frac{1}{\sqrt{2\pi}\sigma_i}\exp\!\left(-\frac{(w_i-\mu_i)^2}{2\sigma_i^2}\right),
\qquad w_1 = \mathbf v\cdot\hat{\mathbf n} > 0,
$$

where $\mathbf w$ holds the (normal, tangent 1, tangent 2) components and
$\boldsymbol\mu$, $\boldsymbol\sigma$ are the surface-source values above.
The density is zero when:

- the collision point is not on an emitting surface, or lies on the EMSES
  non-emitting cover (tag -1);
- the velocity points into the surface ($w_1 \le 0$);
- the ray toward the sun hits another boundary (shadow).

The occlusion list is built with `cover_all`, as in EMSES. The sun direction is the
reverse of the ray direction of `new_ParallelCamera_optimized(thz, thxy, domain)`,
the camera EMSES uses. An unset `ray_zenith_angle_deg(ispec)` (9999) falls back to
`vdthz(1)` and an unset `ray_azimuth_angle_deg(ispec)` to `vdthxy(1) - 180`, both
from species 1. The camera rotation uses a different sign convention from
`direction_from_angles`, so the direction is taken from the camera itself.
Periodic raycasts (`periodic_raycast`) and time-varying emission
(`emission_time_mode`) are not reproduced.

## Coordinate and unit conventions

- All positions are in EMSES grid units ($[0, nx] \times [0, ny] \times
  [0, nz]$).
- Velocities are in EMSES natural units; the magnitude convention matches
  `vdri`, `path`, `peth` from `/intp/`.
- Boundaries codes (`npbnd`): `0` periodic, `1` reflective, `2` open
  (Maxwellian source candidate), `3` absorbing.
- Species index `ispec` is **1-based in Fortran** and **0-based in the
  Python API** — the wrapper adjusts automatically.

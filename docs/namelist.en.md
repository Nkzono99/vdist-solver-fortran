# Namelist reference

> Lang: [日本語](namelist.md) | **English**

`vdist-solver-fortran` reads a subset of the EMSES `plasma.inp` file. The
Python wrapper writes the relevant groups to a temporary file
(`plasma-vdsolverf.inp`) before calling the Fortran shared library — only
keys listed in [`tmpolary_input.TMP_INP_KEYS`](../vdsolverf/emses/tmpolary_input.py)
are forwarded.

Unsupported keys are silently ignored. Adding a new key requires editing both
`TMP_INP_KEYS` (Python side) and `m_namelist` (Fortran side).

## Groups

| Group | Purpose |
|---|---|
| `/esorem/` | `emflag` — electromagnetic mode switch |
| `/plasma/` | plasma parameters (`wp`, `wc`, background B-field angles `phixy`, `phiz`) |
| `/tmgrid/` | time step `dt`, grid sizes `nx`, `ny`, `nz` |
| `/system/` | species count `nspec`, boundary codes `npbnd(3, nspec)` |
| `/intp/` | per-species kinetics (`qm`, `path`, `peth`, `vdri`, `vdthz`, `vdthxy`, `spa`, `spe`, `speth`, `vdx`, `vdy`, `vdz`) |

All angles (`phiz`, `phixy`, `vdthz`, `vdthxy`, `speth`) are in degrees and follow
the EMSES conventions: the background field points along
`(sin(phiz)cos(phixy), sin(phiz)sin(phixy), cos(phiz))`, `spa`/`path` are parallel
to it, `spe` lies in the perpendicular plane at angle `speth`, and `vdri` has zenith
angle `vdthz` from +z and azimuth `vdthxy` from +x; `vdx`/`vdy`/`vdz` are a Cartesian
drift (see [Physics](physics.en.md)).
| `/ptcond/` | inner-boundary geometry and reflection properties |
| `/emissn/` | emission-surface and raycast photoelectron settings |

The most consequential groups for this solver are `/ptcond/` (which
boundaries exist) and `/emissn/` (how probability is assigned at collisions).

## /ptcond/ — Inner boundaries

Two entry points coexist: `boundary_type` / `boundary_types(:)` (typed
surface library) and `geotype` (simpler primitive shapes).

### boundary_type

Single boundary type. Accepted values:

```
'none'
'flat-surface'
'rectangle-hole' | 'cylinder-hole' | 'hyperboloid-hole' | 'ellipsoid-hole'
'rectangle[xyz]' | 'circle[x/y/z]' | 'cuboid' | 'disk[x/y/z]'
'complex'          ! pick multiple from boundary_types(:)
```

### Parameters by boundary kind

| Kind | Keys |
|---|---|
| `flat-surface` / `*-hole` | `zssurf` (surface height, grid units); holes use `[x/y/z][l/u]pc` |
| `complex` | `boundary_types(ntypes)` — list of types (see above) |
| `rectangle` | `rectangle_shape(ntypes, 6)` = `(xmin, xmax, ymin, ymax, zmin, zmax)` |
| `rectangle-corners` | `rectangle_corners(3, 4, ntypes)` = the four corners of an oriented rectangle |
| `circle[x/y/z]` | `circle_origin(ntypes, 3)`, `circle_radius(ntypes)` |
| `cuboid` | `cuboid_shape(ntypes, 6)` = `(xmin, xmax, ymin, ymax, zmin, zmax)` |
| `disk[x/y/z]` | `disk_origin(ntypes, 3)`, `disk_height(ntypes)`, `disk_radius(ntypes)`, `disk_inner_radius(ntypes)` |
| Global rotation | `boundary_rotation_deg(3)` (degrees) |

### geotype (simplified primitives)

| Key | Meaning |
|---|---|
| `npc` | Number of geotype objects |
| `geotype(npc)` | `0`–`1` cuboid, `2` cylinder, `3` sphere |
| Cuboid | `xlpc`, `xupc`, `ylpc`, `yupc`, `zlpc`, `zupc` |
| Cylinder | `bdyalign` (1=X, 2=Y, 3=Z), `bdyedge(1:2)` (axial bounds), `bdyradius`, `bdycoord(1:2)` (center on axis) |
| Sphere | `bdyradius`, `bdycoord(1:3)` |

## /emissn/ — Particle emission

```
nflag_emit(nspec)  = 0 outer-boundary injection, 1 surface emission,
                     2 photoelectron (3 = SEE is not supported)
nepl(nspec)        = # of explicit emission surfaces for this species
nemd(nepl)         = surface normal; sign = direction, magnitude = axis
                     (1=X, 2=Y, 3=Z)
curf(nspec)        = emission current density (patches use it only when positive)
curfs(nepl)        = per-surface current (overrides curf when positive)
abvdem(nspec)      = extra mean velocity along the emission normal
xmine, xmaxe, ...  = emission-surface bounding box per nepl
thetaz, thetaxy    = incidence direction of photoelectron patches (scalars [deg])
decrease_current_according_to_incident_angle
                   = scale photoelectron patch currents by the incidence cosine
                     (default .true.)
```

See [Physics](physics.en.md#emission-patch-normalisation) for the patch mean and
thermal velocities, the current selection and the incidence scaling. Settings whose
EMSES behaviour cannot be replayed raise `ValueError`: patches without a positive
`curf` or `curfs`, patch-emitting species with `vdri ≠ 0`, SEE, reflection,
loss-cone, periodic raycasts, `omniemit`, and two-component emission with
`mainprop < 1`.

### Raycast photoelectron

When `use_raycast = .true.` (Python converts the EMSES key `pe_ray_cast`)
**and** `nflag_emit(ispec) == 2`, internal-
boundary collisions are evaluated by the raycast photoelectron probability
defined in [`m_photoelectron_raycast`](../src/emses/photoelectron_raycast.f90)
and modelled in [Physics](physics.en.md#raycast-photoelectron). Without the
flag, the solver falls back to the default ZeroProbability (absorbing).

| Key | Default | Meaning |
|---|---|---|
| `use_raycast` | `.false.` | Master switch. Enables raycast probability for photoelectron species. |
| `ray_zenith_angle_deg(nspec)` | `9999d0` | EMSES camera zenith angle. The sentinel `9999d0` falls back to `vdthz(1)` (species 1). |
| `ray_azimuth_angle_deg(nspec)` | `9999d0` | EMSES camera azimuth. The sentinel falls back to `vdthxy(1) - 180`. |

The sun direction comes from the same finbound camera EMSES uses. A shaded collision
point has probability 0; otherwise the surface-source density with the hit-face
normal is returned ([Physics](physics.en.md#raycast-photoelectron)).

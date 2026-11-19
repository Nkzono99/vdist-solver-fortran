#!/usr/bin/env bash
# Compare vdist's EMSES-derived code with an MPIEMSES3D checkout.
#
# usage: scripts/check-emses-sync.sh /path/to/MPIEMSES3D [emses-ref]
#
# 1. The mirrored collision files must be byte-identical to <emses-ref>
#    (default HEAD). A difference means EMSES changed them: copy the new
#    files, run `fpm test`, and update SYNCED_COMMIT below.
# 2. The re-implemented EMSES behaviour (sources, raycast, drift) is listed
#    as the EMSES commits since SYNCED_COMMIT that touch WATCHED_PATHS.
#    Review each one against docs/emses_sync.md and port what applies.
#
# This file is the single record of the last reviewed EMSES commit.
set -euo pipefail

SYNCED_COMMIT=2a9a44a5fd2838f96bc8b9910dc4e2a99f070d53  # MPIEMSES3D v4.23.0

# vdist path : EMSES path (byte-identical copies)
MIRRORED=(
  "src/emses/collision/objects.F90:src/physics/collision/objects.F90"
  "src/emses/collision/surfaces.F90:src/physics/collision/surfaces.F90"
  "src/emses/collision/boundary_types.F90:src/physics/collision/boundary_types.F90"
  "src/emses/collision/boundary_helpers.F90:src/physics/collision/boundary_helpers.F90"
  "src/emses/collision/collision_accelerator.F90:src/physics/collision/collision_accelerator.F90"
)

# EMSES sources whose behaviour vdist re-implements (see docs/emses_sync.md).
WATCHED_PATHS=(
  src/main/config/collision_boundary_input_config.F90
  src/main/initialization/simulation/initialization_simulation_setup.F90
  src/physics/particles/generation/sources/domain_boundary_source.F90
  src/physics/particles/generation/sources/surface_emission_source.F90
  src/physics/particles/generation/surface/direct_surface_emitter.F90
  src/physics/particles/generation/surface/raycast_surface_emitter.F90
  src/physics/particles/generation/surface/raycast_sampler.F90
  src/core/utils/m_velocity_distribution.F90
  src/core/utils/m_thermal_frame.F90
  src/physics/fields/common/m_trilinear_interp.F90
  src/physics/particles/solver/particle_solver.F90
  src/main/config/input/namelist_reader.F90
)

if [ $# -lt 1 ]; then
  sed -n '2,12p' "$0" >&2
  exit 2
fi
emses=$1
ref=${2:-HEAD}
root=$(cd "$(dirname "$0")/.." && pwd)

status=0
echo "== mirrored files (vdist vs EMSES ${ref})"
for pair in "${MIRRORED[@]}"; do
  local_path=${pair%%:*}
  emses_path=${pair##*:}
  if git -C "$emses" show "${ref}:${emses_path}" | cmp -s - "$root/$local_path"; then
    echo "  same     $local_path"
  else
    echo "  DIFFERS  $local_path  <-  $emses_path"
    status=1
  fi
done

echo "== EMSES commits since ${SYNCED_COMMIT} touching re-implemented behaviour"
git -C "$emses" log --oneline "${SYNCED_COMMIT}..${ref}" -- "${WATCHED_PATHS[@]}"

exit $status

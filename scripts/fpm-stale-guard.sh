#!/usr/bin/env bash
# Run an fpm command with incremental-build protection.
#
# usage: scripts/fpm-stale-guard.sh <fpm> <action> [args...]   (run from the repository root)

# fpm (0.12) skips a target whose source digest matches its build/**/*.digest
# unless one of its dependencies is rebuilt in the same run. A run that fails or
# is killed midway has already refreshed some modules, so their unchanged
# dependents would keep the old .mod interface and derived-type layout in every
# later run. Leave a marker while fpm runs; when a previous run left it behind,
# delete the digests written since then so fpm recompiles those targets and
# propagates to their dependents. A failing `fpm test`/`fpm run` payload also
# leaves the marker, which only costs recompiling what that run compiled.
fpm_stale_guard_run() {
  local marker="build/.fpm-run-unfinished" status
  mkdir -p build
  if [ -e "$marker" ]; then
    echo "[fpm-stale-guard] Previous fpm run did not finish; invalidating the objects it refreshed."
    find build -name '*.digest' -newer "$marker" -delete
  fi
  # Backdate the marker so digests written in the same timestamp tick count as newer.
  touch -d "@$(($(date +%s) - 2))" "$marker"
  "$@"
  status=$?
  if [ "$status" -eq 0 ]; then
    rm -f "$marker"
  fi
  return "$status"
}

fpm_stale_guard_run "$@"

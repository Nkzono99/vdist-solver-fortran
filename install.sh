#!/usr/bin/env bash
# Build and install the Fortran shared library used by the Python wrapper.
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${ROOT_DIR}/env"
LIBNAME="vdist-solver-fortran"
PACKAGE_DIR="${ROOT_DIR}/vdsolverf"
cd "${ROOT_DIR}"

detect_host_profile() {
  local host_name
  host_name="$(uname -n 2>/dev/null || true)"
  case "${host_name}" in
    camphor*) echo "camphor" ;;
    *) echo "generic" ;;
  esac
}

load_env_file() {
  local env_file="$1"
  if [[ ! -f "${env_file}" ]]; then
    return 1
  fi
  # shellcheck disable=SC1090
  source "${env_file}"
  return 0
}

load_modules_if_requested() {
  local module_name
  if [[ "${USE_MODULES}" != "1" ]]; then
    return 0
  fi

  if ! type module >/dev/null 2>&1 && [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  fi
  if ! type module >/dev/null 2>&1; then
    echo "[install.sh] warning: module command not found; continuing with current environment." >&2
    return 0
  fi

  for module_name in ${MODULES_TO_UNLOAD:-}; do
    module unload "${module_name}" >/dev/null 2>&1 || true
  done
  for module_name in ${MODULES_TO_LOAD:-}; do
    module load "${module_name}"
  done
}

setup_generic_defaults() {
  : "${USE_MODULES:=0}"
  : "${FC:=gfortran}"
  : "${FFLAGS:=-O3 -fPIC -fopenmp}"
  : "${SHARED_LDFLAGS:=-fopenmp}"
}

shared_suffix() {
  case "$(uname -s)" in
    Linux*) echo "so" ;;
    Darwin*) echo "dylib" ;;
    MINGW*|MSYS*|CYGWIN*) echo "dll" ;;
    *) echo "so" ;;
  esac
}

link_shared_library() {
  local suffix="$1"
  local static_lib="${PREFIX}/lib/lib${LIBNAME}.a"
  local shared_lib="${PREFIX}/lib/lib${LIBNAME}.${suffix}"

  if [[ ! -f "${static_lib}" ]]; then
    echo "[install.sh] error: static library not found: ${static_lib}" >&2
    exit 1
  fi

  mkdir -p "${PREFIX}/lib"
  case "${suffix}" in
    so)
      # shellcheck disable=SC2086
      "${FC}" -shared -o "${shared_lib}" \
        -Wl,--whole-archive "${static_lib}" -Wl,--no-whole-archive ${SHARED_LDFLAGS}
      ;;
    dylib)
      # shellcheck disable=SC2086
      "${FC}" -dynamiclib -install_name "lib${LIBNAME}.dylib" -o "${shared_lib}" \
        -Wl,-all_load "${static_lib}" -Wl,-noall_load ${SHARED_LDFLAGS}
      ;;
    dll)
      # shellcheck disable=SC2086
      "${FC}" -shared -static -o "${shared_lib}" \
        -Wl,--out-implib="${PREFIX}/lib/lib${LIBNAME}.dll.a",--export-all-symbols,--enable-auto-import,--whole-archive \
        "${static_lib}" -Wl,--no-whole-archive ${SHARED_LDFLAGS}
      ;;
    *)
      echo "[install.sh] error: unsupported shared library suffix: ${suffix}" >&2
      exit 1
      ;;
  esac
}

: "${BUILD_PROFILE:=auto}"
REQUESTED_PROFILE="${BUILD_PROFILE}"
if [[ "${BUILD_PROFILE}" == "auto" ]]; then
  BUILD_PROFILE="$(detect_host_profile)"
fi

if [[ -f "${ENV_DIR}/common.env" ]]; then
  load_env_file "${ENV_DIR}/common.env"
fi

if [[ "${BUILD_PROFILE}" != "generic" ]]; then
  if ! load_env_file "${ENV_DIR}/${BUILD_PROFILE}.env"; then
    if [[ "${REQUESTED_PROFILE}" == "auto" ]]; then
      echo "[install.sh] warning: profile env not found (${BUILD_PROFILE}). Fallback to generic." >&2
      BUILD_PROFILE="generic"
    else
      echo "[install.sh] error: profile env not found (${ENV_DIR}/${BUILD_PROFILE}.env)." >&2
      exit 1
    fi
  fi
fi

setup_generic_defaults

: "${PROFILE:=release}"
: "${PREFIX:=${ROOT_DIR}}"
: "${COPY_TO_PACKAGE:=1}"

load_modules_if_requested

export FPM_FC="${FC}"
export FPM_FFLAGS="${FFLAGS}"
export FPM_LDFLAGS="${SHARED_LDFLAGS}"

echo "[install.sh] BUILD_PROFILE=${BUILD_PROFILE}"
echo "[install.sh] FPM_FC=${FPM_FC}"
echo "[install.sh] FPM_FFLAGS=${FPM_FFLAGS}"
echo "[install.sh] PREFIX=${PREFIX}"

if ! command -v "${FPM_FC}" >/dev/null 2>&1; then
  echo "[install.sh] error: compiler not found: ${FPM_FC}" >&2
  exit 1
fi

fpm install --profile "${PROFILE}" --compiler "${FPM_FC}" --flag "${FPM_FFLAGS}" --prefix "${PREFIX}"

suffix="$(shared_suffix)"
link_shared_library "${suffix}"

if [[ "${COPY_TO_PACKAGE}" != "0" ]]; then
  mkdir -p "${PACKAGE_DIR}"
  rm -f "${PACKAGE_DIR}/lib${LIBNAME}.so" \
        "${PACKAGE_DIR}/lib${LIBNAME}.dylib" \
        "${PACKAGE_DIR}/lib${LIBNAME}.dll"
  cp "${PREFIX}/lib/lib${LIBNAME}.${suffix}" "${PACKAGE_DIR}/"
  echo "[install.sh] copied: ${PACKAGE_DIR}/lib${LIBNAME}.${suffix}"
fi

echo "[install.sh] done: ${PREFIX}/lib/lib${LIBNAME}.${suffix}"

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

from setuptools import find_packages, setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.develop import develop as _develop
from setuptools.command.install import install as _install

try:
    from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
except Exception:
    try:
        from wheel.bdist_wheel import bdist_wheel as _bdist_wheel
    except Exception:
        _bdist_wheel = None


ROOT_DIR = Path(__file__).resolve().parent
INSTALL_PREFIX = Path(
    os.environ.get("VDSOLVERF_PIP_PREFIX", ROOT_DIR / "build" / "pip-install")
)
LIBNAME = "vdist-solver-fortran"
PACKAGE_DIR = ROOT_DIR / "vdsolverf"
_BUILT_LIBRARY = None


def _library_suffix() -> str:
    if sys.platform.startswith("linux"):
        return ".so"
    if sys.platform == "darwin":
        return ".dylib"
    if sys.platform.startswith("win"):
        return ".dll"
    raise RuntimeError(f"Unsupported platform for vdist-solver-fortran: {sys.platform}")


def _package_library() -> Path:
    return PACKAGE_DIR / f"lib{LIBNAME}{_library_suffix()}"


def _which(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def _clean_package_libraries() -> None:
    for path in PACKAGE_DIR.glob(f"lib{LIBNAME}.*"):
        path.unlink()


def _run_make_install(install_profile: str) -> None:
    command = [
        "make",
        "install",
        f"INSTALL_PROFILE={install_profile}",
        f"PREFIX={INSTALL_PREFIX}",
    ]
    subprocess.check_call(command, cwd=ROOT_DIR)


def _build_with_make() -> Path:
    if not _which("make"):
        print("\nERROR: 'make' is required to build vdist-solver-fortran.\n", file=sys.stderr)
        sys.exit(1)

    _clean_package_libraries()

    install_profile = os.environ.get("INSTALL_PROFILE", "auto")
    try:
        _run_make_install(install_profile)
    except subprocess.CalledProcessError as exc:
        use_fallback = os.environ.get("VDSOLVERF_PIP_FALLBACK_GENERIC", "1") == "1"
        explicit_profile = "INSTALL_PROFILE" in os.environ
        if use_fallback and install_profile == "auto" and not explicit_profile:
            print(
                "\nWARN: auto profile build failed; retrying with INSTALL_PROFILE=generic.\n",
                file=sys.stderr,
            )
            try:
                _run_make_install("generic")
            except subprocess.CalledProcessError as retry_exc:
                print(
                    "\nERROR: failed to build Fortran shared library via make.\n"
                    "       Ensure fpm, make, and a Fortran compiler are available in PATH.\n",
                    file=sys.stderr,
                )
                raise SystemExit(retry_exc.returncode) from retry_exc
        else:
            print(
                "\nERROR: failed to build Fortran shared library via make.\n"
                "       Ensure fpm, make, and a Fortran compiler are available in PATH.\n",
                file=sys.stderr,
            )
            raise SystemExit(exc.returncode) from exc

    libpath = _package_library()
    if not libpath.exists():
        print(f"\nERROR: expected shared library not found: {libpath}\n", file=sys.stderr)
        sys.exit(1)

    return libpath


def _ensure_built_library() -> Path:
    global _BUILT_LIBRARY
    if _BUILT_LIBRARY is not None and _BUILT_LIBRARY.exists():
        return _BUILT_LIBRARY
    _BUILT_LIBRARY = _build_with_make()
    return _BUILT_LIBRARY


class build_py(_build_py):
    def run(self) -> None:
        _ensure_built_library()
        super().run()


class install(_install):
    def run(self) -> None:
        _ensure_built_library()
        super().run()


class develop(_develop):
    def run(self) -> None:
        _ensure_built_library()
        super().run()


if _bdist_wheel is not None:

    class bdist_wheel(_bdist_wheel):
        def finalize_options(self) -> None:
            super().finalize_options()
            self.root_is_pure = False

else:
    bdist_wheel = None


cmdclass = {
    "build_py": build_py,
    "install": install,
    "develop": develop,
}
if bdist_wheel is not None:
    cmdclass["bdist_wheel"] = bdist_wheel


setup(
    cmdclass=cmdclass,
    packages=find_packages(),
    package_data={"vdsolverf": [f"lib{LIBNAME}.*"]},
)

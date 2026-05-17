---
name: build-lib
description: Build the Fortran core with fpm and refresh the platform-specific shared library under vdsolverf/ so the Python ctypes wrapper loads current code. Use after changing src/**/*.f90, fpm.toml, Makefile, or the Python package build path, and before Python smoke tests that need the compiled library.
---

# Build Library

Use this skill when Fortran changes must be reflected in the Python wrapper's
shared library.

## Workflow

1. Check the build inputs before running anything expensive:
   ```bash
   git status --short
   fpm --version
   ```
2. Build the Fortran project with position-independent code:
   ```bash
   fpm build --flag "-fPIC"
   ```
3. Refresh the Python package library through the repository Makefile:
   ```bash
   make
   ```
4. Confirm that the expected OS-specific library exists:
   ```bash
   ls -l vdsolverf/libvdist-solver-fortran.*
   ```
5. Run the import smoke test with the repository virtualenv:
   ```bash
   .venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
   ```

## Notes

- Linux expects `vdsolverf/libvdist-solver-fortran.so`; macOS expects
  `.dylib`; Windows expects `.dll`.
- `setup.py` shells out to `make`, so Python package builds require the
  Fortran toolchain too.
- If the compile fails, report the first Fortran module and diagnostic that
  caused the failure. Avoid pasting the full build log unless the user asks.

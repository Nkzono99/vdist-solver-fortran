---
name: sync-wrapper-interface
description: Audit and update the Python ctypes wrapper in vdsolverf/emses/wrapper.py against Fortran bind(c) entry points in src/emses/emses_solver.f90. Use whenever a C-facing Fortran argument, type, value/reference convention, array shape, return buffer, or public wrapper function changes.
---

# Sync Wrapper Interface

Use this skill to keep the Python-to-Fortran ABI coherent. A mismatch here can
compile cleanly but crash at runtime.

## Scope

- Python: `vdsolverf/emses/wrapper.py`
  - `dll.<symbol>.argtypes`
  - `dll.<symbol>.restype`
  - the actual `dll.<symbol>(...)` call
  - numpy buffer dtype, ndim, shape, and ordering
- Fortran: `src/emses/emses_solver.f90`
  - `bind(c)` routine argument lists
  - `iso_c_binding` scalar kinds
  - `value` versus reference arguments
  - intent and array shape declarations

## Workflow

1. List C-facing entry points and wrapper signatures:
   ```bash
   rg -n "bind\\(c|argtypes|restype|dll\\." src/emses/emses_solver.f90 vdsolverf/emses/wrapper.py
   ```
2. For every `bind(c)` symbol, build a quick table with argument order, Fortran
   declaration, Python `ctypes` type, and the object passed at the call site.
3. Check these correspondences:

   | Python | Fortran |
   | --- | --- |
   | `c_int`, `c_double` | `integer(c_int), value`, `real(c_double), value` |
   | `POINTER(c_int)` with `byref(...)` | non-`value` `integer(c_int)` |
   | `c_char_p` buffer plus length | `character(1, c_char) :: path(*)` plus length |
   | `np.ctypeslib.ndpointer(dtype=np.float64, ndim=N)` | `real(c_double) :: array(...)` with N dimensions |

4. Verify numpy array layout assumptions. This code often declares Fortran
   dimensions in reverse order to match C-contiguous numpy buffers.
5. After any edit, run:
   ```bash
   .venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
   fpm test --target test_public_api
   ```

## Rules

- Update `argtypes`, conversion variables, call-site argument order, and return
  buffer slicing together.
- Set `restype` on the same symbol being called.
- When public Python names change, update `vdsolverf/emses/__init__.py`, docs,
  and README examples in the same change.

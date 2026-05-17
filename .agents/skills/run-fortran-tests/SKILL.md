---
name: run-fortran-tests
description: Run this repository's fpm Fortran tests, either the full auto-discovered test/test_*.f90 suite or a focused --target test, and summarize failures. Use after Fortran changes, Python/Fortran interface changes, probability or boundary logic edits, or before committing solver behavior changes.
---

# Run Fortran Tests

Use this skill to validate the Fortran side of `vdist-solver-fortran`.

## Workflow

1. Inspect available test programs when choosing a focused run:
   ```bash
   rg --files test
   ```
2. For broad validation, run:
   ```bash
   fpm test
   ```
3. For a focused check, run a single fpm target, for example:
   ```bash
   fpm test --target test_public_api
   ```
4. If output is long, summarize the important failure lines and rerun the
   failing target without truncation when the first error is unclear.

## Failure Triage

- Compile failure: identify the source file/module and the first compiler
  diagnostic.
- Runtime assertion failure: report the test program, label, actual value,
  expected value, and likely implementation file.
- Public API failure: inspect `test/test_public_api.f90` and the module public
  lists before changing exports.

## Test Edits

- Add new Fortran regression tests as separate `test/test_*.f90` programs when
  the behavior is independently runnable.
- Reuse `test/test_helpers.f90` helpers for numeric assertions.
- Keep focused physics tests small; broad simulation coverage belongs in a
  separate, intentionally slower target.

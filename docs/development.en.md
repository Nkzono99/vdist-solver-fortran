# Development guide

> Lang: [日本語](development.md) | **English**

Guide for making changes to `vdist-solver-fortran`. Complements the
agent-oriented notes in
[AGENTS.md](../AGENTS.md) / [CLAUDE.md](../CLAUDE.md).

## Environment

- `gfortran` and `fpm` must be on `PATH`. `make` builds call
  `fpm install --profile=release` then re-link the archive as a shared
  library.
- Python 3.12 via `.venv/` at the repository root:

  ```bash
  /usr/bin/python3.12 -m venv .venv
  .venv/bin/pip install --upgrade pip setuptools wheel
  .venv/bin/pip install f90nml emout numpy scipy tqdm
  ```

  Never use system `python3` — it is 3.6 on this cluster and does not
  support `typing.Literal`.

## Build cycle

```bash
fpm build              # Compile the static archive
fpm test               # Compile and run every test/*.f90 program
make                   # End-to-end: fpm install -> link shared lib into vdsolverf/
```

After touching Fortran, always run `fpm test` and then
`make` so the Python-facing `.so` stays in sync.

For Python-only changes a quick import check is usually enough:

```bash
.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
```

## Test organization

Tests live in [`test/`](../test/). `fpm` auto-detects any `.f90` program
under that tree, so adding a file is enough.

- [`test_helpers.f90`](../test/test_helpers.f90) — module with
  `assert_close`, `assert_close_vec`, `assert_equal_int`, `assert_true`.
- Unit tests: `test_particle`, `test_field`, `test_probabilities`,
  `test_maxwell_flux`,
  `test_photoelectron_raycast`.
- Integration: `test_solver_probability` exercises the full backtrace →
  collision → probability path with a hand-built simulator.
- Compile-time surface guards: `test_public_api` imports each public
  symbol that the C API / umbrella modules expose.

When adding a test, follow the existing pattern: a single `program` file
per concern, a top-of-program roll-call that prints `all tests passed.` on
success, each subroutine focused on one behaviour.

## Extending the solver

### New probability function

1. Declare `type, extends(t_Probability) :: t_MyProbability` and a
   factory function in a dedicated module (under `src/core/` if generic,
   `src/emses/` if EMSES-specific).
2. If the type owns heap state (e.g. a boundary list), add a branch to
   the `select type` inside
   [`destroy_simulator`](../src/emses/emses_simulator_builder.f90) to
   release it.
3. Wire it into `register_inner_boundary_probability` (for inner
   boundaries) or `add_probability_boundaries` (outer planes).
4. Add a `test_my_probability.f90` exercising the `at()` method.

### New namelist key

See [Architecture › Extension points](architecture.en.md#extension-points).
The two-sided update (`src/emses/namelist.f90` + `TMP_INP_KEYS`) is the
critical pair to keep in sync.

### New C entry point

See [Architecture › Extension points](architecture.en.md#extension-points)
and the
[`sync-wrapper-interface`](../.claude/skills/sync-wrapper-interface/SKILL.md)
skill which mechanically checks Fortran ↔ ctypes argument alignment.

## Release

Use the
[`release`](../.claude/skills/release/SKILL.md) skill: it bumps
`pyproject.toml` and `fpm.toml` together, writes a `CHANGELOG.md` entry,
builds a GitHub release body under `.release-notes/`, and creates the tag
after a successful `fpm test`.

## Checklist before opening a PR

- [ ] `fpm test` passes locally.
- [ ] `.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"` passes.
- [ ] Touched any `bind(c)` signature or `argtypes`? Rerun
      [`sync-wrapper-interface`](../.claude/skills/sync-wrapper-interface/SKILL.md).
- [ ] Added or removed a namelist key? Update `TMP_INP_KEYS` and
      [Namelist reference](namelist.en.md).
- [ ] Updated docs on both language variants when content changes user-
      visible behaviour (README, usage, namelist, physics).

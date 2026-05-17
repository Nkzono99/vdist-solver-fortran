---
name: release
description: Prepare a vdist-solver-fortran release by synchronizing versions, changelog text, release notes, validation, commit, tag, and optional GitHub release publication. Use when the user asks for a patch/minor/major release, a specific X.Y.Z version, release notes, or a GitHub release.
---

# Release

Use this skill for release preparation. Publishing tags or GitHub releases still
requires explicit user confirmation, even when local commits are allowed.

## Inputs

- Accept `patch`, `minor`, `major`, or an explicit `X.Y.Z`.
- If no version is given, inspect commits since the last tag and propose the
  smallest reasonable bump.

## Workflow

1. Inspect current version state and release range:
   ```bash
   rg -n '^version|^version:' pyproject.toml fpm.toml CITATION.cff
   git tag --sort=-creatordate
   git describe --tags --abbrev=0
   git log --oneline "$(git describe --tags --abbrev=0)..HEAD"
   ```
2. Ensure `pyproject.toml` and `fpm.toml` versions stay identical.
   `CITATION.cff` should be checked and updated when the project policy expects
   it to track releases.
3. Write or update `CHANGELOG.md` using Keep a Changelog-style headings:
   `Added`, `Changed`, `Deprecated`, `Removed`, `Fixed`, `Security`.
4. Write `.release-notes/vX.Y.Z.md` with a short highlight section and a
   compare link:
   `https://github.com/Nkzono99/vdist-solver-fortran/compare/vOLD...vX.Y.Z`
5. Validate before committing:
   ```bash
   fpm test
   .venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
   ```
6. Commit the release prep only after validation:
   ```bash
   git add pyproject.toml fpm.toml CHANGELOG.md .release-notes/vX.Y.Z.md
   git commit -m "Release vX.Y.Z"
   ```
7. Ask before creating tags, pushing, or publishing:
   ```bash
   git tag -a vX.Y.Z -m "Release vX.Y.Z"
   git push origin HEAD
   git push origin vX.Y.Z
   gh release create vX.Y.Z --title "vX.Y.Z" --notes-file .release-notes/vX.Y.Z.md
   ```

## Notes

- Do not copy raw commit logs into release notes; summarize user-visible
  changes.
- If there is no previous tag, ask how the user wants to define the initial
  release range.
- If validation fails, stop release work and report the blocker.

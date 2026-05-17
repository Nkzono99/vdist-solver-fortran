---
name: improve-agent-harness
description: Maintain this repository's agent guidance and local skills under .agents/, while using .claude/ as historical reference when helpful. Use when the user asks to improve AGENTS.md, add or revise repository skills, reduce stale agent instructions, align Claude/Codex guidance, or fix harness friction.
---

# Improve Agent Harness

Use this skill to keep the repository's agent-facing documentation useful,
current, and concise.

## Scope

- Primary Codex-facing guidance:
  - `AGENTS.md`
  - `.agents/skills/**/SKILL.md`
  - `.agents/skills/**/agents/openai.yaml`
- Historical or tool-specific reference:
  - `.claude/skills/**/SKILL.md`
  - `.claude/rules/**/*.md`
  - `CLAUDE.md`
  - `.claude/settings*.json`

## Workflow

1. Read current guidance before editing:
   ```bash
   sed -n '1,240p' AGENTS.md
   find .agents .claude -maxdepth 4 -type f
   ```
2. Identify stale assumptions:
   - old test layout names
   - outdated public API names
   - commands that ignore `.venv/bin/python`
   - duplicated prose that should live in one skill instead of `AGENTS.md`
   - Claude-only metadata copied into Codex skills
3. Keep `AGENTS.md` short. Put repeatable procedures into skills and leave
   `AGENTS.md` as the routing and project baseline.
4. For skill changes, keep YAML frontmatter to only `name` and `description`.
   Put UI metadata in `agents/openai.yaml`.
5. Validate changed skills:
   ```bash
   .venv/bin/python /home/b/b36291/.codex/skills/.system/skill-creator/scripts/quick_validate.py .agents/skills/<skill-name>
   ```
6. If `CLAUDE.md` intentionally mirrors `AGENTS.md`, update both. If they
   diverge by design, state that in the final summary.

## Guidelines

- Prefer specific trigger descriptions over long bodies.
- Do not add scripts, references, or assets directories unless the skill needs
  them.
- Remove stale instructions instead of annotating around them.
- Commit harness-only changes separately from solver implementation changes
  when the user allows commits.

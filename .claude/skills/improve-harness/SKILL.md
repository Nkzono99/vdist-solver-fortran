---
name: improve-harness
description: Claude Codeのハーネス(.claude配下の設定・skills・rules・commands)を自己改善する。`improve-harness`スキル自身も改善対象に含む。ハーネスが期待通り動かない／煩雑／抜け漏れがあるときに呼び出す。
user-invocable: true
---

# improve-harness Skill

Claude Codeのハーネスを自己観察し、反復的に改善するためのスキルです。自分自身(`.claude/skills/improve-harness/SKILL.md`)も改善対象です。

## 対象

`.claude/` 配下の以下を改善対象とします。

- `settings.json` / `settings.local.json` (権限, hooks, env)
- `skills/**/SKILL.md` (`improve-harness`自身を含む)
- `rules/**/*.md`
- `commands/**/*.md`
- ルート `CLAUDE.md` / `AGENTS.md`

## タスク

1. **現状把握**
   - `.claude/` 以下のファイルと `CLAUDE.md` / `AGENTS.md` を読み、現在の構成と意図を把握する。
   - ユーザから与えられたフィードバック・不満点・未対応の要望を整理する (指示があれば優先)。
2. **改善点の抽出**
   - 以下のチェックリストを当てる:
     - [ ] 冗長な`allow`ルールや期限切れのエントリがないか
     - [ ] `ask`/`deny`/`hooks`が最小限か (増えすぎていないか)
     - [ ] 頻繁に permission prompt が出るコマンドが `allow` に無いか
     - [ ] 各 SKILL の `description` が具体的で自動発火に十分か
     - [ ] `CLAUDE.md` / `AGENTS.md` が 200 行を超えていないか
     - [ ] `.venv` 等プロジェクト固有の作業環境が README / 開発ガイドに反映されているか
     - [ ] Python ラッパーと Fortran I/F に関するガイダンス (`rules/fortran-python-interop.md`) が最新か
     - [ ] 自身 (`improve-harness`) の指示が機能を網羅しているか — もし手順が抜けていれば追加する
3. **提案と適用**
   - 変更内容を 2-3 行のサマリで先に提示する。
   - ユーザが承認したら Edit / Write で適用し、差分を短く説明する。
   - 変更は小さく、関連ファイル単位でコミットに分ける (`git add <file>` で file ごと)。
4. **検証**
   - 可能なら `.venv/bin/python -c "from vdsolverf.emses import wrapper"` と `fpm test` を実行し回帰が無いことを確認。
   - `settings.json` に構文エラーが無いか `python3 -m json.tool .claude/settings.json` で検査。
5. **自己更新**
   - 今回のセッションで得た知見・不足していた手順を、この SKILL ファイルの「タスク」「注意事項」節に追記する。
   - 追記した行は 1 行ずつ簡潔に。冗長化したら古い内容を削除する。

## 注意事項

- `ask`/`deny`/`hooks` は最小限に保つ。ユーザの明示的な要望がない限り新規追加しない。
- 破壊的な設定変更 (既存の `allow` を大幅削除、hooks の追加など) は必ずユーザ確認を取る。
- `CLAUDE.md` と `AGENTS.md` は同一内容を維持する (片方だけ更新しない)。
- venv の Python は `.venv/bin/python`。システムの `python3` (3.6) では動かない機能がある。

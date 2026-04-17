---
name: release
description: 新バージョンをリリースする。pyproject.tomlとfpm.tomlのversionを同時に更新し、CHANGELOGとGitHubリリース本文を書き起こし、タグを打って（任意で）push/リリースを作成する。引数例：`patch` / `minor` / `major` / `1.5.0`。
argument-hint: [patch|minor|major|X.Y.Z]
user-invocable: true
---

# release Skill

このリポジトリの新バージョンをリリースする一連の手順をまとめるスキルです。バージョンが `pyproject.toml` と `fpm.toml` の二箇所に分散しているため同時更新が必須。リリースノートは CHANGELOG.md と GitHub リリース本文の両方に揃えます。

## 入力

- 引数 (`patch` / `minor` / `major` / `X.Y.Z`) が与えられれば採用。
- 省略された場合は直近のコミットログから推定案を提示し、ユーザ承認を得る。

## タスク

### 1. 現状把握

```bash
grep -n 'version' pyproject.toml fpm.toml CITATION.cff
git tag --sort=-creatordate | head -5
git log --oneline $(git describe --tags --abbrev=0)..HEAD
```

- `pyproject.toml` と `fpm.toml` の `version` が一致しているか確認。ズレていたら先に整える。
- `CITATION.cff` の `version` は過去ズレていた実績があるので念のため確認 (整合させるかはユーザに確認)。
- 直近タグから HEAD までのコミットを取得してリリース範囲を把握。

### 2. 新バージョンの決定

- `patch`: `X.Y.Z` → `X.Y.(Z+1)` (バグ修正のみ)
- `minor`: `X.Y.Z` → `X.(Y+1).0` (後方互換ありの機能追加)
- `major`: `X.Y.Z` → `(X+1).0.0` (破壊的変更)
- 明示バージョン指定はそのまま採用。ただし既存タグと衝突しないか `git tag | grep` で確認。

### 3. バージョンを同時更新

`pyproject.toml` と `fpm.toml` の `version = "OLD"` を Edit で差し替える。差し替え後に grep で両方 NEW に揃ったか確認する。

```bash
grep -n '^version' pyproject.toml fpm.toml
```

必要に応じて `CITATION.cff` の `version:` も更新 (ただしユーザ承認後)。

### 4. リリースノートを書く

2 つの成果物を用意する。

#### 4-1. `CHANGELOG.md` 追記

ファイルが無ければ新規作成。書式:

```markdown
# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/) and
the project adheres to [Semantic Versioning](https://semver.org/).

## [X.Y.Z] - YYYY-MM-DD

### Added
- ...

### Changed
- ...

### Fixed
- ...

### Removed
- ...
```

- 各項目は 1 行 1 変更。`file:line` や PR 番号があれば付記。
- 該当カテゴリが無ければその見出しごと省略。
- 直近タグから HEAD までの `git log --oneline` をベースに、ユーザ向けに意味のある粒度へ要約する (コミットそのままコピーしない)。

#### 4-2. GitHub リリース本文 (`.release-notes/vX.Y.Z.md`)

ディレクトリが無ければ作成。CHANGELOG の該当セクションを抜粋＋冒頭に 2-3 行のハイライトを付ける。

```markdown
## Highlights

- <ひとこと要約 1>
- <ひとこと要約 2>

<CHANGELOG の該当セクションを貼り付け>

**Full Changelog:** https://github.com/Nkzono99/vdist-solver-fortran/compare/vPREV...vX.Y.Z
```

### 5. 検証

```bash
fpm test 2>&1 | tail -20
.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
python3 -m json.tool .claude/settings.json > /dev/null
```

失敗したらリリース中止してユーザに報告。

### 6. コミット・タグ

ユーザに最終確認を取った上で:

```bash
git add pyproject.toml fpm.toml CHANGELOG.md .release-notes/vX.Y.Z.md
# (CITATION.cff を更新した場合はそれも add)
git commit -m "Release vX.Y.Z"
git tag -a vX.Y.Z -m "Release vX.Y.Z"
```

コミットメッセージには `Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>` を含める。

### 7. push / GitHub Release 作成 (ユーザ承認後)

```bash
git push origin main
git push origin vX.Y.Z
gh release create vX.Y.Z --title "vX.Y.Z" --notes-file .release-notes/vX.Y.Z.md
```

push と release 作成は破壊的・公開に影響するため**必ず**事前承認を取る。

## 注意事項

- **`pyproject.toml` と `fpm.toml` の version を両方同時に更新**する。片方だけだとパッケージと Fortran ライブラリのバージョンが食い違う。
- `CITATION.cff` の version は過去にずれていた実績があるため、整合させるかは毎回ユーザに確認する。
- タグ名は `vX.Y.Z` 形式で揃える (既存タグに合わせる)。
- リリース範囲のコミット列が長い場合、ユーザに「何を目玉としたいか」聞いてから Highlights を書くと迷子にならない。
- `gh release create` 実行前に本文プレビューをユーザに見せる。
- CHANGELOG のセクション順序は Added / Changed / Deprecated / Removed / Fixed / Security (Keep a Changelog 準拠)。

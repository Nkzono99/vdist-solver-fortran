# AGENTS.md

このリポジトリで作業するエージェント向けの入口です。詳細な反復手順は
`.agents/skills/` に置き、ここにはプロジェクト固有の前提と判断基準だけを
残します。

## プロジェクト概要

- プロジェクト名: `vdist-solver-fortran`
- 目的: Fortran で実装された速度分布ソルバを Python から利用するための
  ラッパー。
- 主要構成:
  - Fortran コア: `src/core/`, `src/emses/`, `src/utils/`
  - Python ラッパー: `vdsolverf/`
  - テスト: `test/test_*.f90`
  - ビルド: `fpm.toml`, `Makefile`, `setup.py`
  - エージェント用スキル: `.agents/skills/`

## スキルの使い分け

- `$build-lib`: Fortran 変更後に共有ライブラリを再ビルドし、
  `vdsolverf/libvdist-solver-fortran.*` を更新する。
- `$run-fortran-tests`: `fpm test` または `fpm test --target ...` で
  Fortran 回帰を確認する。
- `$sync-wrapper-interface`: `vdsolverf/emses/wrapper.py` と
  `src/emses/emses_solver.f90` の `bind(c)` ABI を突き合わせる。
- `$release`: バージョン更新、CHANGELOG、リリースノート、タグ、GitHub
  Release を準備する。
- `$improve-agent-harness`: `AGENTS.md` と `.agents/skills/` 自体を更新する。

`.claude/` 配下の rules/skills は過去の Claude 用ハーネスです。参考には
できますが、Codex 向けの一次情報は `AGENTS.md` と `.agents/skills/` です。

## 重要ファイル

- `vdsolverf/emses/wrapper.py`
  - `ctypes` で共有ライブラリをロードし、Fortran の C API を呼び出す。
- `src/emses/emses_solver.f90`
  - Python から入る `bind(c)` エントリポイントを定義する。
- `src/emses/emses_simulator_builder.f90`
  - EMSES 入力からシミュレータ、境界、確率関数を構築する。
- `vdsolverf/emses/tmpolary_input.py`
  - `emout` から一時 namelist (`plasma-vdsolverf.inp`) を生成する。
- `vdsolverf/core/particles.py`, `vdsolverf/core/phase_grid.py`
  - Python 側の粒子・位相空間グリッドのデータモデル。

## 開発環境

- Python はリポジトリ直下の `.venv/` を使う。
  - `.venv/bin/python`
  - `.venv/bin/pip`
- システムの `python3` は古い環境の可能性があるため、Python ラッパー確認
  には使わない。
- Python smoke test:
  ```bash
  .venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
  ```
- Fortran test:
  ```bash
  fpm test
  ```

## ビルドの前提

- 共有ライブラリ名は OS ごとに固定:
  - Linux: `vdsolverf/libvdist-solver-fortran.so`
  - macOS: `vdsolverf/libvdist-solver-fortran.dylib`
  - Windows: `vdsolverf/libvdist-solver-fortran.dll`
- `Makefile` は `fpm install` 後に共有ライブラリを生成し、`vdsolverf/`
  へコピーする。
- `setup.py` の custom build は `make` を実行するため、Python パッケージ
  のビルドにも Fortran ツールチェーンが必要。

## 変更時のルール

1. Python と Fortran の I/F 整合を最優先する。
   - `wrapper.py` の `argtypes`, `restype`, dtype, 配列次元、呼び出し順を
     変える場合は、対応する `bind(c)` 側も同時に確認する。
2. 公開 API の破壊的変更は、ユーザが明示的に許可した場合だけ行う。
   - 対象例: `vdsolverf/core/__init__.py`, `vdsolverf/emses/__init__.py`,
     `get_backtrace`, `get_backtraces`, `get_probabilities`。
3. プラットフォーム分岐を壊さない。
   - `system` 判定 (`linux`/`darwin`/`windows`) と `CDLL`/`WinDLL` の
     ロード方式を維持する。
4. EMSES 入力の namelist キー変更は慎重に扱う。
   - `TMP_INP_KEYS`, `convert_from_geotype`, `src/emses/namelist.f90`,
     docs の整合を取る。
5. 既存のスペル揺れは単独で直さない。
   - 例: `TempolaryInput`, `tmpolary_input.py`。改名は影響範囲が広い。

## 推奨チェック

- Python だけの変更:
  - `.venv/bin/python -c "from vdsolverf.core import Particle, PhaseGrid; from vdsolverf.emses import wrapper; print('ok')"`
- Fortran または ABI の変更:
  - `$sync-wrapper-interface`
  - `fpm test`
  - `$build-lib`
  - Python smoke test
- 外部ユーザー向け挙動の変更:
  - `README.md`, `README.en.md`, `docs/usage*.md` の Usage 例を更新する。
  - `ispec` の意味など README 記載済み仕様と実装の不一致を作らない。

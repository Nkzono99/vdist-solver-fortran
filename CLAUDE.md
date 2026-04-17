# AGENTS.md

このリポジトリで作業するエージェント向けのガイドです。`README.md` と実装コード（Python + Fortran）を前提にしています。

## プロジェクト概要

- プロジェクト名: `vdist-solver-fortran`
- 目的: Fortran で実装された速度分布ソルバを Python から利用するためのラッパー。
- コア構成:
  - Fortran本体: `src/`
  - Pythonラッパー: `vdsolverf/`
  - ビルド: `fpm.toml` + `Makefile` + `setup.py`

## ディレクトリ構成（重要箇所）

- `src/core/*.f90` / `src/emses/*.f90`
  - 物理計算・トレース計算のFortran実装。
- `vdsolverf/emses/wrapper.py`
  - `ctypes` 経由で共有ライブラリをロードし、Fortran関数を呼び出す主要ラッパー。
- `vdsolverf/emses/tmpolary_input.py`
  - `emout` の入力情報から一時namelist (`plasma-vdsolverf.inp`) を生成。
- `vdsolverf/core/particles.py`, `vdsolverf/core/phase_grid.py`
  - Python側のデータモデル（粒子、位相空間グリッド）。

## ビルドと実行に関する注意

- 共有ライブラリ名は OS ごとに固定:
  - Linux: `vdsolverf/libvdist-solver-fortran.so`
  - macOS: `vdsolverf/libvdist-solver-fortran.dylib`
  - Windows: `vdsolverf/libvdist-solver-fortran.dll`
- `Makefile` は `fpm install` 後に共有ライブラリを生成し、`vdsolverf/` へコピーする。
- `setup.py` の custom build は `make` を実行するため、Pythonパッケージのビルドには Fortran ツールチェーンが必要。

## 開発環境 (Python)

- Python 3.12 の venv をリポジトリ直下 `.venv/` に配置する。
  - 初期化: `/usr/bin/python3.12 -m venv .venv`
  - 依存導入: `.venv/bin/pip install -e .` もしくは `.venv/bin/pip install f90nml emout numpy scipy tqdm`
- Python コマンドは原則 `.venv/bin/python` / `.venv/bin/pip` を使用する（システムの `python3` は 3.6 で `Literal` 未対応のため利用不可）。
- Fortran テスト: `fpm test`
- Python smoke test: `.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"`

## 変更時のガイドライン

1. **PythonとFortranのI/F整合を最優先**
   - `vdsolverf/emses/wrapper.py` の `argtypes`/`dtype`/配列次元を変更する場合、対応するFortran側の引数定義と同時に確認する。
2. **公開APIの互換性を維持**
   - `vdsolverf/core/__init__.py` と `vdsolverf/emses/__init__.py` の公開関数・クラスは既存利用者への影響が大きい。破壊的変更は避ける。
3. **プラットフォーム分岐を壊さない**
   - `wrapper.py` の `system` 判定（linux/darwin/windows）とロード方式（`CDLL`/`WinDLL`）を維持する。
4. **EMSES入力の取り扱いに注意**
   - `tmpolary_input.py` で扱う namelist キーは計算結果に直結。キーの追加・削除時は `TMP_INP_KEYS` と `convert_from_geotype` の整合を取る。
5. **命名・スタイルは既存に合わせる**
   - 既存実装にはスペル揺れ（例: `TempolaryInput`）がある。リネームは影響範囲が広いため、単独での改名は避ける。

## 推奨チェック項目（変更後）

- Python側だけの変更でも、少なくとも次を確認:
  - import が壊れていないこと（`vdsolverf.core`, `vdsolverf.emses`）
  - `wrapper.py` のプラットフォーム別ライブラリパスが有効なままであること
- Fortran変更時は、`Makefile` / `fpm.toml` のビルドフロー影響を確認。

## README整合のルール

- 外部ユーザー向け挙動を変更した場合は `README.md` の Usage 例も更新する。
- `ispec` の意味（0: electron, 1: ion, 2: photoelectron）など、README記載済みの仕様と実装の不一致を作らない。

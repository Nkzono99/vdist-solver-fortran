# 開発ガイド

> Lang: **日本語** | [English](development.en.md)

`vdist-solver-fortran` に変更を加えるためのガイドです。
[AGENTS.md](../AGENTS.md) / [CLAUDE.md](../CLAUDE.md) の
エージェント向けメモを補完します。

## 環境

- `gfortran` と `fpm` が `PATH` 上にあること。`make` は
  `fpm install --profile=release` を内部で呼び出したあと、アーカイブを
  共有ライブラリとして再リンクします。
- リポジトリ直下の `.venv/` に Python 3.12:

  ```bash
  /usr/bin/python3.12 -m venv .venv
  .venv/bin/pip install --upgrade pip setuptools wheel
  .venv/bin/pip install f90nml emout numpy scipy tqdm
  ```

  システムの `python3` は 3.6 で `typing.Literal` が使えないため
  利用不可。

## ビルドサイクル

```bash
fpm build              # 静的アーカイブをコンパイル
fpm test               # test/*.f90 の全プログラムをビルドして実行
make                   # エンドツーエンド: fpm install → 共有ライブラリを vdsolverf/ に配置
```

Fortran を触ったら必ず `fpm test` を回してから `make` を走らせ、Python
から参照する `.so` を同期させてください。

Python 側のみの変更なら import チェックで十分なことが多いです:

```bash
.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
```

## テスト構成

テストは [`test/`](../test/) に置きます。`fpm` がサブツリー内の `.f90`
プログラムを自動検出するので、ファイルを追加するだけで組み込まれます。

- [`test_helpers.f90`](../test/test_helpers.f90) — `assert_close`、
  `assert_close_vec`、`assert_equal_int`、`assert_true` を提供する
  モジュール。
- ユニット: `test_particle`、`test_field`、`test_probabilities`、
  `test_maxwell_flux`、
  `test_photoelectron_raycast`。
- 統合: `test_solver_probability` は手組みシミュレータでバックトレース →
  衝突 → 確率評価の一連を通す。
- コンパイル時 API サーフェスガード: `test_public_api` が C API /
  アンブレラモジュールの公開シンボルをすべて import。

テスト追加時は既存パターンに合わせてください: 関心ごとに 1 program
ファイル、先頭でサブルーチンを列挙、成功時に `all tests passed.` を
出力、各サブルーチンで 1 振る舞いだけを検証。

## ソルバの拡張

### 新しい確率関数の追加

1. 専用モジュール (汎用なら `src/core/`、EMSES 固有なら `src/emses/`) に
   `type, extends(t_Probability) :: t_MyProbability` とファクトリを定義。
2. heap 上の状態 (境界リストなど) を持つなら
   [`destroy_simulator`](../src/emses/emses_simulator_builder.f90) の
   `select type` 分岐を追加して解放。
3. 内部境界用なら `register_inner_boundary_probability` に、外側平面用
   なら `add_probability_boundaries` に配線。
4. `test_my_probability.f90` で `at()` の挙動を検証。

### 新しい namelist キー

[アーキテクチャ › 拡張ポイント](architecture.md#拡張ポイント) を参照。
`src/emses/namelist.f90` と `TMP_INP_KEYS` のペア更新が最重要。

### 新しい C 入口

[アーキテクチャ › 拡張ポイント](architecture.md#拡張ポイント) と、
Fortran ↔ ctypes 引数整合を機械的にチェックする
[`sync-wrapper-interface`](../.claude/skills/sync-wrapper-interface/SKILL.md)
スキルを活用。

## リリース

[`release`](../.claude/skills/release/SKILL.md) スキルを使用: `pyproject.toml`
と `fpm.toml` のバージョンを同時に上げ、`CHANGELOG.md` に追記し、
`.release-notes/` に GitHub リリース本文を書き、`fpm test` が通ったあと
にタグを打ちます。

## PR 作成前のチェックリスト

- [ ] `fpm test` がローカルで通過。
- [ ] `.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"` が通過。
- [ ] `bind(c)` シグネチャや `argtypes` を触った場合は
      [`sync-wrapper-interface`](../.claude/skills/sync-wrapper-interface/SKILL.md)
      を再実行。
- [ ] namelist キーを追加・削除した場合は `TMP_INP_KEYS` と
      [Namelist リファレンス](namelist.md) を更新。
- [ ] ユーザから見える挙動を変更した場合は、日本語・英語の両言語版
      (README、usage、namelist、physics) を更新。

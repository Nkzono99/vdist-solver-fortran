# アーキテクチャ

> Lang: **日本語** | [English](architecture.en.md)

`vdist-solver-fortran` は Fortran の数値計算コアを `ctypes` ブリッジで
Python に公開した構成です。本ドキュメントはリポジトリレイアウトと、
実行時に実際に走るモジュール群の対応を示します。

## トップレベルのレイアウト

```
src/                  Fortran ソース (fpm でビルド)
  core/               汎用物理とソルバインフラ
  emses/              EMSES 固有の境界構築と C API
  utils/              小さな共通ヘルパ
  vdsolverf.f90       アンブレラモジュール (m_vdsolverf)。C API を再公開
vdsolverf/            Python パッケージ
  core/               データクラス (Particle, DustParticle, PhaseGrid)
  emses/              ctypes ラッパ、一時入力ビルダ、geotype ヘルパ
fpm.toml              fpm ビルド設定 (共有ライブラリ、テスト自動検出)
Makefile              `fpm install` と OS 別の共有ライブラリリンクをラップ
setup.py              Python ビルド。内部で `make` を呼ぶ
test/                 fpm テストプログラムと m_test_helpers モジュール
docs/                 本ディレクトリ
```

## Fortran モジュール依存関係

```
m_vdsolverf                           (src/vdsolverf.f90)
  └── m_emses_solver                  (src/emses/emses_solver.f90)      -- C API (bind(c))
        └── m_emses_simulator_builder (src/emses/emses_simulator_builder.f90)
              ├── m_allcom            (src/emses/allcom.f90)            -- namelist グローバル変数
              ├── m_namelist          (src/emses/namelist.f90)          -- ファイル入出力
              ├── m_emses_boundaries  (src/emses/collision/*.F90)       -- 境界構築
              ├── m_photoelectron_raycast (src/emses/photoelectron_raycast.f90)
              └── m_vdsolverf_core    (src/core/vdsolverf_core.f90)     -- core の集約再公開
                    ├── m_particle
                    ├── m_field
                    ├── m_probabilities
                    ├── m_simulator
                    ├── m_dust_charge_simulator
                    └── m_solver
```

C シンボルを公開するのは `m_emses_solver` だけです。ビルダモジュールが
シミュレータの構築 (raycast 確率の配線を含む) と `destroy_simulator` に
よるクリーンアップを担当します。C API から辿れるすべてはアンブレラを
通るので、Python 利用者から見える入口は 3 つ (`get_backtraces`,
`get_probabilities`, `get_backtrace_dust`) だけです。

## Python パッケージ

```
vdsolverf.core         Particle, DustParticle, PhaseGrid データクラス
vdsolverf.emses.wrapper
  _load_dll(...)       プラットフォームに応じた共有ライブラリを解決
  get_backtrace(...)   単一粒子用のコンビニエンスラッパ
  get_backtraces(...)  多粒子 ctypes 呼び出し
  get_probabilities(...)
  get_dust_backtrace(...)
  create_relocated_ebvalues / create_relocated_current_values
                       emout から EB / 電流場配列を組み立てる
vdsolverf.emses.tmpolary_input
  TempolaryInput       emout から最小限の plasma-vdsolverf.inp を書き出して
                       終了時に削除するコンテキストマネージャ
  TMP_INP_KEYS         Fortran に渡す namelist キーのホワイトリスト
vdsolverf.emses.geotype
  geotype プリミティブを boundary_type / boundary_shape タプルに変換
```

## 境界を越える流れ

1. Python が EB 場 (ダストモードでは電流場も) を `bind(c)` サブルーチン
   で宣言された形状に合わせて `numpy` 配列に詰めます。引数アライメント
   規則は
   [`.claude/rules/fortran-python-interop.md`](../.claude/rules/fortran-python-interop.md)
   を参照。
2. `TempolaryInput` がフィルタ済み namelist を
   `data.directory / plasma-vdsolverf.inp` に書き出します。
3. `_load_dll` が `libvdist-solver-fortran.so` / `.dylib` / `.dll` を解決。
4. ctypes 呼び出しで 3 つの `bind(c)` 入口のいずれかに入り、シミュレータ
   を構築、バックトレースループを回し、終了時に `destroy_simulator` を
   呼びます。
5. `TempolaryInput.__exit__` で一時ファイルが削除されます。

## 拡張ポイント

- **新しい確率関数**: `src/core/` または `src/emses/` の新モジュールで
  `type, extends(t_Probability)` を定義し、内部境界用なら
  `register_inner_boundary_probability`、外側平面用なら
  `add_probability_boundaries` で allocate。heap 上の状態を持つなら
  `destroy_simulator` の `select type` 分岐も追加。
- **新しい namelist キー**: `src/emses/namelist.f90` の該当グループに
  追記し、`vdsolverf/emses/tmpolary_input.py` の `TMP_INP_KEYS` にも
  追加、[Namelist リファレンス](namelist.md) を更新。
- **新しい C 入口**: `src/emses/emses_solver.f90` に `bind(c)` サブルー
  チンを追加し、`wrapper.py` でシグネチャ (`argtypes`、`restype`、
  呼び出しサイト) をミラーし、
  [`sync-wrapper-interface`](../.claude/skills/sync-wrapper-interface/SKILL.md)
  で整合を確認。

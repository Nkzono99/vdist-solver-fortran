# ドキュメント

> Lang: **日本語** | [English](README.en.md)

`vdist-solver-fortran` のリファレンス資料です。プロジェクトの概要と
インストール手順は [トップレベルの README](../README.md) を参照してください。

## 目次

| ドキュメント | 内容 |
|---|---|
| [使用方法](usage.md) | Python 側のレシピ集 (単一バックトレース、多粒子バックトレース、位相空間確率ソルバ) |
| [速度範囲自動推定](autorange.md) | `estimate_velocity_range_map` によるセルごとの速度範囲推定、validation、diagnostics、性能メモ |
| [Namelist リファレンス](namelist.md) | サポートする EMSES `plasma.inp` グループ (`/ptcond/`, `/emissn/`) と raycast 光電子用パラメータ |
| [物理モデル](physics.md) | リウビル定理、Maxwellian 表面放出、raycast 光電子モデル |
| [アーキテクチャ](architecture.md) | リポジトリ構成、Fortran/Python の境界、モジュールの責務 |
| [開発ガイド](development.md) | ビルド、テスト、拡張方法 (`fpm`、`.venv`、確率関数の追加、CI チェックリスト) |

## 読者別の入口

- **EMSES の後処理を回す研究者** — [使用方法](usage.md)、
  [速度範囲自動推定](autorange.md) と
  [Namelist リファレンス](namelist.md) から。
- **Fortran / Python を編集する貢献者** — [アーキテクチャ](architecture.md) を
  ざっと確認してから [開発ガイド](development.md)。
- **物理モデルを確認したい読者** — [物理モデル](physics.md) が各確率関数の
  計算内容を導出しています。

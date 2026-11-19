# EMSES との同期

> Lang: **日本語** | [English](emses_sync.en.md)

本ソルバは MPIEMSES3D の出力を再生するため、EMSES の形状・粒子源・運動の式を
EMSES と同じにしておく必要があります。EMSES のコードはライブラリとして共有せず、
次の 2 通りで追従します。

- **形状はそのまま写す。** 内部境界の生成コードは EMSES の
  `src/physics/collision/` のファイルを byte 単位で同じ内容のまま置きます。
- **粒子源と運動は本ソルバの側で書く。** 配列やデータの持ち方が違うため、
  EMSES の式を本ソルバのコードで再実装し、テストで EMSES の定義と照合します。

最後に照合した EMSES の commit は
[`scripts/check-emses-sync.sh`](../scripts/check-emses-sync.sh) の
`SYNCED_COMMIT` が唯一の記録です。

## 写しているファイル

| vdist | MPIEMSES3D |
|---|---|
| `src/emses/collision/objects.F90` | `src/physics/collision/objects.F90` |
| `src/emses/collision/surfaces.F90` | `src/physics/collision/surfaces.F90` |
| `src/emses/collision/boundary_types.F90` | `src/physics/collision/boundary_types.F90` |
| `src/emses/collision/boundary_helpers.F90` | `src/physics/collision/boundary_helpers.F90` |
| `src/emses/collision/collision_accelerator.F90` | `src/physics/collision/collision_accelerator.F90` |

これらのファイルは編集しません。修正が必要なら EMSES で直してから写します。
本ソルバ固有なのは次の 2 つです。

- `src/emses/collision/collision_boundary_input_config.F90`: EMSES の
  `src/main/config/collision_boundary_input_config.F90` と同じ型を、
  `m_allcom` の namelist 変数から埋めます。
- `src/emses/collision/emses_boundaries.F90`: 確率関数の tag を付けます。

gfortran の既定（1 行 132 文字）でそのままコンパイルできるよう、
EMSES 側のこれらのファイルも 132 文字以内に保ちます。

## 再実装している挙動

| 挙動 | 本ソルバ | EMSES の参照元 |
|---|---|---|
| 外側境界の注入ドリフト・磁場座標の熱分布 | `vdri_vector`、`vth_vector`、`thermal_basis` (`src/emses/allcom.f90`) | `initialization_simulation_setup.F90`（`calculate_drift_velocities`）、`domain_boundary_source.F90`、`m_thermal_frame.F90` |
| 表面源の平均・熱幅、パッチ電流と入射角 | `surface_source_parameters`、`emission_vdri_vector`、`patch_flux` | `surface_emission_source.F90`、`direct_surface_emitter.F90` |
| raycast 光電子の流束・法線・接線・太陽方向 | `src/emses/photoelectron_raycast.f90`、`resolve_sun_direction` | `raycast_surface_emitter.F90`、`raycast_sampler.F90` |
| 片側フラックスモーメント | `src/utils/maxwell_flux_erf.f90` | `m_velocity_distribution.F90`（`velocity_flux_rate`） |
| Boris 更新と格子補間 | `src/core/simulator.f90`、`src/core/field.f90` | `particle_solver.F90`、`m_trilinear_interp.F90` |
| namelist のキー | `src/emses/namelist.f90`、`vdsolverf/emses/tmpolary_input.py` | `namelist_reader.F90` |

式の内容は [物理モデル](physics.md#対応するソース分布) にまとめています。
本ソルバが再現できない EMSES の設定は、黙って無視せず `ValueError` で拒否します。
Python 側（`vdsolverf/emses/source_support.py`）は EMSES の元の namelist を、
Fortran 側（`check_emses_input`）は本ソルバ用の一時 namelist を調べます。

## 追従の手順

```sh
scripts/check-emses-sync.sh /path/to/MPIEMSES3D origin/main
```

1. 「mirrored files」に `DIFFERS` が出たら、EMSES のファイルを写して `fpm test` を
   実行します。新しい形状が増えた場合は `check_geometry`
   (`src/emses/emses_input_validation.f90`) と namelist のキーも更新します。
2. 「EMSES commits since ...」に出た commit を読み、上の表の挙動に関わる変更を
   本ソルバへ移します。挙動を変えたら対応するテストを追加します。
3. すべて確認したら `SYNCED_COMMIT` を照合した EMSES の commit に更新します。

KUDPC ではビルドとテストを計算ノードで実行してください。

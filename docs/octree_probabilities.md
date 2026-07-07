# Octree 速度空間確率ソルバ

> Lang: **日本語** | [English](octree_probabilities.en.md)

`get_probabilities_octree` は、指定した空間点ごとに速度空間だけを adaptive
octree で細分化し、各サンプル点の `get_probabilities` 相当の確率を Fortran
側で評価する API です。全空間に一様な dense 6D 配列を作らず、確率が変化する
速度領域を重点的に評価したい場合に使います。

## 基本例

```python
import emout
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_probabilities_octree

data = emout.Emout("my-run")

phase_grid = PhaseGrid(
    x=(120, 180, 16),
    y=64,
    z=(300, 420, 16),
    vx=(-8.0e6, 8.0e6, 2),
    vy=(-4.0e6, 4.0e6, 2),
    vz=(-8.0e6, 8.0e6, 2),
)

result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    phase_grid=phase_grid,
    dt=data.inp.dt,
    max_step=30_000,
    use_adaptive_dt=False,
    scout_bins=(11, 9, 11),
    max_depth=5,
    max_samples_per_cell=80_000,
    max_leaves_per_cell=4096,
    n_threads=32,
)
```

`PhaseGrid` の `x/y/z` が評価する空間点、`vx/vy/vz` が各空間点の root
velocity box になります。octree では `vx/vy/vz` の bin 数は使わず、範囲だけを
使います。

## 任意の空間点とセル別速度範囲

`PhaseGrid` の代わりに、空間点配列とセル別の速度範囲を渡せます。

```python
positions = [
    [120.5, 64.5, 320.5],
    [121.5, 64.5, 320.5],
]

velocity_bounds = [
    [-8e6, 8e6, -4e6, 4e6, -8e6, 8e6],
    [-6e6, 6e6, -3e6, 3e6, -9e6, 7e6],
]

result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    position=positions,
    velocity_bounds=velocity_bounds,
    dt=data.inp.dt,
    max_step=30_000,
)
```

`velocity_bounds` は `(6,)`, `(3, 2)`, `(nspatial, 6)`,
`(nspatial, 3, 2)` を受け取れます。`(6,)` または `(3, 2)` を複数空間点に
渡した場合は、同じ velocity box が全点に使われます。

## 戻り値

戻り値は `VelocityOctreeResult` です。dense 6D 配列ではなく、compact された
サンプル列と octree box 列を持ちます。

| 属性 | 形状 | 意味 |
|---|---:|---|
| `spatial_points` | `(nspatial, 3)` | 評価した空間点 |
| `velocities` | `(nsample, 3)` | 評価した速度サンプル |
| `probabilities` | `(nsample,)` | サンプルごとの確率。未到達は `nan` |
| `spatial_index` | `(nsample,)` | 各サンプルが属する空間点 index |
| `leaf_bounds` | `(nleaf, 6)` | 評価した octree box の速度範囲 |
| `leaf_value_min/max` | `(nleaf,)` | box 内サンプル確率の min/max |
| `leaf_depth` | `(nleaf,)` | octree depth |
| `leaf_sample_start/count` | `(nleaf,)` | box に対応するサンプル範囲 |
| `status` | `(nspatial,)` | 空間点ごとの status |
| `sample_count`, `leaf_count` | `(nspatial,)` | 空間点ごとの出力数 |

status は以下です。

| 値 | 意味 |
|---:|---|
| `0` | OK |
| `1` | root box 内に有効確率が見つからない |
| `2` | `max_samples_per_cell` に到達 |
| `3` | edge validation により root box を拡張しきった |
| `4` | `max_leaves_per_cell` に到達 |
| `5` | velocity bounds が不正 |

## 可視化

サンプル点を空間点ごとに取り出せます。

```python
import numpy as np

i = 0
mask = result.spatial_index == i
v = result.velocities[mask]
p = result.probabilities[mask]

valid = np.isfinite(p)
# 例: vx-vz 平面へ scatter
ax.scatter(v[valid, 0], v[valid, 2], c=p[valid], s=2)
```

octree box の射影を見たい場合は `leaf_spatial_index` と `leaf_bounds` を使います。

```python
leaf_mask = result.leaf_spatial_index == i
boxes = result.leaf_bounds[leaf_mask]
value_range = result.leaf_value_max[leaf_mask] - result.leaf_value_min[leaf_mask]
```

## 探索パラメータ

| 引数 | 目安 |
|---|---|
| `scout_bins` | root box の初期サンプル数。狭い複数ローブを拾いたい場合は大きくする。 |
| `max_depth` | octree の最大細分化深さ。深くすると高解像度になるが指数的に増える。 |
| `refine_threshold_rel` | box 内の `pmax - pmin` が `pmax` に対してこの値以上なら split。 |
| `edge_threshold_rel` | root box の端に有意確率があるかを判定する閾値。 |
| `expand_factor` / `max_expansions` | 端に確率が残る root box を拡張する設定。 |
| `max_samples_per_cell` | 空間点ごとのサンプル容量。status `2` が多い場合は増やす。 |
| `max_leaves_per_cell` | 空間点ごとの octree box 容量。status `4` が多い場合は増やす。 |
| `n_threads` | Fortran OpenMP thread 数。空間点方向に並列化する。 |

全 box が depth `D` まで split される最悪ケースでは、評価 box 数は
`(8 ** (D + 1) - 1) / 7` です。root 以外の box は現在 `3x3x3` サンプルで評価
するため、容量は余裕を持って設定してください。

## 注意点

- この API は確率を直接評価するため、`estimate_velocity_range_map` のような
  envelope 誤差の伝搬はありません。
- 一方で、初期 `scout_bins` が狭いローブを完全に外すと、そのローブは split
  されません。複数の小さな Maxwellian ローブが想定される場合は
  `scout_bins` を増やすか、手動で velocity box を分けて複数回評価します。
- `leaf_bounds` は最終 terminal leaf だけでなく、評価済みの中間 box も含みます。
- 返るのは sparse なサンプル/box 表現です。保存したい場合は
  `np.savez` などで `result.velocities`, `result.probabilities`,
  `result.leaf_bounds` を保存してください。

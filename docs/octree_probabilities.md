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

`PhaseGrid` の3要素タプルは `(start, end, count)` です。例えば
`x=(120, 180, 16)` は 120 から 180 までの 16 点を `np.linspace` で作る指定で、
刻み幅 16 ではありません。

上の `max_depth=5`, `max_samples_per_cell=80_000`,
`max_leaves_per_cell=4096` は、全 box を depth 5 まで完全に評価する容量ではなく、
実務用の budget-limited 例です。status `2` または `4` が多い場合は容量を
増やすか、velocity box を分割して評価してください。

## アルゴリズム概要

octree は空間方向を細分化しません。各空間点に対して、指定された
`vx/vy/vz` の root velocity box だけを 3D octree として探索します。
各速度サンプルでは、通常の `get_probabilities` と同じく粒子をバックトレースし、
境界到達時の確率関数を評価します。

1. `velocity_bounds` から root box を作る。
2. root box 全体を `scout_bins` の規則格子で評価する。
3. root box の境界面に有意な確率が残る場合、`expand_factor` で box を拡張する。
   これは `max_expansions` 回まで繰り返す。
4. root box を queue に入れ、box ごとにサンプルする。
   root では `scout_bins`、子 box では現在 `3x3x3` 点を評価する。
5. box 内の `pmax > 0` かつ `(pmax - pmin) / pmax >= refine_threshold_rel`
   なら、`max_depth` まで 8 子 box に split する。
6. `max_samples_per_cell` または `max_leaves_per_cell` に到達した空間点は、
   それまでの partial result と status を返す。

このため、狭い速度ローブや複数の Maxwellian-like lobe がある場合でも、
root の `scout_bins` がその近傍を一度でも拾えば局所的に細分化されます。
逆に、root scout grid の隙間に完全に入るほど狭いローブは見逃され得ます。
その場合は `scout_bins` を増やすか、速度範囲を複数の box に分けて
`get_probabilities_octree` を複数回実行してください。

## KUDPC での実行例

`get_probabilities_octree` は Fortran/OpenMP で重いバックトレースを実行するため、
KUDPC の login node では直接長時間実行せず、`tssrun` または batch job 内で
実行してください。`n_threads` は確保した CPU 数に合わせます。

```bash
tssrun -p gr20001g --rsc p=1:t=32:c=32:m=64G -t 6:00 \
  bash -lc 'export OMP_NUM_THREADS=32; .venv/bin/python run_octree.py'
```

Python 側では `n_threads` を明示するか、`OMP_NUM_THREADS` に任せます。

```python
result = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    phase_grid=phase_grid,
    dt=data.inp.dt,
    max_step=30_000,
    n_threads=32,
)
```

大きなケースでは、まず少数の空間点と浅い `max_depth` で `status` と
`sample_count` を確認してから、空間点数・depth・capacity を増やしてください。

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
| `leaf_bounds` | `(nleaf, 6)` | 評価した octree box の速度範囲。分割済みの中間 box も含む |
| `leaf_value_min/max` | `(nleaf,)` | box 内サンプル確率の min/max |
| `leaf_depth` | `(nleaf,)` | octree depth |
| `leaf_sample_start/count` | `(nleaf,)` | compact 後の `velocities/probabilities` に対する box ごとのサンプル範囲 |
| `status` | `(nspatial,)` | 空間点ごとの status |
| `sample_count`, `leaf_count` | `(nspatial,)` | 空間点ごとの出力数 |
| `metadata` | `dict` | 実行時パラメータ、`actual_sample_count`、`actual_leaf_count` |

status は以下です。

| 値 | 意味 |
|---:|---|
| `0` | OK |
| `1` | root box 内に正の確率 signal が見つからない |
| `2` | `max_samples_per_cell` に到達 |
| `3` | root box が少なくとも1回拡張された。端の signal が残っているとは限らない |
| `4` | `max_leaves_per_cell` に到達 |
| `5` | velocity bounds が不正 |

非ゼロ status の空間点にも scout sample や部分的な sample が返ることがあります。
可視化や後段処理では、先に status を確認してください。

```python
import numpy as np

bad = np.flatnonzero(result.status != 0)
if bad.size:
    print("non-OK octree points:", bad[:10], "status:", result.status[bad[:10]])
```

## 可視化

サンプル点を空間点ごとに取り出せます。

```python
import numpy as np

i = 0
if result.status[i] != 0:
    print("warning: non-OK octree status", result.status[i])

v = result.velocities_for_spatial(i)
p = result.probabilities_for_spatial(i)

valid = np.isfinite(p)
# 例: vx-vz 平面へ scatter
ax.scatter(v[valid, 0], v[valid, 2], c=p[valid], s=2)
```

octree box の射影を見たい場合は `leaf_spatial_index` と `leaf_bounds` を使います。
ここでの `leaf_*` は出力互換のための名前で、terminal leaf だけでなく、split 判定
前に評価された中間 box も含みます。

```python
leaf_mask = result.leaf_mask(i)
boxes = result.leaf_bounds[leaf_mask]
value_range = result.leaf_value_max[leaf_mask] - result.leaf_value_min[leaf_mask]
```

ある leaf に属するサンプルは、compact 済み配列へ直接 slice できます。

```python
leaf_index = np.flatnonzero(leaf_mask)[0]
sample_slice = result.leaf_sample_slice(leaf_index)
leaf_velocities = result.velocities[sample_slice]
leaf_probabilities = result.probabilities[sample_slice]
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

## status 別の対処

| status | 典型的な状況 | 対処 |
|---:|---|---|
| `0` | 評価完了 | `velocities_for_spatial()` などで可視化・集約する |
| `1` | root box 内に正の確率がない | 速度範囲、空間点、`ispec`、`dt/max_step` を確認する |
| `2` | sample 容量に到達 | `max_samples_per_cell` を増やす、`max_depth` を下げる、速度 box を分ける |
| `3` | root box が拡張された | `leaf_bounds` や端の sample を確認し、必要なら初期速度範囲を広げる |
| `4` | leaf/box 容量に到達 | `max_leaves_per_cell` を増やす、`refine_threshold_rel` を大きくする |
| `5` | `vmin >= vmax` など不正な bounds | `velocity_bounds` の順序と shape を確認する |

`status != 0` でも partial sample は返ることがあります。解析に使う場合は、
まず status ごとの cell 数を集計し、問題のある cell を可視化してから本番設定を
決めるのが安全です。

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

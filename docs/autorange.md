# セルごとの速度範囲自動推定

> Lang: **日本語** | [English](autorange.en.md)

`estimate_velocity_range_map` は、EMSES の source 分布から決定論的な
support 粒子を生成し、Fortran 側で forward trace して、空間セルごとに
異なる `vx/vy/vz` 範囲を推定する API です。手動で全空間共通の速度 box を
決める代わりに、セルごとの `VelocityRangeMap` を作り、その範囲から
`get_probabilities` 用の粒子列を生成できます。

## 基本フロー

```python
import emout
from vdsolverf.emses import estimate_velocity_range_map, get_probabilities

data = emout.Emout("my-run")

range_map = estimate_velocity_range_map(
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=0.25,
    max_step=30_000,
    use_adaptive_dt=True,
    coverage_sigma=4.0,
    safety_factor=1.25,
    n_threads=4,
)

particles, index = range_map.create_particles(velocity_bins=(16, 8, 16))

probabilities, ret_particles = get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    use_adaptive_dt=False,
    n_threads=4,
)

probability_grid = index.reshape(probabilities)
```

`probability_grid` の形状は `(nz, ny, nx, nvz, nvy, nvx)` です。速度範囲が
推定できなかったセルは、`create_particles()` の既定ではスキップされます。

## validation 付きの使い方

sigma point 的な support 粒子は厳密な envelope ではないため、必要なら粗い
速度格子で `get_probabilities` を走らせ、速度 box の端に有意な確率が残る
セルだけ拡張します。

```python
from vdsolverf.emses import (
    estimate_velocity_range_map,
    validate_and_expand_velocity_range_map,
)

range_map = estimate_velocity_range_map(
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=0.25,
    max_step=30_000,
    use_adaptive_dt=True,
    coverage_mode="relative_density",
    eps_rel=1e-6,
    safety_factor=1.25,
    n_threads=4,
)

range_map = validate_and_expand_velocity_range_map(
    range_map,
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=data.inp.dt,
    max_step=30_000,
    coarse_bins=(8, 4, 8),
    edge_threshold=1e-3,
    expand_factor=1.5,
    max_iter=2,
    n_threads=4,
)
```

`status == 3` のセルは validation により拡張されたセルです。validation は
追加で `get_probabilities` を実行するので、本番解像度より粗い
`coarse_bins` を指定します。

## 主な引数

| 引数 | 目安 |
|---|---|
| `dt` | forward trace の刻み。`use_adaptive_dt=True` では 1 step の移動量を抑えるため、`0.25` から `0.5` 程度を初期値にする。 |
| `max_step` | source から対象領域を十分覆うまでの最大 step 数。反射後に戻る粒子も見たい場合は長めにする。 |
| `use_adaptive_dt` | セル飛びを避けるため、範囲推定では `True` を基本にする。 |
| `coverage_sigma` | source Maxwellian の support 半径。探索用は `3.0`、標準は `4.0`、保守的には `5.0`。 |
| `coverage_mode` / `eps_rel` | `coverage_mode="relative_density"` なら `eps_rel` から `sqrt(-2 log eps_rel)` を使う。`eps_rel=1e-6` は約 `5.26 sigma`。 |
| `safety_factor` | deposit 後の min/max を中心から拡張する係数。標準は `1.25`、取りこぼしが疑わしい場合は `1.5` 以上。 |
| `source_samples_per_cell` | source 面の接線方向 sampling 数。`1` が最速。source 面上の空間変化が強い場合は増やす。 |
| `velocity_bins` | `create_particles()` でセルごとに作る速度格子数。返り値の順序は `(nvx, nvy, nvz)`。 |
| `n_threads` | Fortran 側の OpenMP thread 数。未指定時は `OMP_NUM_THREADS`、なければ `1`。 |
| `collect_moments` | `True` で `mean_v` / `cov_v` を計算する。メモリ使用量が増えるため既定は `False`。 |

## 戻り値と diagnostics

`estimate_velocity_range_map` は `VelocityRangeMap` を返します。主要配列の
形状はすべて `(nz, ny, nx)` です。

| 属性 | 意味 |
|---|---|
| `vx_min`, `vx_max`, `vy_min`, `vy_max`, `vz_min`, `vz_max` | セルごとの速度範囲。hit がないセルは `nan`。 |
| `count` | support 粒子の hit 数。物理密度ではなく sampling 診断値。 |
| `weight_sum` | support 粒子重みの和。現行 MVP では診断用途。 |
| `mean_v`, `cov_v` | `collect_moments=True` のときだけ有効な速度 moment 診断。 |
| `status` | `0`: OK、`1`: LOW_COUNT、`2`: FALLBACK/no hit、`3`: EDGE_EXPANDED。 |
| `confidence` | `count / minimum_count` を 1 で切った簡易スコア。 |
| `metadata` | `coverage_sigma`、`safety_factor` など推定時の設定。 |

`VelocityRangeMap.create_particles()` は `(particles, index)` を返します。
`particles` は `get_probabilities` にそのまま渡せる flatten 済み粒子列です。
`index.reshape(probabilities)` を使うと、flatten された確率を
`(nz, ny, nx, nvz, nvy, nvx)` に戻せます。

## 反射と適用範囲

forward trace は実際の EM field と境界を使って粒子を進めるため、電位分布に
よる反射や、境界条件としての反射が軌道に現れれば、その後に通過したセルの
速度範囲にも反映されます。反射を拾えない典型例は、`max_step` が短い、
`dt` が大きくセルを飛ぶ、source support が粗く反射ローブを代表できていない
場合です。

現行実装が source envelope として扱うのは、EMSES の open boundary と
emission surface です。raycast 光電子や内部境界で生成される二次電子の
source envelope はまだ自動生成対象ではありません。そのような成分を使う場合は、
手動の速度範囲指定や validation 結果を併用してください。

## 性能メモ

範囲推定本体は Fortran 側で実行され、source patch 単位で OpenMP 並列化されます。
deposit は `atomic` ではなく thread-local accumulator に書き込み、最後に
merge します。このため atomic 競合は避けられますが、メモリ使用量は概ね
`セル数 * thread数` に比例します。`collect_moments=True` はさらに moment 用の
thread-local 配列を持つため、まずは `False` で使うのが安全です。

付属の簡易ベンチマークは次のように実行できます。

```bash
mpiifort -O3 -qopenmp -Iinclude benchmarks/bench_emses_autorange.f90 \
  -Llib -lvdist-solver-fortran -Wl,-rpath,"$PWD/lib" \
  -o /tmp/bench_emses_autorange

/tmp/bench_emses_autorange 4 32 32 32 96 5
```

引数は `n_threads lx ly lz max_step reps` です。KUDPC の login node では直接
実行せず、`tssrun` または batch job 内で実行してください。

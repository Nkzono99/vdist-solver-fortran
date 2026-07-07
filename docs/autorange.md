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

注意: `estimate_velocity_range_map` は既定で `use_adaptive_dt=False` です。
引数行をコメントアウトすると通常の時間刻み解釈になります。セル飛びを抑えたい
場合は `use_adaptive_dt=True` を明示してください。その場合の `dt` は EMSES の
`data.inp.dt` ではなく、1 step で進む格子距離の目安です。たとえば
`dt=0.004`, `max_step=30000` では、直線的に進んでも約 120 grid 分しか
source から届きません。

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
| `dt` | forward trace の刻み。`use_adaptive_dt=True` では 1 step の格子移動量の目安になり、EMSES の `data.inp.dt` ではない。`0.25` から `0.5` 程度を初期値にする。 |
| `max_step` | source から対象領域を十分覆うまでの最大 step 数。反射後に戻る粒子も見たい場合は長めにする。 |
| `use_adaptive_dt` | 既定は `False`。セル飛びを避ける場合は `True` を明示する。 |
| `coverage_sigma` | source Maxwellian の support 半径。探索用は `3.0`、標準は `4.0`、保守的には `5.0`。 |
| `coverage_mode` / `eps_rel` | `coverage_mode="relative_density"` なら `eps_rel` から `sqrt(-2 log eps_rel)` を使う。`eps_rel=1e-6` は約 `5.26 sigma`。 |
| `safety_factor` | deposit 後の min/max を中心から拡張する係数。標準は `1.25`、取りこぼしが疑わしい場合は `1.5` 以上。 |
| `source_samples_per_cell` | source 面の接線方向 sampling 数。`1` が最速。source 面上の空間変化が強い場合は増やす。 |
| `velocity_bins` | `create_particles()` でセルごとに作る速度格子数。返り値の順序は `(nvx, nvy, nvz)`。 |
| `n_threads` | Fortran 側の OpenMP thread 数。未指定時は `OMP_NUM_THREADS`、なければ `1`。 |
| `collect_moments` | `True` で `mean_v` / `cov_v` を計算する。メモリ使用量が増えるため既定は `False`。`False` の場合は巨大な moments 配列を確保せず、`range_map.mean_v` / `cov_v` は NaN view になる。 |
| `show_progress` | `True` で `get_probabilities` と同じ形式の progress bar を表示する。バッチログを抑えたい場合は `False`。 |
| `accumulator_cache_size` | 各 OpenMP thread が flush まで保持する hit cell 数。既定は `20000`。大きいほど flush 同期は減るがメモリは増える。 |

## メモリと並列化

`estimate_velocity_range_map` は、deposit ごとの atomic を避けるために各
thread に bounded sparse cache を持ちます。cache が
`accumulator_cache_size` に達すると、cell lock stripe で global range map
へまとめて reduction します。

メモリ使用量はおおよそ次の形です。

```text
global dense range arrays
+ EB field
+ n_threads * accumulator_cache_size * cache_entry_size
```

以前のような `全空間セル数 * n_threads` の dense thread-local 配列は確保し
ません。大規模格子で 112 thread を使う場合も、まずは既定の
`accumulator_cache_size=20000` から始め、flush が支配的なら増やします。

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

## 速度範囲と 6 次元分布へのアクセス

`VelocityRangeMap` は 6 次元分布そのものではなく、各空間セルの速度範囲を
保持します。保存される速度軸情報は、セルごとの min/max です。

```python
range_map.vx_min      # shape: (nz, ny, nx)
range_map.vx_max
range_map.vy_min
range_map.vy_max
range_map.vz_min
range_map.vz_max
range_map.valid_mask  # shape: (nz, ny, nx)
```

あるセルの速度軸は、可視化時に指定する `velocity_bins` から復元します。

```python
iz, iy, ix = 10, 20, 30
nvx, nvy, nvz = 16, 8, 8

vx = np.linspace(range_map.vx_min[iz, iy, ix], range_map.vx_max[iz, iy, ix], nvx)
vy = np.linspace(range_map.vy_min[iz, iy, ix], range_map.vy_max[iz, iy, ix], nvy)
vz = np.linspace(range_map.vz_min[iz, iy, ix], range_map.vz_max[iz, iy, ix], nvz)
```

`get_probabilities` の結果を 6 次元配列として見る場合は、
`create_particles()` が返す `index` を使います。

```python
velocity_bins = (16, 8, 8)  # (nvx, nvy, nvz)
particles, index = range_map.create_particles(velocity_bins=velocity_bins)

probabilities, _ = get_probabilities(..., particles=particles)
prob_grid = index.reshape(probabilities)

print(prob_grid.shape)
# (nz, ny, nx, nvz, nvy, nvx)

cell_probability = prob_grid[iz, iy, ix]
# shape: (nvz, nvy, nvx)
```

`range_map.save()` で保存されるのは `vx_min/vx_max` などの範囲と diagnostics
です。`prob_grid` や `vx/vy/vz` の全格子点は保存されません。したがって、
ロード後も同じ `velocity_bins` を指定すれば同じ速度軸を復元できます。

## セル単位アクセス

単一セルだけを詳しく見たい場合は、`range_map[iz, iy, ix]` で
`VelocityRangeCell` view を取得できます。

```python
cell = range_map[iz, iy, ix]

print(cell.position)     # セル中心 [x, y, z]
print(cell.vmin)         # [vx_min, vy_min, vz_min]
print(cell.vmax)         # [vx_max, vy_max, vz_max]
print(cell.count)
print(cell.status)

vx, vy, vz = cell.velocity_axes((16, 8, 8))
particles, index = cell.create_particles(velocity_bins=(16, 8, 8))
probabilities, _ = get_probabilities(..., particles=particles)

cell_probability = index.reshape(probabilities)
# shape: (nvz, nvy, nvx)
```

無効セルでは `cell.valid == False` になり、`cell.create_particles()` は既定で
空の粒子列を返します。

## 保存とロード

`estimate_velocity_range_map` で作った `VelocityRangeMap` は元の
`data.directory`、`ispec`、`istep` を保持します。そのため、引数なしの
`save()` で次の既定名に保存できます。

```python
path = range_map.save()
print(path)
# data.directory / "vdsolverf-velocity-range-map-ispec0-istep-1.npz"
```

ロードするときは、同じ `directory` / `ispec` / `istep` から既定パスを
組み立てられます。

```python
from vdsolverf.core import VelocityRangeMap

range_map = VelocityRangeMap.load(
    directory=data.directory,
    ispec=0,
    istep=-1,
)
```

任意のパスに保存する場合は `range_map.save("path/to/range-map.npz")`、
任意のパスから読む場合は `VelocityRangeMap.load("path/to/range-map.npz")`
を使います。

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

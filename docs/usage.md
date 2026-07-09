# 使用方法

> Lang: **日本語** | [English](usage.en.md)

本ドキュメントのサンプルはすべて、`emout` で読み込める EMSES 計算ディレクトリ
(例: `./my-run/`) が手元にあり、アクティブな Python 環境に
`vdist-solver-fortran` がインストール済みであることを前提にしています。

> [!Note]
> `ispec` は Python API での 0-オリジン粒子種インデックスです。
> 慣習: `0` = 電子、`1` = イオン、`2` = 光電子。
> 光電子の評価には `/emissn/` で `use_raycast = .true.` が必要です
> ([Namelist リファレンス](namelist.md#raycast-光電子)を参照)。

## 単一バックトレース

位相空間中の一点から粒子を時間逆行させ、確率タグ付き境界に到達したか、
そのときの確率を返します。

```python
import emout
import matplotlib.pyplot as plt
from vdsolverf.core import Particle
from vdsolverf.emses import get_backtrace

data = emout.Emout("my-run")

particle = Particle(position=[32, 32, 400], velocity=[0, 0, -10])

ts, probability, positions, velocities = get_backtrace(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particle=particle,
    dt=data.inp.dt,
    max_step=300_000,
    output_interval=1,
    use_adaptive_dt=False,
)

plt.plot(positions[:, 0], positions[:, 2])
plt.gcf().savefig("backtrace.png")
```

## 多粒子バックトレース

OpenMP で複数粒子を並列にトレースし、確率で重み付けして軌跡を重ね描きします。

```python
import emout
import matplotlib.pyplot as plt
import numpy as np
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_backtraces

data = emout.Emout("my-run")

NVX, NVZ = 50, 50
phase_grid = PhaseGrid(
    x=32, y=32, z=130,
    vx=(-100, 100, NVX),
    vy=0,
    vz=(-400, -360, NVZ),
)

particles = phase_grid.create_particles()

ts, probabilities, positions, velocities, last_indexes = get_backtraces(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=10_000,
    output_interval=1,
    use_adaptive_dt=False,
    n_threads=4,
)

maxp = np.nanmax(probabilities)
for p, pos, li in zip(probabilities, positions, last_indexes):
    if np.isnan(p):
        continue
    alpha = min(1.0, p / maxp)
    plt.scatter(pos[:li, 0], pos[:li, 2], s=0.1, color="black", alpha=alpha)

plt.gcf().savefig("backtraces.png")
```

`last_indexes[i]` は粒子 `i` の有効サンプル数で、それ以降の位置は 0 で
パディングされます。

## 位相空間確率ソルバ

軌跡を記録せず、各粒子の確率と最終ステップの位相空間点のみ返します。

```python
import emout
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_probabilities

data = emout.Emout("my-run")

phase_grid = PhaseGrid(
    x=32, y=32, z=130,
    vx=(-100, 100, 50),
    vy=0,
    vz=(-400, -360, 50),
)

particles = phase_grid.create_particles()

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
```

## Octree で速度空間を adaptive に探索する

速度空間の全点を一様格子で解く代わりに、各空間点で速度 box を octree に分割し、
確率が変化する領域を重点的に評価できます。重い計算は Fortran 側で実行され、
空間点方向に OpenMP 並列化されます。

```python
from vdsolverf.core import PhaseGrid
from vdsolverf.emses import get_probabilities_octree

phase_grid = PhaseGrid(
    x=(120, 180, 16),
    y=64,
    z=(300, 420, 16),
    vx=(-8.0e6, 8.0e6, 2),
    vy=(-4.0e6, 4.0e6, 2),
    vz=(-8.0e6, 8.0e6, 2),
)

octree = get_probabilities_octree(
    directory=data.directory,
    ispec=0,
    istep=-1,
    phase_grid=phase_grid,
    dt=data.inp.dt,
    max_step=30_000,
    scout_bins=(11, 9, 11),
    max_depth=5,
    n_threads=32,
)

if octree.status[0] != 0:
    print("warning: non-OK octree status", octree.status[0])

v = octree.velocities_for_spatial(0)
p = octree.probabilities_for_spatial(0)
```

戻り値は dense 6D 配列ではなく、compact な sample/leaf 配列です。複数の狭い
速度ローブを拾いたい場合は `scout_bins` と `max_samples_per_cell` を増やします。
アルゴリズム、KUDPC での `tssrun` 実行例、status 別の対処は
[Octree 速度空間確率ソルバ](octree_probabilities.md) を参照してください。

## セルごとの速度範囲を推定してから確率を解く

EMSES の open boundary と emission surface の source 分布から support
粒子を Fortran 側で forward trace し、空間セルごとの速度範囲を推定できます。
戻り値の `VelocityRangeMap` はセルごとに異なる `vx/vy/vz` 範囲を持ち、
`create_particles()` で `get_probabilities` に渡す粒子列へ flatten できます。

```python
from vdsolverf.emses import estimate_velocity_range_map, get_probabilities

range_map = estimate_velocity_range_map(
    directory=data.directory,
    ispec=0,
    istep=-1,
    dt=0.25,
    max_step=30_000,
    use_adaptive_dt=True,
    coverage_sigma=4.0,
    safety_factor=1.25,
    collect_moments=False,
    show_progress=True,  # 既定値。ログを抑える場合は False。
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
)

probability_grid = index.reshape(probabilities)
```

`range_map.count` は envelope support point の hit 数で、物理密度ではなく
診断値です。`count == 0` のセルは `create_particles()` の既定では
スキップされます。速度の `mean_v` / `cov_v` 診断も必要な場合は
`collect_moments=True` を指定します。詳しい引数の選び方と validation は
[セルごとの速度範囲自動推定](autorange.md) を参照してください。

## MPI 粒子並列

既存の `vdsolverf.emses.get_*` はそのまま OpenMP/スレッド並列の入口です。
MPI を使う場合は optional backend を明示します。

```python
from vdsolverf.emses.mpi import get_probabilities

probabilities, ret_particles = get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    n_threads=2,        # rank 内のスレッド数
)
```

この形は `srun -n 8 python script.py` のように Python スクリプト自体を
MPI 起動した場合に使います。`mpi4py` は optional 依存なので、
必要な環境だけ `pip install "vdist-solver-fortran[mpi]"` で追加します。

通常の Python プロセスから Slurm に投げる場合は launcher を使えます。

```python
from vdsolverf.emses.mpi import srun_get_probabilities

probabilities, ret_particles = srun_get_probabilities(
    directory=data.directory,
    ispec=0,
    istep=-1,
    particles=particles,
    dt=data.inp.dt,
    max_step=30_000,
    ntasks=8,
    n_threads=2,
    cpus_per_task=2,
)
```

## 共有ライブラリのパスを明示指定する

通常はラッパが OS を自動判定して同梱の共有ライブラリを読み込みますが、
明示オーバーライドも可能です。

```python
from vdsolverf.emses import get_backtrace

get_backtrace(..., system="linux", library_path="/custom/path/libvdist-solver-fortran.so")
```

`system` に指定できる値: `"auto"` (デフォルト)、`"linux"`、`"darwin"`、
`"windows"`。

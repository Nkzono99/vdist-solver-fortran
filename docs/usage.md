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

## ダスト粒子の帯電込みバックトレース

ダストソルバは EMSES の電流場を使い、ダスト粒子の軌跡と帯電を同時積分します。

```python
from vdsolverf.core import DustParticle
from vdsolverf.emses import get_dust_backtrace

dust = DustParticle(
    charge=-4e3, mass=1.0, radius=1.0,
    position=[32, 32, 400],
    velocity=[0, 0, -10],
)

ts, charges, positions, velocities = get_dust_backtrace(
    directory=data.directory,
    istep=-1,
    dust=dust,
    dt=data.inp.dt,
    max_step=10_000,
    use_adaptive_dt=False,
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

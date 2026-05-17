# vdist-solver-fortran

> Lang: **日本語** | [English](README.en.md)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14018863.svg)](https://doi.org/10.5281/zenodo.14018863)
[![CI](https://github.com/Nkzono99/vdist-solver-fortran/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Nkzono99/vdist-solver-fortran/actions/workflows/ci.yml)
[![PyPI version](https://img.shields.io/pypi/v/vdist-solver-fortran)](https://pypi.org/project/vdist-solver-fortran/)

Fortran で実装した速度分布ソルバを Python から利用するためのパッケージです。

コアは EMSES 計算結果の上でバックトレース (時間逆行) と、境界に紐付けた
ソース分布 (Maxwellian、raycast 光電子、吸収) を組み合わせて位相空間
確率密度を評価します。共有ライブラリを `ctypes` ラッパ経由で Python
から駆動する構成です。

## 必要環境

- `gfortran`
- `make`
- `fpm`
- Python 3.7 以上 (開発は `.venv/` 内の 3.12 で実施)

## インストール

PyPI からのインストールを推奨します。pip のビルド中に `make install` が
走り、Fortran 共有ライブラリをビルドして Python package に同梱します。

> [!Note]
> macOS でもビルド自体は通る想定ですが CI では検証していません。

```bash
python -m pip install -U pip setuptools wheel
python -m pip install vdist-solver-fortran
```

開発版を GitHub から直接入れることもできます。

```bash
python -m pip install "git+https://github.com/Nkzono99/vdist-solver-fortran.git"
```

pip 経由のビルドでは既定で `INSTALL_PROFILE=auto` を使います。必要なら
`INSTALL_PROFILE=generic` や `INSTALL_PROFILE=camphor` を環境変数で
指定してください。

```bash
INSTALL_PROFILE=generic python -m pip install vdist-solver-fortran
```

## クイックスタート

```python
import emout
from vdsolverf.core import Particle
from vdsolverf.emses import get_backtrace

data = emout.Emout("EMSES-simulation-directory")

ts, probability, positions, velocities = get_backtrace(
    directory=data.directory,
    ispec=0,                             # 0 電子, 1 イオン, 2 光電子
    istep=-1,
    particle=Particle([32, 32, 400], [0, 0, -10]),
    dt=data.inp.dt,
    max_step=300_000,
    output_interval=1,
    use_adaptive_dt=False,
)
```

光電子の評価には EMSES の namelist `/emissn/` に `use_raycast = .true.`
を設定する必要があります。詳細は
[Namelist リファレンス](docs/namelist.md#raycast-光電子) を参照。

## ドキュメント

| ドキュメント | 概要 |
|---|---|
| [使用方法](docs/usage.md) | Python レシピ集: 単一/多粒子バックトレース、位相空間確率 |
| [Namelist リファレンス](docs/namelist.md) | サポートしている `plasma.inp` のグループとパラメータ |
| [物理モデル](docs/physics.md) | リウビル定理、Maxwellian 放出、raycast 光電子 |
| [アーキテクチャ](docs/architecture.md) | Fortran / Python レイアウト、モジュール依存、拡張ポイント |
| [開発ガイド](docs/development.md) | ビルド、テスト、貢献手順 |

## サンプルノートブック

- [位相確率分布ソルバと多粒子バックトレース](https://nbviewer.org/github/Nkzono99/examples/blob/main/examples/vdist-solver-fortran/example.ipynb)

## ライセンス

Apache License 2.0。詳細は [LICENSE](LICENSE)。

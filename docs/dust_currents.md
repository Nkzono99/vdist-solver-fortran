# ダスト帯電電流モーメント

> Lang: **日本語** | [English](dust_currents.en.md)

`vdsolverf.emses.get_dust_currents` は、各空間点で速度 octree を一度だけ評価し、
`psi_grid` 全点の OML 積分を Fortran 内で計算します。サンプル配列を Python に
返す必要がないため、DRIFT の電流テーブル生成にも利用できます。

## 定義と単位

入力・出力は [`get_probabilities_octree`](octree_probabilities.md) と同じ
**EMSES 正規化単位**です。Python の `ispec` は 0 始まりで、C 境界で 1 始まりに
変換します。`qm` は選択した species の namelist 値を使います。

```text
F(psi) = integral p(v) |v| max(0, 1 - 2 qm psi / |v|^2) d^3v
j_s(psi) = q_s n_s F(psi)
I_s(psi) = pi a^2 j_s(psi)
```

`p(v)` は、境界源から backtrace で伝搬した密度 1 の PDF です。密度、電荷、
断面積、SI 単位への変換は呼び出し側で掛けます。`psi` は局所プラズマ電位からの
ダスト電位差で、`qm` と速度に整合する正規化を使います。速度は静止ダストの系で
与え、ダスト速度によるシフトは行いません。境界源につながらない捕捉粒子と、
ダストからプラズマへの逆作用は含みません。OML はダスト半径が局所 Debye 長より
十分小さく、捕集領域内の背景電位変化が小さいことを仮定します。表面近傍の捕集・
放出には別途物理モデルが必要です。

非ドリフト・等方 Maxwellian では、ソルバの熱速度は
`path = peth = vth_1d = sqrt(kT/m)`（1 成分の標準偏差）です。
`sqrt(2 kT/m)` ではありません。

```text
chi = qm psi / vth_1d^2
<v> = sqrt(8/pi) vth_1d
F = <v> exp(-chi)    (chi > 0)
F = <v> (1 - chi)    (chi <= 0)
```

## 使用例と出力

```python
from vdsolverf.emses import get_dust_currents
from vdsolverf.io import write_field_bundle

bundle = write_field_bundle("my-run", -1, [0, 1], "my-run/field_bundle")
result = get_dust_currents(
    "my-run", ispec=0, istep=-1,
    positions=[[12.5, 8.5, 4.5], [13.5, 8.5, 4.5]],
    velocity_bounds=[-5., 5., -5., 5., -5., 5.],
    psi_grid=[-1., 0., 1.], dt=0.1, max_step=30_000,
    n_threads=8, field_bundle=bundle,
)
j = species_charge * species_density * result.flux
print(result.status, result.sample_count, result.terminal_box_count)
```

`positions` は `(3,)` または `(nspatial, 3)`、`velocity_bounds` は `(6,)`、
`(3, 2)`、`(nspatial, 6)`、`(nspatial, 3, 2)` を受け付け、単一 box は全点へ
broadcast します。`psi_grid` は有限値からなる空でない 1 次元配列で、順序を保持
します。`DustCurrentResult` は `vdsolverf.core` からも import できます。

| 属性 | Python shape | 内容 |
|---|---|---|
| `positions` | `(nspatial, 3)` | 評価位置 |
| `psi_grid` | `(npsi,)` | 電位配列 |
| `flux`, `flux_error` | `(nspatial, npsi)` | 密度 1 のモーメントとサンプル変動推定値 |
| `status` | `(nspatial,)` | Octree status |
| `sample_count` | `(nspatial,)` | 中間 box を含む保存サンプル数 |
| `terminal_box_count` | `(nspatial,)` | 子を持たず積分した box 数 |
| `metadata` | `dict` | 入力設定、単位、積分則 |

場の有効・無効、platform/library 指定、一時 namelist の suffix は既存 octree と
同じです。既定値は `scout_bins=(9,9,9)`、`max_depth=6`、
`max_samples_per_cell=50000`、`max_leaves_per_cell=4096`、
`refine_threshold_rel=edge_threshold_rel=1e-4`、`expand_factor=1.5`、
`max_expansions=3`、`max_probability_types=100`、`use_adaptive_dt=False`。
`n_threads=None` は `OMP_NUM_THREADS`（未指定なら 1）に従います。

## 積分則と誤差推定

root の scout grid は box の両端を含む `linspace`、子 box は両端を含む
`3 x 3 x 3` 格子です（vx が最速）。各 **terminal box 自身のサンプル**に
テンソル積の台形重みを掛けます。各軸の端点重みは `1/2`、内部は `1` で、
格子間隔を掛け合わせます。分割済みの親は積分しません。積分領域は
**自動拡張後の最終 root** 全体です。無効な backtrace は寄与 0 とします。

引力側の被積分関数は原点に積分可能な `1/|v|` 特異性を持ちます。厳密に速度 0 の
格子点は、境界に届かない backtrace も含めて寄与を 0 とします。有限値を返すための
選択であり、特異セルの解析補正ではありません。強い引力では解像度の収束確認が必要です。

```text
flux_error = sum_terminal_boxes (max sampled integrand - min sampled integrand) * box_volume
```

通常の有限入力では非負・有限です。これはサンプル上の変動量で、**統計的信頼区間や
厳密な誤差上限ではありません**。未検出の狭い lobe、原点特異性、root 外の tail、
backtrace の時間刻み誤差は推定しません。

電流評価では分割前に子 8 個の box・サンプル容量を確保します。
容量不足なら親を terminal のまま保持し、領域欠落を避けます。有限予算では
各 terminal box 自身の格子と corner のみの台形積分の差（全 psi の最大値）を
誤差指標として、差が大きい box を優先して細分化します。確率に基づく分割条件と、
出力値に使う各 box 自身の格子の台形則は維持します。
status `2/4` のときも root 全体を覆う粗い積分を返しますが、収束確認が必要です。
root 格子自体がサンプル容量に収まらない場合は status `2`、flux/error と count は 0。
既存のサンプル返却 octree の途中終了動作は変更していません。

| Status | 意味 |
|---:|---|
| 0 | 指定した細分化を完了 |
| 1 | 正の境界源確率なし |
| 2 | サンプル容量が細分化または root 評価を制限 |
| 3 | Root 拡張あり。端の signal が残る場合もある |
| 4 | Box 容量が細分化を制限 |
| 5 | 不正な速度 bounds。出力 0 |

容量 status は拡張 status より優先し、両容量が制限した場合は box 容量を優先します。
最後の許容拡張後に root を再スコアしない点は既存 octree と同じです。

## Octree の terminal 情報

`get_probabilities_octree(...).leaf_is_terminal` は `(nleaf,)` の bool 配列です。
`result.leaf_mask(i) & result.leaf_is_terminal` で空間点 i の terminal box を選択
できます。既存 `leaf_*` の意味は変わらず、中間 box も含みます。既存 C symbol は
維持し、新しい `get_probabilities_octree_with_terminals` が探索中の flag を返します。
完了した木では terminal 体積和が最終 root 体積に一致します。既存 octree の
status `2/4` は途中終了なので、この分割を保証しません。
保存済み result や旧ライブラリでは、返された box の親子包含から flag を導出します。
途中終了で未返却の子は推定できません（旧ライブラリでは metadata の
`terminal_flags='returned_tree'`）。

両 API の `field_bundle=` は directory、`bundle.json` の path、`FieldBundle` を
受け付けます。元 directory、step、species、場の flag は bundle と一致する必要が
あり、不一致は `ValueError` です。bundle 呼び出しでは EMSES を開かず、namelist も
再生成しません。

## Field bundle 形式 0.1.0

`vdsolverf.io.write_field_bundle(directory, istep, ispecs, output_dir, *,
use_electric_field=True, use_magnetic_field=True)` は `FieldBundle` を返します。
`read_field_bundle(path)` は metadata のみ読み、`bundle.get_ebvalues(ispec)` または
`bundle[ispec]` が `np.load(..., mmap_mode='r')` で遅延読み込みします。

| ファイル | 内容 |
|---|---|
| `plasma-vdsolverf.inp` | geometry 変換・添字付きキーも含め、`TempolaryInput` と同一の namelist |
| `ebvalues_ispec{N}.npy` | Python の 0 始まり species N。shape `(nz+1, ny+1, nx+1, 9)` |
| `bundle.json` | `schema="vdsolverf.field_bundle"`、`version="0.1.0"`、絶対 path `directory`、`istep`、`nx/ny/nz`、`ispecs`、両 field flag、UTC ISO 8601 の `created_at`、`vdsolverf_version` |

配列は **npy v1.0、little-endian float64 (`descr='<f8'`)、C order**
(`fortran_order=False`)。先頭 6 byte が `0x93 NUMPY`、続いて version byte `1,0`、
2 byte little-endian unsigned header length、改行で終わる padding 付き ASCII の
Python dict（`descr`、`fortran_order`、`shape`）です。binary data は成分、x、y、z の
順に速く変化します。pickle・圧縮なし。Fortran 側では `(9,nx+1,ny+1,nz+1)` と
解釈します。header length を読み、固定 data offset を仮定しないでください。
standalone batch の Fortran reader は v1.0 と v2.0 に対応します。

## Standalone batch

`vdsolverf-batch [--resume|--no-resume] batch.nml` は各 rank で bundle と simulator を
一度だけ読み込みます。実行には Python / mpi4py は不要です。MPI は専用ビルドでのみ
有効になり、既定の `make` / `make install` / pip 用共有ライブラリは MPI 不要です。
依存追加なしの Fortran npy reader/writer は v1.0/v2.0、C order、little-endian の
`<f8` / `<i4` / `<i8` を扱い、別 dtype、big endian、Fortran order は拒否します。

```bash
make batch                       # gfortran、OpenMP 付き serial 版
make batch-mpi                   # mpif90 を優先、なければ mpiifx
# 別の wrapper: make batch-mpi MPI_FC=mpif90 MPI_FLAGS='-O3 -fopenmp'
build/batch-serial/bin/vdsolverf-batch batch/batch.nml
```

MPI 版は `build/batch-mpi/bin/vdsolverf-batch`。別 prefix で両実行ファイルを保持します。
`BATCH_FLAGS` / `MPI_FLAGS`、`BATCH_PREFIX` / `MPI_PREFIX` で上書きできます。
KUDPC ではビルドも実行も計算ノードで行います。

この環境では `mpif90` は GCC 8 + Intel MPI を使います。この組合せを検証構成と
します。Intel Fortran (`BATCH_FC=ifx` / `MPI_FC=mpiifx`) の既定フラグは
`-O3 -qopenmp -fp-model precise` です。以前の非ゼロ場での数値不一致は
ifx 2023.2 が生成する共有作業記述子の競合で、境界衝突時の allocatable な
確率ラッパのコピーも原因でした。確率関数を読み取り専用のまま直接呼び、
多相粒子は成分ごとにコピーし、solver の構築・解放を排他処理にしました。
回帰テストは19地点を8 threads で2回評価し、octree 全サンプルと最終軌道を
含む全出力のビットを単一 thread の結果と比較します。

```bash
fpm test --flag "-fopenmp" --target test_thread_determinism
fpm test --compiler ifx --flag "-O3 -qopenmp -fp-model precise" --target test_thread_determinism
```

ビット一致が必要ならコンパイラと浮動小数点フラグを揃えてください。
異なるコンパイラや fast-math の比較には数値許容誤差が必要です。
autorange の任意の平均・共分散は、ロックで保護された順序不定の加算を使います。
詳細は [autorange の再現性](autorange.md#再現性) を参照してください。

```python
import numpy as np
from vdsolverf.io import write_field_bundle, write_batch, read_batch_results

bundle = write_field_bundle("my-run", -1, [0], "field_bundle")
targets = np.array([
    [12.5, 8.5, 4.5, -5., 5., -5., 5., -5., 5.],
    [13.5, 8.5, 4.5, -5., 5., -5., 5., -5., 5.],
])
batch_nml = write_batch("batch", bundle, 0, targets, [-1., 0., 1.],
                        dt=0.1, max_step=30_000, chunk_size=16)
# 外部ジョブで vdsolverf-batch batch/batch.nml を実行後:
result = read_batch_results("batch")  # DustCurrentResult、元の target 順
```

`write_batch` は directory または `FieldBundle` を受け、`batch.nml` の Path を
返します。helper に渡す bundle path は呼び出し時の cwd から解決し、namelist
からの相対 path で保存します。`output_dir` は namelist からの相対 path または
絶対 path。`dt` の省略は bundle の `/tmgrid/ dt`、`n_threads` の省略または `None`
は実行時の OpenMP 設定を使います。

```text
batch/
  batch.nml
  targets.npy    # float64 C order (n,9): x,y,z,vxmin,vxmax,vymin,vymax,vzmin,vzmax
  psi.npy        # float64 C order (npsi,): EMSES 電位、入力順を保持
  out/
```

両配列は空にできません。位置と psi は有限値を要求し、不正な速度 box は既存 evaluator
と同じ status `5` を返します。`ispec` は Python と同じ 0 始まりです。全入出力は前述の
EMSES 単位・密度 1 の flux を使います。

`batch.nml` の group は `&vdsolverf_batch` です。

| キー | 既定値 / 内容 |
|---|---|
| `bundle_dir`, `ispec`, `dt` | 手書き入力では必須。bundle は namelist からの相対 path 可 |
| `max_step`, `use_adaptive_dt`, `max_probability_types` | `30000`, `.false.`, `100` |
| `scout_bins`, `max_depth` | `9,9,9`, `6` |
| `max_samples_per_cell`, `max_leaves_per_cell` | `50000`, `4096` |
| `refine_threshold_rel`, `edge_threshold_rel` | 両方 `1e-4` |
| `expand_factor`, `max_expansions` | `1.5`, `3` |
| `chunk_size`, `output_dir` | `16`, `'out'` |
| `n_threads` | OpenMP runtime の既定（通常 `OMP_NUM_THREADS`、OpenMP 無効なら 1） |
| `use_electric_field`, `use_magnetic_field` | 両方 `.true.`。bundle の flag と完全一致が必要 |

field flag の不一致は停止します。必要な flag で別 bundle を生成してください。

### Chunk と再開

chunk の開始 target index は **0 始まり**。最後の chunk だけ短くできます。

| ファイル (`start=0` の例) | dtype / C-order shape |
|---|---|
| `chunk_000000000_flux.npy` | `<f8` `(count,npsi)` |
| `chunk_000000000_flux_error.npy` | `<f8` `(count,npsi)` |
| `chunk_000000000_status.npy` | `<i4` `(count,)` |
| `chunk_000000000_sample_count.npy` | `<i8` `(count,)` |
| `chunk_000000000_terminal_box_count.npy` | `<i8` `(count,)` |
| `chunk_000000000.done` | 空の完了 marker。5 配列を閉じた後、最後に atomic rename |

既定の resume（`--resume` も可）は marker と 5 配列が揃う chunk を skip します。
marker がなければ途中の `.npy` / `.tmp` が残っていても再計算します。
`--no-resume` は各 marker を外して全 chunk を再計算します。bundle、targets、psi、
solver 設定、chunk_size は再開中に変更しないでください。変更時は新しい output
directory または `--no-resume` を使います。入力 hash による鮮度判定は行いません。
同じ output directory への同時実行はサポートしません。

全 rank 完了後、rank 0 が `summary.json` を atomic rename で保存します。
`status_counts`（全出力の status 別件数）、`n_targets`、`npsi`、`wall_time_seconds`
（読み込み・計算・待ちを含む経過時間）、`ranks`、`threads`、`chunks_computed` /
`chunks_skipped`（今回の実行分）を記録します。`read_batch_results` は全 marker、shape、
dtype を検査し、部分結果は返しません。summary は `result.metadata['summary']` に入り、
全 chunk 完了直後の停止で summary だけ未生成でも集約できます。batch の count は
int64 を保持し、ctypes の既存 int32 count は変えません。

### MPI と KUDPC

MPI rank が 2 以上なら rank 0 は配布専用、他 rank は worker です。worker は完了後に
次の chunk を要求し、rank 0 が未処理 index または終了指示を返します。rank 0 は
octree を実行しないため重い target で配布が停止しません。1 rank の MPI 版と serial
版は rank 0 が全 chunk を計算します。各 worker 内では既存 evaluator の OpenMP が
target を動的に処理します。

```bash
# 計算ノードでビルド。controller は一度だけ起動
tssrun -p <partition> -t 0:40:00 --rsc p=1:t=4:c=4 make batch-mpi

# PATH に build/batch-mpi/bin を追加してから投入
tssrun -p <partition> -t 2:00:00 --rsc p=16:t=8:c=8 srun vdsolverf-batch batch/batch.nml
```

2026-10-15 の KUDPC smoke test は次の形で成功しました（リポジトリ直下から）。
GCC 8.5 + Intel MPI の 4 ranks / 各 4 threads、19 targets / 10 chunks で serial と
全配列が bit 一致し、同じ入力の `--resume` は全 10 chunks を skip しました。
PMI 環境の変更や hydra への切り替えは不要でした。

```bash
tssrun -p gr20001b -t 0:40:00 --rsc p=4:t=4:c=4 \
  srun build/batch-mpi/bin/vdsolverf-batch --no-resume build/batch-smoke/mpi/batch.nml
```

再現用入力は `.venv/bin/python test/batch_mpi_smoke.py prepare build/batch-smoke`、
serial 計算後の比較は同 script の `compare build/batch-smoke` で行えます。

`n_threads` は確保した thread/core 数以下にします。場は rank ごとに保持し、dispatcher
も 1 rank 分のメモリを使います。chunk を小さくすると負荷を均しやすい一方、通信と
ファイル数は増えます。既定の 16 は出発点です。chunk 内の target 数より多い thread
は有効に使えません。scratch は worker の thread 数と octree 予算に比例し、全 target
数には比例しません。
既存 octree は評価中の target ごとに solver/simulator を複製するため、場のコピーも
同時評価 target 数に応じた peak memory に含めて見積もってください。

## 性能と検証

OpenMP は空間点を動的分担し、速度探索を全 `psi` で共用します。scratch は
`n_threads * max_samples_per_cell` と `n_threads * max_leaves_per_cell` に比例し、
全空間点のサンプル総数には比例しません。出力は `nspatial * npsi` に比例します。
拡張判定用 scout は既存 octree と同じく `sample_count` に含みません。

bundle で場の読み込み・relocation を再利用できます。物理・熱速度の推定や autorange
で bounds を選び、容量、scout 間隔、tail、時間刻み、深さの収束を個別に確認してください。
既定設定の Maxwellian 回帰が 1% でも、任意 VDF や bounds に同精度を保証しません。
KUDPC の backtrace は計算ノードで実行し、thread 数を資源に合わせます。

```bash
tssrun -p gr20001b -t 0:40:00 --rsc p=1:t=8:c=8 .venv/bin/python make_currents.py
```

## C ABI

`evaluate_dust_currents` は返値なしの `bind(c)` で、引数順は次のとおりです。

```text
inppath, length, lx, ly, lz, ebvalues, ispec, nspatial,
spatial_points, velocity_bounds, npsi, psi_grid, dt, max_step,
use_adaptive_dt, max_probability_types, scout_nvx, scout_nvy, scout_nvz,
max_depth, max_samples_per_cell, max_leaves_per_cell,
refine_threshold_rel, edge_threshold_rel, expand_factor, max_expansions,
flux, flux_error, status, sample_count, terminal_box_count, n_threads
```

整数 scalar は `integer(c_int), value`、実数 scalar は `real(c_double), value`。
例外は optional の **参照渡し** `n_threads`（ctypes では `byref(c_int(...))`）。
path は `character(c_char)` byte 列と byte 長。実数配列は `real(c_double)`、整数配列は
`integer(c_int)` です。Fortran shape は位置 `(3,nspatial)`、bounds `(6,nspatial)`、
psi `(npsi)`、moment `(npsi,nspatial)`、diagnostics `(nspatial)`。
Python は逆順 shape の C-contiguous float64/int32 配列を使います。
新 octree symbol は既存引数列の末尾に容量 `nspatial*max_leaves_per_cell` の
`integer(c_int)` flag 出力（0/1）を追加します。既存 symbol の ABI は変更しません。

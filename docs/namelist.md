# Namelist リファレンス

> Lang: **日本語** | [English](namelist.en.md)

`vdist-solver-fortran` は EMSES の `plasma.inp` のうち一部のグループだけを
読み込みます。Python ラッパが必要なグループだけを一時ファイル
(`plasma-vdsolverf.inp`) に書き出してから Fortran の共有ライブラリを
呼び出す仕組みで、
[`tmpolary_input.TMP_INP_KEYS`](../vdsolverf/emses/tmpolary_input.py) に
列挙されたキーだけが渡されます。

未対応のキーは黙って無視されます。新しいキーを追加したいときは Python 側の
`TMP_INP_KEYS` と Fortran 側の `m_namelist` を両方編集してください。

## グループ一覧

| グループ | 用途 |
|---|---|
| `/esorem/` | `emflag` — 静電 / 電磁モードの切替 |
| `/plasma/` | プラズマパラメータ (`wp`、`wc`、背景磁場角 `phixy`、`phiz`) |
| `/tmgrid/` | 時間刻み `dt`、グリッドサイズ `nx`、`ny`、`nz` |
| `/system/` | 粒子種数 `nspec`、境界コード `npbnd(3, nspec)` |
| `/intp/` | 粒子種ごとの運動パラメータ (`qm`、`path`、`peth`、`vdri`、`vdthz`、`vdthxy`、`spa`、`spe`、`speth`) |
| `/ptcond/` | 内部境界の形状・反射特性 |
| `/emissn/` | 放出面設定と raycast 光電子設定 |

本ソルバの確率計算に最も影響するのは `/ptcond/` (どこに境界があるか) と
`/emissn/` (衝突時にどう確率を割り当てるか) の 2 グループです。

## /ptcond/ — 内部境界

エントリは 2 系統あります。`boundary_type` / `boundary_types(:)` (型付き
境界ライブラリ) と、`geotype` (単純な基本形状) です。

### boundary_type

単一境界タイプ。取り得る値:

```
'none'
'flat-surface'
'rectangle-hole' | 'cylinder-hole' | 'hyperboloid-hole' | 'ellipsoid-hole'
'rectangle[xyz]' | 'circle[x/y/z]' | 'cuboid' | 'disk[x/y/z]'
'complex'          ! boundary_types(:) で複数を指定したいとき
```

### 種別ごとのパラメータ

| 種別 | キー |
|---|---|
| `flat-surface` / `*-hole` | `zssurf` (表面高さ、グリッド単位)。穴タイプは `[x/y/z][l/u]pc` |
| `complex` | `boundary_types(ntypes)` — タイプ名の配列 (上記参照) |
| `rectangle` | `rectangle_shape(ntypes, 6)` = `(xmin, xmax, ymin, ymax, zmin, zmax)` |
| `circle[x/y/z]` | `circle_origin(ntypes, 3)`、`circle_radius(ntypes)` |
| `cuboid` | `cuboid_shape(ntypes, 6)` = `(xmin, xmax, ymin, ymax, zmin, zmax)` |
| `disk[x/y/z]` | `disk_origin(ntypes, 3)`、`disk_height(ntypes)`、`disk_radius(ntypes)`、`disk_inner_radius(ntypes)` |
| 全体回転 | `boundary_rotation_deg(3)` (degrees) |

### geotype (簡易形状)

| キー | 意味 |
|---|---|
| `npc` | geotype オブジェクト数 |
| `geotype(npc)` | `0`–`1` 直方体、`2` 円柱、`3` 球 |
| 直方体 | `xlpc`, `xupc`, `ylpc`, `yupc`, `zlpc`, `zupc` |
| 円柱 | `bdyalign` (1=X, 2=Y, 3=Z)、`bdyedge(1:2)` (軸方向範囲)、`bdyradius`、`bdycoord(1:2)` (軸上中心座標) |
| 球 | `bdyradius`、`bdycoord(1:3)` |

## /emissn/ — 粒子放出

```
nflag_emit(nspec)  = 0 吸収、1 表面放出、2 光電子
nepl(nspec)        = この種に対する明示放出面の数
nemd(nepl)         = 放出面法線 (符号 = 正負方向、絶対値 = 軸 1=X, 2=Y, 3=Z)
curf(nspec)        = 粒子種ごとの放出電流密度
curfs(nepl)        = 放出面ごとの電流密度 (指定時は curf より優先)
xmine, xmaxe, ...  = 放出面のバウンディングボックス (nepl ごと)
thetaz, thetaxy    = 放出粒子の熱速度軸の傾き (nepl ごと [deg])
```

### Raycast 光電子

`use_raycast = .true.` **かつ** `nflag_emit(ispec) == 2` のとき、内部境界
への衝突は
[`m_photoelectron_raycast`](../src/emses/photoelectron_raycast.f90) で
定義された raycast 光電子確率で評価されます。モデルの数理は
[物理モデル](physics.md#raycast-光電子) を参照。フラグ未設定時は既定の
ZeroProbability (吸収) にフォールバックします。

| キー | 既定値 | 意味 |
|---|---|---|
| `use_raycast` | `.false.` | マスタースイッチ。光電子種に対する raycast 確率を有効化 |
| `ray_zenith_angle_deg(nspec)` | `9999d0` | 太陽方向の天頂角オーバーライド。センチネル `9999d0` は `vdthz(ispec)` にフォールバック |
| `ray_azimuth_angle_deg(nspec)` | `9999d0` | 太陽方向の方位角オーバーライド。センチネル時は `vdthxy(ispec)` にフォールバック |

これらの角度から計算したドリフトベクトルの「反対方向」に衝突点から
レイを飛ばし、内部境界に遮蔽されれば確率 0、そうでなければ
`locs = vdri_vector(ispec)` と `scales = vth_vector(ispec)` で定義される
シフト Maxwell PDF (半空間正規化のため 2 倍) を返します。

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
| `/intp/` | 粒子種ごとの運動パラメータ (`qm`、`path`、`peth`、`vdri`、`vdthz`、`vdthxy`、`spa`、`spe`、`speth`、`vdx`、`vdy`、`vdz`) |

角度 (`phiz`、`phixy`、`vdthz`、`vdthxy`、`speth`) はすべて度単位で、EMSES の規約に
従います。背景磁場は `(sin(phiz)cos(phixy), sin(phiz)sin(phixy), cos(phiz))` 方向、
`spa`/`path` はそれに平行、`spe` は垂直面内で角度 `speth` の方向、`vdri` は +z からの
天頂角 `vdthz`、+x からの方位角 `vdthxy`、`vdx`/`vdy`/`vdz` は Cartesian のドリフトです
([Physics](physics.md) 参照)。
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
| `rectangle-corners` | `rectangle_corners(3, 4, ntypes)` = 任意向きの長方形の 4 頂点 |
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
nflag_emit(nspec)  = 0 外側境界から注入、1 表面放出、2 光電子 (3 = SEE は非対応)
nepl(nspec)        = この種に対する明示放出面の数
nemd(nepl)         = 放出面法線 (符号 = 正負方向、絶対値 = 軸 1=X, 2=Y, 3=Z)
curf(nspec)        = 粒子種ごとの放出電流密度 (パッチでは正のときだけ使う)
curfs(nepl)        = 放出面ごとの電流密度 (正のとき curf より優先)
abvdem(nspec)      = 表面放出の法線方向に加える平均速度
xmine, xmaxe, ...  = 放出面のバウンディングボックス (nepl ごと)
thetaz, thetaxy    = 光電子パッチの入射方向 (スカラー [deg])
decrease_current_according_to_incident_angle
                   = 光電子パッチの電流を入射角の余弦で減らす (既定 .true.)
```

パッチの平均速度・熱速度、電流の選び方と入射角の扱いは
[物理モデル](physics.md#放出パッチの規格化) を参照。EMSES の挙動を再現できない
設定 (`curf`・`curfs` がどちらも正でないパッチ、`vdri ≠ 0` のパッチ放出種、
SEE・反射・loss-cone・周期 raycast・`omniemit`・`mainprop < 1` の二成分放出) は
`ValueError` で拒否します。

### Raycast 光電子

`use_raycast = .true.` (EMSES 入力の `pe_ray_cast` は Python が変換します)
**かつ** `nflag_emit(ispec) == 2` のとき、内部境界
への衝突は
[`m_photoelectron_raycast`](../src/emses/photoelectron_raycast.f90) で
定義された raycast 光電子確率で評価されます。モデルの数理は
[物理モデル](physics.md#raycast-光電子) を参照。フラグ未設定時は既定の
ZeroProbability (吸収) にフォールバックします。

| キー | 既定値 | 意味 |
|---|---|---|
| `use_raycast` | `.false.` | マスタースイッチ。光電子種に対する raycast 確率を有効化 |
| `ray_zenith_angle_deg(nspec)` | `9999d0` | EMSES カメラの天頂角。センチネル `9999d0` は `vdthz(1)` (粒子種 1) にフォールバック |
| `ray_azimuth_angle_deg(nspec)` | `9999d0` | EMSES カメラの方位角。センチネル時は `vdthxy(1) - 180` にフォールバック |

太陽方向は EMSES と同じ finbound のカメラから求めます。衝突点が日陰なら確率 0、
そうでなければ交差面の法線を使った表面源の密度を返します
([物理モデル](physics.md#raycast-光電子))。

# 物理モデル

> Lang: **日本語** | [English](physics.en.md)

本ソルバは、時間逆行でのバックトレース軌道計算と、境界にタグ付けされた
ソース分布の評価を組み合わせて、位相空間の確率密度を計算します。中心的な
道具はリウビルの定理と、[`src/core/probabilities.f90`](../src/core/probabilities.f90) および
[`src/emses/photoelectron_raycast.f90`](../src/emses/photoelectron_raycast.f90)
で定義される `t_Probability` 実装群です。

## リウビルの定理

位置と時間についてなめらかな力場中の無衝突系では、1 粒子分布関数
$f(\mathbf{x}, \mathbf{v}, t)$ は粒子軌道に沿って一定です:

$$
\frac{df}{dt} = \partial_t f + \mathbf{v} \cdot \nabla_{\mathbf{x}} f + \mathbf{a} \cdot \nabla_{\mathbf{v}} f = 0.
$$

観測点 $(\mathbf{x}_0, \mathbf{v}_0, t_0)$ から軌道を時間逆行させ、
ソース分布 $f_\Sigma$ が既知の境界 $\Sigma$ 上の点 $(\mathbf{x}_\Sigma,
\mathbf{v}_\Sigma)$ に到達したならば、

$$
f(\mathbf{x}_0, \mathbf{v}_0, t_0) = f_\Sigma(\mathbf{x}_\Sigma, \mathbf{v}_\Sigma).
$$

[`src/core/solver.f90`](../src/core/solver.f90) の
`solver_calculate_probability` はまさにこの操作を実装しています: 運動
方程式を逆向きに積分し、境界に衝突した瞬間、その境界に紐付いた確率関数
を参照します。

## Boris 更新と逆追跡

磁場回転には `s = 2*t/(1 + sum(t*t))` を使います。分母は磁場に対応する
ベクトル `t` のノルム二乗を含む共通のスカラーで、電場がない場合は
磁場の向きによらず速さと磁場方向の速度成分を保存します。

MPIEMSES3D の通常の Boris 更新は、現在位置の場で速度を更新してから、
更新後の速度で位置を進めます。逆追跡ではこの順序を逆にし、まず現在の速度で
位置を戻し、その位置の場を使って速度更新を逆に行います。周期境界を
またいだ場合は、場を評価する位置を周期領域内に戻します。

同じ静的な場・同じ時間刻みを使い、境界衝突のないステップを往復した場合、
位置と速度は丸め誤差の範囲で元に戻ります。可変刻みや衝突処理を含む
軌道全体について、この往復性を保証するものではありません。

境界衝突の位置・速度・時刻は同じ線分上の割合 `r` で補間します。
時刻は `t_before + r * dt_step` で、`dt_step` の符号と可変刻みを反映します。

## 対応するソース分布

### Zero (吸収)

`t_ZeroProbability%at` は常に 0。逆トレースを「吸い込み」にしたい壁に使います。

### シフト Maxwell 分布

`t_MaxwellianProbability%at` は 3D Maxwell–Boltzmann 密度:

$$
f_M(\mathbf{v}) = C \prod_{i=1}^{3} \frac{1}{\sqrt{2\pi}\sigma_i}
  \exp\!\left(-\frac{u_i^2}{2 \sigma_i^2}\right),
$$

$\boldsymbol\mu$ = `locs` はドリフト速度、$\boldsymbol\sigma$ = `scales`
は主軸方向の熱速度、$C$ = `coefficient` は乗算係数です。
$\mathbf u = R^T(\mathbf v-\boldsymbol\mu)$ として、速度を分布の基準座標へ
戻して評価します。省略時の `basis` ($R$) は単位行列です。

以下の EMSES 源は MPIEMSES3D が実際に粒子を生成する式に合わせています。
対応する EMSES のソースと追従手順は [EMSES との同期](emses_sync.md) にまとめています。

**外側境界からの注入**（`npbnd = 2`、`nflag_emit = 0`）は
`locs = vdri_vector(ispec)`、`scales = [peth, peth, path]`、`basis = thermal_basis()` です。
EMSES の初期配置と外側境界の注入（4.23.0 以降、`m_thermal_frame.F90`）と同じく、
熱分布は磁場座標で定義され、共分散は $R\,\mathrm{diag}(peth^2, peth^2, path^2)R^T$ です。
斜めの異方分布では Cartesian 成分の相関を保ち、等方分布は回転で変化しません。
ドリフトは EMSES の `vdtx/vdty/vdtz` と同じく

$$
\boldsymbol\mu = (vdx, vdy, vdz) + R\,(spe\cos(speth),\ spe\sin(speth),\ spa)
+ vdri\,(\sin\vartheta\cos\varphi,\ \sin\vartheta\sin\varphi,\ \cos\vartheta)
$$

です（$\vartheta$ = `vdthz`、$\varphi$ = `vdthxy`、角度はすべて度単位）。
$R$ = `thermal_basis()` は EMSES の背景磁場フレーム（`t11`..`t33`）を列とし、
第 3 列が Python ラッパーの詰める背景磁場
$\hat{\mathbf B}_0 = (\sin\phi_z\cos\phi_{xy},\ \sin\phi_z\sin\phi_{xy},\ \cos\phi_z)$
に一致します。`vth_vector(ispec)` はこの共分散の対角成分の平方根（Cartesian 各成分の
標準偏差）で、EMSES の境界テーブルの幅と同じです。autorange の support にも使います。

### 放出パッチの規格化

表面源（放出パッチと raycast 光電子）は、局所座標（外向き法線、接線 1、接線 2）で
平均 `[spa + abvdem, spe, spe]`、標準偏差 `[path, peth, peth]` を使います。
EMSES は法線速度を平均 `spa + abvdem` の flux-weighted 片側 Gaussian の表から、
両接線を同じ Gaussian（平均 `spe`）の表から独立に取り、`speth` は使いません
（`surface_emission_source.F90`、`direct_surface_emitter.F90`）。接線の符号は
パッチの向き（`nemd` の符号）によらず正のままです。

$u = spa + abvdem$、$\sigma = path$ とすると、密度 1 の Gaussian の
一方向フラックスモーメントは

$$
I(u,\sigma)=\frac{u}{2}\,\mathrm{erfc}\left(-\frac{u}{\sqrt{2}\sigma}\right)
+\frac{\sigma}{\sqrt{2\pi}}\exp\left(-\frac{u^2}{2\sigma^2}\right)
$$

です。$1+\mathrm{erf}$ の代わりに $\mathrm{erfc}$ を使い、ドリフトが表面に
向かう場合（$u \ll -\sigma$）の桁落ちを避けています。係数は `C = J/I` とし、
外向き（$\mathbf v\cdot\hat{\mathbf n} > 0$）の半空間だけで $C f_M$ を返します。
この分布の外向き一次モーメントが $J$ を再現します。

$J$ は EMSES と同じ順で決めます: `curfs(iepl) > 0` ならその値、そうでなく
`curf(ispec) > 0` ならその値。どちらも正でない場合、EMSES は `flpf`・`flpfs`
や `npin` 由来の密度を使いますが、本ソルバはそれらを読まないので入力を拒否します。
光電子パッチ（`nflag_emit = 2` で `pe_ray_cast = .false.`）では、EMSES と同様に
`decrease_current_according_to_incident_angle`（既定 `.true.`）のとき
$J$ に $\mathrm{sign}(nemd)\,\hat d_{|nemd|}$ を掛け、負なら 0 にします。
$\hat{\mathbf d}$ は `thetaz`・`thetaxy`（スカラー、度単位）から作る入射方向の
単位ベクトルです。EMSES の `vdri` によるパッチ放出の平行移動は半空間をずらすため
再現できず、`vdri ≠ 0` のパッチ放出種は拒否します。

`curf/curfs` は EMSES 単位の電流密度です。MPIEMSES3D はこれをマクロ粒子電荷の
絶対値で割って生成粒子数を求めます。本ソルバの $C$ は電荷換算密度（実粒子の
電荷絶対値が 1 なら数密度）であり、連続分布にマクロ粒子電荷・重みを掛けては
いけません。電荷絶対値が異なる場合の数密度はその値で割って求めます。
この源は密度 1 の PDF ではありません。
[ダスト電流の単位](dust_currents.md#定義と単位)も参照してください。

### Raycast 光電子

`use_raycast = .true.`（EMSES 入力の `pe_ray_cast`）と `nflag_emit(ispec) == 2`
の組み合わせで有効化されます。EMSES は太陽方向から平行カメラでレイを撃ち、
最初に当たった放出面（`material%tag /= -1`）の点から、交差面の法線
$\hat{\mathbf n}$ と `raycast_sampler.F90` の接線基底で上の表面源の速度を
1 粒子ずつ生成します。射影面積あたりのレイ数が一定なので、面積あたりの放出流束は
$|curf|\max(\hat{\mathbf n}\cdot\hat{\mathbf s}, 0)$ です（$\hat{\mathbf s}$ は
表面から太陽への単位ベクトル）。

本ソルバは衝突点から太陽側へ少しずらした点からカメラと同じ向きにレイを戻して
交差面と法線を求め、

$$
f_{\rm PE}(\mathbf v) = \frac{|curf|\,\max(\hat{\mathbf n}\cdot\hat{\mathbf s},0)}{I}
\prod_{i=1}^{3}\frac{1}{\sqrt{2\pi}\sigma_i}\exp\!\left(-\frac{(w_i-\mu_i)^2}{2\sigma_i^2}\right),
\qquad w_1 = \mathbf v\cdot\hat{\mathbf n} > 0
$$

を返します。$\mathbf w$ は（法線、接線 1、接線 2）成分、$\boldsymbol\mu$ と
$\boldsymbol\sigma$ は上の表面源の値です。次の場合は 0 です。

- 衝突点が放出面上にない、または EMSES の非放出カバー（tag -1）上にある。
- 速度が表面へ向かう（$w_1 \le 0$）。
- 太陽方向へのレイが他の境界に当たる（日陰）。

遮蔽判定の境界リストは EMSES と同じく `cover_all` 付きで作ります。
太陽方向は EMSES と同じ `new_ParallelCamera_optimized(thz, thxy, domain)` の
レイ方向の逆です。`ray_zenith_angle_deg(ispec)` が未設定（9999）なら `vdthz(1)`、
`ray_azimuth_angle_deg(ispec)` が未設定なら `vdthxy(1) - 180` を使います
（どちらも粒子種 1 の値）。カメラの回転は `direction_from_angles` と符号規約が
異なるため、カメラそのものから求めます。周期 raycast（`periodic_raycast`）と
時間変化する放出（`emission_time_mode`）は再現しません。

## 座標と単位の規約

- 位置は EMSES グリッド単位 ($[0, nx] \times [0, ny] \times [0, nz]$)。
- 速度は EMSES の自然単位系。スケールは `/intp/` の `vdri`, `path`,
  `peth` と整合。
- 境界コード `npbnd`: `0` 周期、`1` 反射、`2` 開放 (Maxwellian 源候補)、
  `3` 吸収。
- 粒子種インデックス `ispec` は **Fortran 側で 1-オリジン**、**Python
  API で 0-オリジン**。ラッパが自動変換。

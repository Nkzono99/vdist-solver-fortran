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

EMSES の外部源と raycast 光電子では、`locs = vdri_vector(ispec)`、
`scales = [peth(ispec), peth(ispec), path(ispec)]` を使います。
`thermal_basis()` はドリフトと同じ `phiz`、`phixy` 回転を主軸へ適用します。
共分散は $R\,\mathrm{diag}(\sigma_i^2)R^T$ で、斜めの異方分布の相関を保持し、
等方分布は回転で変化しません。`vth_vector(ispec)` はこの共分散の対角成分の
平方根（Cartesian 各成分の標準偏差）を返し、autorange の support に使います。
放出パッチは `nemd` の符号を含む法線方向に `path`、両接線方向に `peth`
を標準偏差として使います（標準偏差自体に方向の符号は付きません）。

### 放出パッチの規格化

放出方向（物質からプラズマへ）のドリフトを
$u=\mathrm{sign}(nemd)\,\mu_{|nemd|}$、$\sigma=path$ とすると、
密度 1 の Gaussian の一方向フラックスモーメントは

$$
I(u,\sigma)=\frac{u}{2}\left[1+\mathrm{erf}\left(\frac{u}{\sqrt{2}\sigma}\right)\right]
+\frac{\sigma}{\sqrt{2\pi}}\exp\left(-\frac{u^2}{2\sigma^2}\right)
$$

です。係数は `C = curf/I`（非負の `curfs` があれば優先）とし、分布の
一方向モーメントが入力値を再現します。MPIEMSES3D は法線速度を
$v_n f(v_n)$ に比例する flux-weighted half-Gaussian（標準偏差 `path`）、
接線速度を標準偏差 `peth` の Gaussian から生成します。本ソルバで伝搬する
位相空間分布は $C f_M$ であり、速度をもう一度掛けることはありません。
`curf/curfs` は EMSES 単位の電流密度の絶対値です。MPIEMSES3D はこれを
マクロ粒子電荷の絶対値で割って生成粒子数を求めます。本ソルバの $C$ は
電荷換算密度（実粒子の電荷絶対値が 1 なら数密度）であり、連続分布に
マクロ粒子電荷・重みを掛けてはいけません。電荷絶対値が異なる場合の数密度は
その値で割って求めます。この源は密度 1 の PDF ではありません。

erf ユーティリティの幅は $\sqrt{2}\sigma$ です。引数逆転と幅を修正すると、
ゼロドリフトでは $C=curf\sqrt{2\pi}/\sigma$ となります。従来も正だった法線幅が
同じ場合（例: ±Z パッチ）、放出確率は $\sqrt{2\pi}\simeq2.5066$ 倍です。
X/Y パッチでは従来の符号付き・入れ替わった熱幅も修正するため、以前ゼロだった
源が非ゼロになるケースがあり、共通の有限倍率だけでは表せません。
[ダスト電流の単位](dust_currents.md#定義と単位)も参照してください。

### Raycast 光電子

`use_raycast = .true.` と `nflag_emit(ispec) == 2` の組み合わせで有効化
されます。解釈: 内部表面に衝突したバックトレース粒子は、その面から放出
されたばかりの光電子である可能性がある。以下の 2 つの物理条件がゲート
となります。

**1. 外向き放出の半空間制約**: 光電子は表面から物質の外側へ向かって放出
されます。外向き法線 $\hat{\mathbf n}$ を太陽方向 (表面から太陽への単位
ベクトル) で近似すると、粒子速度は $\mathbf{v} \cdot \hat{\mathbf n} > 0$
を満たす必要があります。満たさなければ密度 0。

**2. 日射の遮蔽チェック**: 光電子が存在するためには、その表面に日射が
届いている必要があります。衝突点から $\hat{\mathbf n}$ 方向にレイを発射
し、遮蔽候補となる境界リスト (内部の表面やオブジェクト) と干渉させます。
$t > 0$ で何かに当たれば日陰となり密度 0。

両条件を満たす場合の確率は

$$
f_{\rm PE}(\mathbf{v}) = 2 \cdot C \prod_{i=1}^{3}
  \frac{1}{\sqrt{2\pi}\sigma_i}
  \exp\!\left(-\frac{u_i^2}{2 \sigma_i^2}\right),
$$

$\boldsymbol\mu$ = `vdri_vector(ispec)`, $\boldsymbol\sigma$ =
`[peth(ispec), peth(ispec), path(ispec)]`、$\mathbf u$ は上と同じ主軸への
逆回転です。係数 2 は $\mu_\parallel = 0$ のときの半空間正規化
として厳密です。シフト分布では 1 次近似となり、$|\boldsymbol\mu \cdot
\hat{\mathbf n}|$ が熱速度に対して支配的でない範囲で十分精確です。

太陽方向は
[`src/emses/emses_simulator_builder.f90`](../src/emses/emses_simulator_builder.f90)
の `resolve_sun_direction(ispec)` で決定します:

1. `ray_zenith_angle_deg(ispec) < 9000d0` ならそれを、そうでなければ
   `vdthz(ispec)` を使用。方位角も同様。
2. $[0, 0, 1]$ から開始し、$y$ 軸まわりに $-\zeta$、$z$ 軸まわりに
   $\phi$ だけ回転。
3. 得られた単位ベクトルの符号を反転して返す (放出面から太陽への向き;
   ユーザ指定の慣習)。

## 座標と単位の規約

- 位置は EMSES グリッド単位 ($[0, nx] \times [0, ny] \times [0, nz]$)。
- 速度は EMSES の自然単位系。スケールは `/intp/` の `vdri`, `path`,
  `peth` と整合。
- 境界コード `npbnd`: `0` 周期、`1` 反射、`2` 開放 (Maxwellian 源候補)、
  `3` 吸収。
- 粒子種インデックス `ispec` は **Fortran 側で 1-オリジン**、**Python
  API で 0-オリジン**。ラッパが自動変換。

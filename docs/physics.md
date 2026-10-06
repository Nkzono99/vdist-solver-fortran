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

## 対応するソース分布

### Zero (吸収)

`t_ZeroProbability%at` は常に 0。逆トレースを「吸い込み」にしたい壁に使います。

### シフト Maxwell 分布

`t_MaxwellianProbability%at` は 3D Maxwell–Boltzmann 密度:

$$
f_M(\mathbf{v}) = C \prod_{i=1}^{3} \frac{1}{\sqrt{2\pi}\sigma_i}
  \exp\!\left(-\frac{(v_i - \mu_i)^2}{2 \sigma_i^2}\right),
$$

$\boldsymbol\mu$ = `locs` はドリフト速度、$\boldsymbol\sigma$ = `scales`
は熱速度、$C$ = `coefficient` は乗算係数です。ビルダは EMSES パラメータ
から
[`src/emses/allcom.f90`](../src/emses/allcom.f90) の `vdri_vector(ispec)`
と `vth_vector(ispec)` を用いて値を供給します。

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
  \exp\!\left(-\frac{(v_i - \mu_i)^2}{2 \sigma_i^2}\right),
$$

$\boldsymbol\mu$ = `vdri_vector(ispec)`, $\boldsymbol\sigma$ =
`vth_vector(ispec)`。係数 2 は $\mu_\parallel = 0$ のときの半空間正規化
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

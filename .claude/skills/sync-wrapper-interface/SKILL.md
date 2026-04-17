---
name: sync-wrapper-interface
description: `vdsolverf/emses/wrapper.py`の`argtypes`とFortran側(`src/emses/emses_solver.f90`等)の`bind(c)`インタフェースが一致しているかを検証する。Fortran側の引数を変更した際に呼び出す。
user-invocable: true
---

# sync-wrapper-interface Skill

Python → Fortran のブリッジが壊れていないか、引数の型・数・順序を突き合わせて検証します。

## 対象

- Python 側: `vdsolverf/emses/wrapper.py` の各 `*_dll` 関数における `dll.<name>.argtypes = [...]`
- Fortran 側: `src/emses/emses_solver.f90` の `subroutine / function` で `bind(c, name="<name>")` 指定されたもの

## タスク

1. wrapper.py から `dll.<fname>.argtypes` と `dll.<fname>(...)` の呼び出しを Grep で列挙。
2. 各 `<fname>` について Fortran 側の `bind(c, name="<fname>")` を Grep で探し、引数リストを表にする。
3. 型対応表で突き合わせ:
   | Python (ctypes)                            | Fortran (iso_c_binding)            |
   |---|---|
   | `c_int` / `POINTER(c_int)`                 | `integer(c_int), value` / `integer(c_int)` |
   | `c_double`                                 | `real(c_double), value`            |
   | `c_char_p`                                 | `character(1, c_char)` 配列        |
   | `np.ctypeslib.ndpointer(dtype=np.float64, ndim=N)` | `real(c_double), intent(in/out/inout) :: x(...)` shape が一致 |
4. 不一致があれば該当行を `file:line` で提示し、修正案を示す (ただし確認後に適用)。

## 注意事項

- `value` 修飾子の有無は Python 側が値渡し / `POINTER` のどちらかと対応する。
- 配列の次元順 (C-contiguous / Fortran order) は numpy 側で `np.asfortranarray` 等を挟まない限り、Fortran で `(3*nspecies, lx+1, ly+1, lz+1)` のように宣言して整合させている点に注意。
- 引数追加時は wrapper.py の argtypes / 実呼び出し両方、Fortran 側の bind(c) 引数両方、そして必要なら Python の `c_int` 変換コードも忘れない。

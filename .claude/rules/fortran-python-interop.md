# Glob: vdsolverf/emses/wrapper.py, src/emses/*.f90

## Fortran ↔ Python I/F 整合ルール

`vdsolverf/emses/wrapper.py` は ctypes 経由で Fortran の共有ライブラリを呼び出します。以下を常に守ります。

### 変更手順 (引数追加/変更時)

1. **両側を同時に見る**
   - Fortran: `src/emses/emses_solver.f90` 等の `bind(c, name="<fname>")` 付きルーチン
   - Python: `wrapper.py` の `dll.<fname>.argtypes = [...]` と対応する `dll.<fname>(...)` 呼び出し
2. **引数リスト・順序・型を完全一致**させる。numpy の `ndim` と Fortran 側の配列次元数 (例: `real(c_double), intent(in) :: x(nx, ny, nz, 6)` なら `ndim=4`) を揃える。
3. **値渡し/参照渡しの区別**
   - Fortran `value` 修飾 ↔ Python スカラー (`c_int`, `c_double`)
   - Fortran 参照渡し ↔ Python `POINTER(c_int)` と `byref(...)`
4. **検証**
   - `sync-wrapper-interface` skill で突き合わせるか、Grep で `bind(c` と `argtypes` を比較する。
   - `.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"` で import だけでも確認。

### してはいけないこと

- wrapper.py の `argtypes` だけ変更し Fortran 側に反映しない (境界不整合でクラッシュする)。
- `restype` の関数名を別関数から流用する (直近のバグで `get_probabilities.restype` を `get_backtraces` 用に誤設定していた例がある。型は常に対応する関数のもの)。
- `mtd_vbnd` 等の配列インデックスを軸ごとに変えずにハードコード (`mtd_vbnd[0]` を X/Y/Z 全部に使うと壁面条件が誤る)。

### 公開 API 変更時の注意

- `get_backtrace` / `get_backtraces` / `get_probabilities` / `get_dust_backtrace` のキーワード引数名は既存利用者に影響するため、リネーム時はユーザに確認する (`AGENTS.md` 参照)。
- 例外: `get_dust_backtrace` の `os` → `system` への改名は他関数との整合を取るための意図的変更 (2026-04-17 実施)。

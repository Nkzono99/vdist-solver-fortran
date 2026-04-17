---
name: run-fortran-tests
description: `fpm test`でFortranのテストを実行する。`test/check.f90`のユニットテスト(ダスト帯電など)の回帰確認に使う。
user-invocable: true
---

# run-fortran-tests Skill

Fortran 側のユニットテストを `fpm test` で実行し、`test/check.f90` の各 `test_*` サブルーチンの成否を報告します。

## タスク

1. `fpm test 2>&1 | tail -40` を実行する。
2. "All tests passed." が表示されるか確認。
3. 失敗時は `FAILED:` 行の直後に出力される `actual` / `expected` 値を示し、該当テスト名と `src/core/dust_charge_simulator.f90` 等の該当関数を Grep して原因を説明する。

## 追加テストを書く場合

- `test/check.f90` の `contains` 節に `test_*` サブルーチンを追加し、`check` プログラム先頭の `call` リストにも追加する。
- 比較は `assert_close(label, actual, expected)` を利用 (許容誤差 1e-12)。
- 物理的に意味のある境界条件 (正/負のダスト電位、光電子電流の閾値など) を優先。

## 注意事項

- 長時間テストが生まれた場合は、`fpm test --target <name>` で個別実行できるようテスト名を実行形式 (`test/` 配下のプログラム) に分割することを検討する。

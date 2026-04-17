---
name: build-lib
description: Fortranコアをfpmでビルドし共有ライブラリ(.so/.dylib/.dll)をvdsolverf/配下に配置する。Fortran側のコードを変更した後にPythonから呼び直したいときに使用。
user-invocable: true
---

# build-lib Skill

Fortran を `fpm` でビルドし、Python ラッパーが参照する共有ライブラリを `vdsolverf/` に同期します。

## 前提

- `.venv/` が用意済みであること (必要なら `improve-harness` や開発ガイドを参照)。
- `fpm` が PATH 上にあること (`which fpm` で確認)。

## タスク

1. `fpm build --flag "-fPIC"` でビルド。
2. Linux 例:
   ```bash
   make
   ```
   (内部で `fpm install` → 共有ライブラリを `vdsolverf/libvdist-solver-fortran.so` に配置)
3. ライブラリが更新されたことを `ls -l vdsolverf/libvdist-solver-fortran.*` で確認。
4. Python 側 smoke test:
   ```bash
   .venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"
   ```

## 注意事項

- `Makefile` の対応する OS ターゲットが無い場合はユーザに相談する (macOS/Windows は未整備)。
- ビルドログは長くなるため `tail -30` など要約して報告。
- ビルド失敗時は `src/` のどのモジュールで落ちたかをエラーメッセージ先頭から特定する。

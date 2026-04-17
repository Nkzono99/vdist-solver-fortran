# Glob: **/*.py, pyproject.toml, setup.py

## Python 実行環境

このリポジトリの Python ラッパーは **Python 3.7 以降**を要求します (内部で `typing.Literal` 等を使用)。

### 必ず `.venv` を使う

- リポジトリ直下 `.venv/` に Python 3.12 の venv を配置しています。
- Python / pip 実行は以下のみを使用:
  - `.venv/bin/python`
  - `.venv/bin/pip`
- システムの `python3` は 3.6 のため `Literal` 未対応で import が失敗します。`python3 -c "import vdsolverf"` は使わない。

### セットアップ手順 (venv が無い場合)

```bash
/usr/bin/python3.12 -m venv .venv
.venv/bin/pip install --upgrade pip setuptools wheel
.venv/bin/pip install f90nml emout numpy scipy tqdm
```

`.venv/` は `.gitignore` に含まれている前提 (まだなら `improve-harness` で追加)。

### よく使うコマンド

- import 確認: `.venv/bin/python -c "from vdsolverf.emses import wrapper; print('ok')"`
- 構文チェックのみ: `.venv/bin/python -c "import ast; ast.parse(open('path.py').read())"`

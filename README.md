はじめに
======

sympy に pull request する前の雑多コードや、個人的メモ書き。これまでローカルで管理してきたが、さすがに辛くなってきたので。

To Do
=====

* `Function`の勉強

チートシート(というか備忘録)
============================

pytestの使い方
--------------

基本

```
pytest
```

失敗したテストだけ

```
pytest --lf
```

変数を表示

```
pytest --showlocals
```

動作環境の変更

```
set SYMPY_GROUND_TYPES=python
```

ログメッセージの表示

```
pytest --log-cli-level=DEBUG
```

カバレッジの取得

```
pytest -v --cov=. --cov-report=html
```

git rebase のしかた
-------------------

1. `master`を最新化
2. 自分のブランチに移動
3. `git rebase master`
4. プッシュ時にforceオプションを付ける

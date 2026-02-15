# GitHubリポジトリ登録手順

このパッケージをGitHubに登録する手順です。

## 方法1: GitHub CLIを使用（推奨）

### 1. GitHub CLIのインストール

macOSの場合:
```bash
brew install gh
```

### 2. GitHub CLIの認証

```bash
gh auth login
```

### 3. リポジトリの作成とプッシュ

```bash
cd /Volumes/Caches/Dropbox/gitbox/genice3-svg
./setup_github.sh
```

または手動で:
```bash
gh repo create vitroid/genice3-svg --public --source=. --remote=origin --push
```

## 方法2: 手動でセットアップ

### 1. GitHubでリポジトリを作成

1. https://github.com/new にアクセス
2. リポジトリ名: `genice3-svg`
3. Publicを選択
4. README、.gitignore、LICENSEは追加しない（既にローカルに存在）
5. "Create repository"をクリック

### 2. リモートを追加してプッシュ

```bash
cd /Volumes/Caches/Dropbox/gitbox/genice3-svg
git remote add origin https://github.com/vitroid/genice3-svg.git
git branch -M main
git push -u origin main
```

## 確認

リポジトリが正常に作成されたか確認:
```bash
git remote -v
```

ブラウザで https://github.com/vitroid/genice3-svg にアクセスして確認してください。





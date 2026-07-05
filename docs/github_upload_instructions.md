# GitHub upload instructions

## 1. Create a new repository on GitHub

Create a new repository named `lif-cell-counter`. Do not initialize it with a README, license, or `.gitignore`.

## 2. Initialize Git locally

```bash
cd lif-cell-counter
git init
git add .
git commit -m "Initial commit: LIF cell counting workflow"
```

## 3. Connect and push

HTTPS:

```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/lif-cell-counter.git
git push -u origin main
```

SSH:

```bash
git branch -M main
git remote add origin git@github.com:YOUR_USERNAME/lif-cell-counter.git
git push -u origin main
```

## 4. Verify after upload

```bash
git clone https://github.com/YOUR_USERNAME/lif-cell-counter.git
cd lif-cell-counter
python -m venv .venv
source .venv/bin/activate
pip install -e .
lif-cell-counter --help
```

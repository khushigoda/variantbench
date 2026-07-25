#!/usr/bin/env bash
# Step-2 preflight: verify LDSC is installed and reachable BEFORE the DAG hits it.
# Run:  bash scripts/preflight_ldsc.sh
set -u
ok=1
echo "=== LDSC preflight ==="

# 1. ldsc repo present?
if [ -f ~/ldsc/ldsc.py ] && [ -f ~/ldsc/munge_sumstats.py ]; then
  echo "  [ok]  ~/ldsc/ldsc.py and munge_sumstats.py found"
else
  echo "  [FAIL] ~/ldsc not found. Clone it:  git clone https://github.com/bulik/ldsc.git ~/ldsc"
  ok=0
fi

# 2. ldsc conda env present?
if conda env list 2>/dev/null | grep -qE '(^|/)ldsc[[:space:]]'; then
  echo "  [ok]  conda env 'ldsc' exists"
  # 3. can it import the LDSC deps (py2: numpy/scipy/pandas/bitarray)?
  if conda run -n ldsc python -c "import numpy,scipy,pandas,bitarray" 2>/dev/null; then
    echo "  [ok]  ldsc env imports numpy/scipy/pandas/bitarray"
  else
    echo "  [FAIL] ldsc env missing deps. Create per LDSC README:"
    echo "         conda create -n ldsc python=2.7 numpy scipy pandas bitarray -y"
    ok=0
  fi
else
  echo "  [FAIL] conda env 'ldsc' not found. Create it:"
  echo "         conda create -n ldsc python=2.7 numpy scipy pandas bitarray -y"
  ok=0
fi

# 4. reference files present?
for f in data/ref/w_hm3.snplist data/ref/eur_w_ld_chr/1.l2.ldscore.gz data/ref/eur_w_ld_chr/22.l2.M_5_50; do
  [ -f "$f" ] && echo "  [ok]  $f" || { echo "  [FAIL] missing $f"; ok=0; }
done

echo "======================"
[ "$ok" = 1 ] && echo "PREFLIGHT PASS — Step 2 is ready to run." || echo "PREFLIGHT FAIL — fix the items above before running Step 2."
exit $((1-ok))

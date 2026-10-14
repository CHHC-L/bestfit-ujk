#!/bin/zsh
# usage: run_vi_case.sh j k M mode [kelley args...]      (set REACH=reach to restrict to reachable patterns)
# Kelley/value-iteration search -> exact C++ verify -> independent Python verify -> gzip certificate;
# logs to results/search/logs/.  PYTHON (default: .venv/bin/python if present, else python3) needs numpy, highspy.
set -u
j=$1; k=$2; M=$3; mode=$4; shift 4
root=${0:A:h:h}
tag=${j}_${k}_${mode}_vi_M${M//,/-}${REACH:+_R}
lpd=$root/lp/$tag
cert=$root/certs/$tag.cert
log=$root/results/search/logs/$tag.log
mkdir -p $root/results/search/logs $root/lp
py=${PYTHON:-$root/.venv/bin/python}
[[ -x $py ]] || py=python3
{
  echo "== $tag  $(date)"
  $py $root/scripts/kelley.py $j $k $M $mode ${REACH:+--reach} "$@"
  $root/bfcert gen $j $k $M $lpd ${REACH:-} > /dev/null
  $root/bfcert verify $j $k $M $cert ${REACH:-}
  if [[ -z "${SKIP_INDEP:-}" ]]; then
    $py $root/scripts/verify_indep.py $j $k $M $lpd $cert ${REACH:-}
  fi
  rm -f $lpd/lp.bin
  [[ -f $cert ]] && gzip -f $cert
} > $log 2>&1
grep -E "^best c|final VI|extreme|VERDICT|independent" $log | sed "s/^/[$tag] /"

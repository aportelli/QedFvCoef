#!/usr/bin/env bash

set -euo pipefail

vel() {
  local x=$1
  local w=$2
  echo "scale=6; ep=e($w*$x); em=e(-$w*$x); ((ep-em)/(ep+em) + 0.2)/1.2" | bc -l
}

if (( $# != 7 )); then
    echo "usage: $(basename "$0") <j> <precision> <points> <length (seconds))> <fps> <slow down factor> <output directory>" 1>&2
    exit 1
fi
j=$1
prec=$2
points=$3
length=$4
fps=$5
w=$6
outdir=$7

mkdir -p "${outdir}"
frames=$((length*fps))
outar="c${j}_${points}pts_data.tar.bz2"
for f in $(seq 1 ${frames}); do
  x=$(echo "scale=10; ${f}/(${frames}+1)" | bc)
  v=$(vel "${x}" "${w}")
  out=$(printf "c${j}_${points}pts_f%06d_%.4f.dat" "${f}" "${v}")
  printf "Frame %4d/%d x= %8.6f v= %12.10f => '%s'\n" "${f}" "${frames}" "${x}" "${v}" "${out}"
  qed-fv-coef-spherescan "${j}" -v "${v}" -q r -e "${prec}" -n "${points}" -o "${outdir}/${out}" -g
done
cwd=$(pwd)
cd "${outdir}"
tar cjf "${outar}" ./*.dat
rm -f ./*.dat
cd "${cwd}"

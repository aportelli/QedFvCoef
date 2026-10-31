#!/usr/bin/env bash

set -euo pipefail

vel() {
  local x=$1
  local w=$2
  echo "scale=6; ep=e($w*$x); em=e(-$w*$x); ((ep-em)/(ep+em) + 0.2)/1.2" | bc -l
}

if (( $# != 8 )); then
    echo "usage: $(basename "$0") <j> <QED> <precision> <points> <length (seconds))> <fps> <slow down factor> <output directory>" 1>&2
    exit 1
fi
j=$1
qed=$2
prec=$3
points=$4
length=$5
fps=$6
w=$7
outdir=$8

mkdir -p "${outdir}"
frames=$((length*fps))
outar="c${j}_${points}pts_data.tar.bz2"
for f in $(seq 1 ${frames}); do
  x=$(echo "scale=10; ${f}/(${frames}+1)" | bc)
  v=$(vel "${x}" "${w}")
  out=$(printf "c${j}_${points}pts_f%06d_%.4f.dat" "${f}" "${v}")
  printf "Frame %4d/%d x= %8.6f v= %12.10f => '%s'\n" "${f}" "${frames}" "${x}" "${v}" "${out}"
  qed-fv-coef-spherescan "${j}" -v "${v}" -q "${qed}" -e "${prec}" -n "${points}" -o "${outdir}/${out}" -g
done
cwd=$(pwd)
cd "${outdir}"
tar cvjf "${outar}" ./*.dat
rm -f ./*.dat
cd "${cwd}"

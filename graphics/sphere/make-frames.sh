#!/usr/bin/env bash

set -euo pipefail

if (( $# != 1 )); then
    echo "usage: $(basename "$0") <data archive>" 1>&2
    exit 1
fi
ar=$1

cwd=$(pwd)
dir=$(dirname "${ar}")
ar=$(basename "${ar}")
outar=${ar//data/frames}
cd "${dir}"
tar xf "${ar}"
cd "${cwd}"
for file in "${dir}"/*.dat; do
  v=$(echo "${file%.dat}" | awk -F '_' '{print $4}')
  echo "${file} => ${file%.dat}.png"
  ./sphereplot.sh "${file}" "${file%.dat}.png" "${v}"
done
cd "${dir}"
tar cjf "${outar}" ./*.png
rm -f ./*.png ./*.dat
cd "${cwd}"

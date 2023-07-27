#!/usr/bin/env bash

set -euo pipefail

if (( $# != 1 )); then
    echo "usage: $(basename "$0") <dir>" 1>&2
    exit 1
fi
dir=$1

for file in "${dir}"/*.dat; do
  echo "${file} => ${file%.dat}.png"
  ./sphereplot.sh "${file}" "${file%.dat}.png"
done
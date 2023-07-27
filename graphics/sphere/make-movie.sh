#!/usr/bin/env bash

set -euo pipefail

if (( $# != 3 )); then
    echo "usage: $(basename "$0") <frame archive> <fps> <out>" 1>&2
    exit 1
fi
ar=$1
fps=$2
out=$3

cwd=$(pwd)
dir=$(dirname "${ar}")
ar=$(basename "${ar}")
cd "${dir}"
tar xf "${ar}"
cd "${cwd}"
ffmpeg -framerate "${fps}" -pattern_type glob -i "${dir}/*.png" -c:a copy -c:v libx264 -vb 20M -preset veryslow "${out}"
rm "${dir}"/*.png
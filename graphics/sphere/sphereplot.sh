#!/usr/bin/env bash

set -euo pipefail

if (( $# != 2 )); then
    echo "usage: $(basename "$0") <data> <output>" 1>&2
    exit 1
fi
data=$1
out=$2

gnuplot << EOF
file = "${data}"

set term png enhanced size 2048,2048 font "source sans pro,36.0"
set output "${out}"
set pm3d depthorder
set pm3d interpolate 2,2
set mapping spherical
set view equal xyz
unset border
unset xtics
unset ytics
unset ztics
set lmargin 0

stats file u 3 nooutput
set cbrange [STATS_min:STATS_max]
 set palette defined ( STATS_min 'red', 0 'white', STATS_max 'blue' )
splot file u 2:(\$1-pi/2):(1):3 notitle w pm3d
EOF

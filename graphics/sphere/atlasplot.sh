#!/usr/bin/env bash

set -euo pipefail

if (( ($# != 2) && ($# != 3) )); then
    echo "usage: $(basename "$0") <data> <output> [<v>]" 1>&2
    exit 1
fi
data=$1
out=$2
if (($# == 3)); then
    vcmd="set label '|v| = $3' at screen 0.8,0.8 font 'source sans pro,48.0'"
else
    vcmd=''
fi

gnuplot << EOF
file = "${data}"

set term png enhanced transparent size 2917,1947 font "source sans pro,36.0"
set output "${out}"
${vcmd}
set pm3d map
set pm3d interpolate 4,4
unset border
unset xtics
unset ytics
unset colorbox
set lmargin screen 0.
set rmargin screen 1.
set tmargin screen 1.
set bmargin screen 0.

load 'gnuplot-palettes/moreland.pal'

splot file u 2:1:3 image notitle

EOF

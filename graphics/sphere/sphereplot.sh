#!/usr/bin/env bash

set -euo pipefail

if (( ($# != 2) && ($# != 3) )); then
    echo "usage: $(basename "$0") <data> <output> [<v>]" 1>&2
    exit 1
fi
data=$1
out=$2
if (($# == 3)); then
    # vcmd="set label '|v| = $3' at screen 0.8,0.8 font 'source sans pro,48.0'"
    vcmd="set label '$3' at screen 0.85,0.85 font 'source sans pro,128.0'"
else
    vcmd=''
fi

gnuplot << EOF
file = "${data}"

set term png enhanced transparent size 2048,2048 font "source sans pro,36.0"
set output "${out}"
set pm3d depthorder
set pm3d interpolate 4,4
set mapping spherical
set view equal xyz
set view ,,1.6
${vcmd}
unset border
unset xtics
unset ytics
unset ztics
set lmargin 0
set rmargin 0
set tmargin 0
set bmargin 0
unset colorbox

stats file u (sin(\$1)*\$3*pi/2.) nooutput
sph_mean = STATS_mean
print "spherical average =", sph_mean, "+/-", STATS_stddev/sqrt(STATS_records)
stats file u 3 nooutput
set cbrange [STATS_min:STATS_max]
set palette defined ( STATS_min 'red', sph_mean 'white', STATS_max 'blue' )
splot file u 2:(\$1-pi/2):(1):3 notitle w pm3d
EOF

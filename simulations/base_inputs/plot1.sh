#!/bin/bash

# Bin to ASCII conversion
~/hGAST/bin2ascii.out < bin2ascii_mikro.inp

# Launch gnuplot - Rotor speed vs. Pitch angle
gnuplot -persist <<EOF
set terminal x11 size 1600,500
set key off
set multiplot layout 1,2 

set title "1: Rotor speed vs. Pitch angle"
set xlabel "Time [s]"
set ylabel "Rotor Speed [rad/s]"
set y2label "Pitch Angle [deg]"
set y2tics
set ytics nomirror
set key right bottom
set xtics 0,60
set mxtics 6
set grid xtics mxtics ytics
plot \
  'qsdof.dat' using 1:3 with lines linecolor rgb "red" title 'Rotor Speed' axes x1y1, \
  '' using 1:5 with lines linecolor rgb "blue" title 'Pitch Angle' axes x1y2

set title "2: Torque Demand vs. Mechanical Power"
set xlabel "Time [s]"
set ylabel "Torque [kNm]"
set y2label "Power [kW]"
set y2tics
set ytics nomirror
set key right bottom
set xtics 0,60
set mxtics 6
set grid xtics mxtics ytics
plot \
  'qsdof.dat' using 1:17 with lines linecolor rgb "red" title 'Torque Demand' axes x1y1, \
  '' using 1:18 with lines linecolor rgb "green" title 'Mechanical Power' axes x1y2

unset multiplot
EOF


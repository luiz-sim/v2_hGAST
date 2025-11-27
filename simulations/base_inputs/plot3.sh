#!/bin/bash
# plot.sh: compute and plot FFT of a specified DOF from qsdof3.dat
# Usage: ./plot.sh [dof_index]
# dof_index: 2=surge,3=sway,4=heave,5=roll,6=pitch,7=yaw (default 6: pitch)

DOF=${1:-6}

# Run FFT
~/fft.py $DOF

# Plot with gnuplot
gnuplot -persist <<EOF
set grid
set xrange [0:0.1]
set xlabel "Frequency [Hz]"
set ylabel "Amplitude"
set xtics 0.01
set title "FFT of DOF $DOF"
plot "fft_dof${DOF}.dat" using 1:2 with lines lw 2 title "DOF $DOF"
EOF


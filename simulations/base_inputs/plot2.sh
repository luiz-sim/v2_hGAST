#!/bin/bash
# Bin to ASCII conversion
~/hGAST/bin2ascii.out < bin2ascii_mikro.inp

# Launch gnuplot - 6 DOFs in subplots
gnuplot -persist <<EOF
set terminal x11 size 1600,500
set key off
set grid

# Layout: 3 columns, 2 rows
set multiplot layout 2,3 rowsfirst title "Floater Motion per DOF"

# Plot 1: Surge (x)
set xlabel "Time [s]"
set ylabel "Surge [m]"
plot 'qsdof3.dat' using 1:2 with lines lc rgb "blue" title "Surge"

# Plot 2: Sway (y)
set xlabel "Time [s]"
set ylabel "Sway [m]"
plot 'qsdof3.dat' using 1:3 with lines lc rgb "green" title "Sway"

# Plot 3: Heave (z)
set xlabel "Time [s]"
set ylabel "Heave [m]"
plot 'qsdof3.dat' using 1:4 with lines lc rgb "cyan" title "Heave"

# Plot 4: Roll (x)
set xlabel "Time [s]"
set ylabel "Roll [deg]"
plot 'qsdof3.dat' using 1:5 with lines lc rgb "red" title "Roll"

# Plot 5: Pitch (y)
set xlabel "Time [s]"
set ylabel "Pitch [deg]"
plot 'qsdof3.dat' using 1:6 with lines lc rgb "magenta" title "Pitch"

# Plot 6: Yaw (z)
set xlabel "Time [s]"
set ylabel "Yaw [deg]"
plot 'qsdof3.dat' using 1:7 with lines lc rgb "orange" title "Yaw"

unset multiplot
EOF


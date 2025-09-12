# Gnuplot script to plot Alpha vs Cl
set terminal pdfcairo enhanced color solid
set output "alpha_vs_Cl.pdf"

# Set title and labels
set title "Alpha vs Cl"
set xlabel "Alpha (degrees)"
set ylabel "Cl (Lift Coefficient)"
set grid

# Data formatting
set datafile separator whitespace

# Line styles
set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1.8
set style line 2 linecolor rgb 'red' linetype 6 linewidth 2 ps 1.8
set style line 3 linecolor rgb 'blue' linetype 6 linewidth 2 ps 1.8
set style line 4 linecolor rgb 'green' linetype 6 linewidth 2 ps 1.8

# Legend and axes
set key top right box
set xrange [-3:3]
set yrange [-0.4:0.6]
set xtic nomirror
set ytic nomirror
set mxtic 4
set mytic 4

# Define pi for sin
pi = 3.14159265358979

# Plot
plot \
  'reference_data.dat' using 5:7 with points pt 6 ps 1.5 lc rgb "blue" title 'Experiment', \
  'forces/0/forceCoeffs.dat' every ::1 using \
    (alpha = 0.016 - 2.51 * sin(39 * $1)) : \
    ($3 * cos(alpha * pi / 180.0) + $2 * sin(alpha * pi / 180.0)) \
    with lines ls 3 title "dbns"
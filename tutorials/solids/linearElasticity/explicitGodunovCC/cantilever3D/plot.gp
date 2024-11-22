# Gnuplot script to plot Time vs Dy
set terminal pngcairo size 1024,768
set output 'time_vs_dy.png'

# Set labels
set xlabel 'Time (s)'
set ylabel 'Dy (m)'
set title 'Time vs Dy Displacement'

# Style and grid
set grid
set style data linespoints

# Plotting the data
plot 'postProcessing/0/solidPointDisplacement_displacement.dat' using 1:3 with linespoints title 'Dy Displacement'

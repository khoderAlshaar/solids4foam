# Gnuplot script to plot Time vs Dx
set terminal pdfcairo enhanced color solid

set output "alpha_vs_Cm_gn.pdf"

# Set title and labels
set title "Alpha vs Cm"
set xlabel "Alpha (degrees)"
set ylabel "Cm"
set grid

# Set data style
set datafile separator whitespace

set style line 1 linecolor rgb 'black' linetype 6 linewidth 2 ps 1.8
set style line 2 linecolor rgb 'red' linetype 6 linewidth 2 ps 1.8
set style line 3 linecolor rgb 'blue' linetype 6 linewidth 2 ps 1.8
set style line 4 linecolor rgb 'green' linetype 6 linewidth 2 ps 1.8
set key t r
set key box
set xrange [-3:3]
set yrange [-0.02:0.02]
set xtic nomirror
set ytic nomirror
set mxtic 4
set mytic 4



# Plot as unconnected open circles
plot 'reference_data.dat' using 5:3  with points pt 6 ps 1.5 lc rgb "blue" title 'Exprement'\
    ,"forces/0/forceCoeffs.dat" every ::1 u (0.016 + 2.51 * sin(39.971 * $1)) : ($4) w l ls 3 title "dbns" \






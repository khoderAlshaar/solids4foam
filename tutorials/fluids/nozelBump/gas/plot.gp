# Gnuplot script to plot Time vs Dx
set terminal pdfcairo enhanced color solid
set output "Ma_0.1_cp.pdf"

#set terminal svg size 800,600 font "Arial,12"
#set output "PxxVsTime.svg"

# Set labels
set xlabel 'X'
set ylabel 'Cp'
#set title 'Ma = 0.1'

set style line 1 linecolor rgb 'black' linetype  6 linewidth 1 ps 0.5
set style line 2 linecolor rgb 'red' linetype 6 linewidth 2 ps 1.8
set style line 3 linecolor rgb 'blue' linetype 6 linewidth 2 ps 1.8
set style line 4 linecolor rgb 'green' linetype 6 linewidth 2 ps 1.8
set key t r
set key box

set grid

set xtic nomirror
set ytic nomirror
set mxtic 4
set mytic 4

# Plotting the data
plot "bottomData.csv" u ($1):($3) w l ls 3 title "HGrid" \
   


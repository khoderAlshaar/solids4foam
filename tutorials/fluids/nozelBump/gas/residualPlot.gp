set datafile separator ","

set terminal pngcairo size 1200,800 font "Arial,14"
set output "errors_log.png"

set title "L2 Error History (Log Scale)"
set xlabel "Time Index"
set ylabel "Error Value (log scale)"

set logscale y
set format y "10^{%L}"
set grid
set key left top

plot \
    "physicalResiduals.csv" using 2:3 with linespoints lw 2 pt 7 title "rho", \
    "physicalResiduals.csv" using 2:4 with linespoints lw 2 pt 7 title "rhoU", \
    "physicalResiduals.csv" using 2:5 with linespoints lw 2 pt 7 title "rhoE"

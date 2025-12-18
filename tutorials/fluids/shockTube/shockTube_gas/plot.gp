# ===================================================
# Enhanced Publication-Quality 2x2 Shock Tube Plot
# ===================================================

set terminal pdfcairo enhanced color solid font "Arial,12" size 7in,8in
set output "shockTube_4plots_enhanced.pdf"

# Global styling
set grid linetype 1 linecolor rgb "#cccccc"
set tics font "Arial,10"
set key bottom left font "Arial,10"


# Define line styles
set style line 1 lc rgb "black" lt 1 lw 2    # Exact solution line
set style line 2 lc rgb "red"   lt 1 lw 2 pt 7 ps 0.7   # Numerical points

# Margins tighter for better space usage
set multiplot layout 2,2 title "Shock Tube Results at t = 0.00075 s" font "Arial,14"

# ===================================================
# 1) PRESSURE
# ===================================================
set xlabel "x (m)"
set ylabel "Pressure (MPa)"
set title "Pressure Profile" font "Arial,12"

plot "exact_solution.dat" \
        using 1:($2/1e6) with lines ls 1 title "Exact P", \
     "postProcessing/sets/0.00075/data_p_rho_T.xy" \
        using 1:($2/1e6) with lines ls 2 title "Numerical P"

# ===================================================
# 2) TEMPERATURE
# ===================================================
set xlabel "x (m)"
set ylabel "Temperature (K)"
set title "Temperature Profile" font "Arial,12"

plot "exact_solution.dat" \
        using 1:4 with lines ls 1 title "Exact T", \
     "postProcessing/sets/0.00075/data_p_rho_T.xy" \
        using 1:4 with lines ls 2 title "Numerical T"

# ===================================================
# 3) DENSITY
# ===================================================
set xlabel "x (m)"
set ylabel "Density (kg/m³)"
set title "Density Profile" font "Arial,12"

plot "exact_solution.dat" \
        using 1:3 with lines ls 1 title "Exact ρ", \
     "postProcessing/sets/0.00075/data_p_rho_T.xy" \
        using 1:3 with lines ls 2 title "Numerical ρ"

# ===================================================
# 4) VELOCITY
# ===================================================
set xlabel "x (m)"
set ylabel "Velocity Ux (m/s)"
set title "Velocity Profile" font "Arial,12"

plot "exact_solution.dat" \
        using 1:5 with lines ls 1 title "Exact Ux", \
     "postProcessing/sets/0.00075/data_U.xy" \
        using 1:2 with lines ls 2 title "Numerical Ux"

unset multiplot

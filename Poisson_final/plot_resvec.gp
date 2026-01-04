### plot_resvec.gp — convergence plot (final report version)

reset
set encoding utf8

# ---------- Inputs ----------
f_alpha = "results/RESVEC_alpha.dat"
f_jac   = "results/RESVEC_jac.dat"
f_gs    = "results/RESVEC_gs.dat"

tol = 1e-3

# Auto-detect nbite from number of records (lines)
stats f_alpha nooutput
nb_alpha = STATS_records
stats f_jac nooutput
nb_jac   = STATS_records
stats f_gs nooutput
nb_gs    = STATS_records

# ---------- Common style ----------
set title "Convergence — Richardson (alpha / Jacobi / Gauss-Seidel)"
set xlabel "Itération k"
set ylabel "Résidu relatif ||r_k||_2 / ||b||_2"

set logscale y
set grid xtics ytics
set key top left

set xrange [0:*]
set yrange [1e-3:1]

# tol horizontal line
set arrow 1 from graph 0, first tol to graph 1, first tol nohead lw 1 dt 2
set label 1 sprintf("tol = %.0e", tol) at graph 0.02, first (tol*1.05)

# Vertical lines at nbite (optional but nice for report)
set arrow 2 from first nb_alpha, graph 0 to first nb_alpha, graph 1 nohead lw 1 dt 3
set arrow 3 from first nb_gs,    graph 0 to first nb_gs,    graph 1 nohead lw 1 dt 3

# nbite labels — avoid overlap by separating y and shifting x a bit
# (alpha and jacobi usually end at same maxit, so give different y)
x_right = nb_alpha - 1200     # 更向左
x_gs    = nb_gs + 150            # shift right a little (still inside plot)


set label 10 sprintf("nbite(alpha)  = %d", nb_alpha) at first x_right, first (tol*3.20)   
set label 11 sprintf("nbite(Jacobi) = %d", nb_jac)   at first x_right, first (tol*2.40)   

set label 12 sprintf("nbite(GS)     = %d", nb_gs)    at first x_gs,    first (tol*1.15)

# ---------- Plot command (wrapped as a macro) ----------
PLOT_CMD = sprintf("plot \
  '%s' using 0:1 with lines lw 2 title 'Alpha optimal', \
  '%s' using 0:1 with lines lw 2 title 'Préconditionnement Jacobi', \
  '%s' using 0:1 with lines lw 2 title 'Préconditionnement Gauss-Seidel'", \
  f_alpha, f_jac, f_gs)

# ---------- Output PNG ----------
set terminal pngcairo size 1400,900 enhanced font "Sans,18"
set output "results/convergence_resvec.png"
eval(PLOT_CMD)
unset output

# ---------- Output PDF ----------
set terminal pdfcairo size 16cm,10cm enhanced font "Sans,12"
set output "results/convergence_resvec.pdf"
replot
unset output

print sprintf("Wrote: results/convergence_resvec.png and results/convergence_resvec.pdf")

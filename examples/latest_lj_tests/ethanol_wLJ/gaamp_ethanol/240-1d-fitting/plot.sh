set terminal png medium font "arial,26" size 1200,800
set title "Dihedral  1 - 2 - 3 - 9"
set xlabel "phi-1" font "arial,26"
set ylabel "E" font "arial,26"
set output "plot-1d-1.png"
plot "fitting-1d-1.dat" using 1:2 t "QM" with linespoints pt 6 lw 2 ps 2 lc rgb "black", "fitting-1d-1.dat" using 1:3 t "Fitted" with linespoints pt 1 lw 2 ps 2 lc 1, "org-1d-1.dat" using 1:3 t "org" with linespoints pt 2 lw 2 ps 2 lc 2

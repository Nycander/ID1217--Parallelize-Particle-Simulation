set xlabel "Number of Particles"
set ylabel "Time (s)"
set terminal pdfcairo #enhanced "Helvetica" 16


set log y
plot "plot_dir/openmp.plt" using 1:2 with lines\
    title "openmp"

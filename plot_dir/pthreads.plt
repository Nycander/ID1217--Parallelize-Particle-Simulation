set xlabel "Number of Particles"
set ylabel "Time (s)"
set terminal pdfcairo #enhanced "Helvetica" 16


set log y
set log x
plot "plot_dir/pthreads.dat" using 1:2 with lines\
    title "pthreads"

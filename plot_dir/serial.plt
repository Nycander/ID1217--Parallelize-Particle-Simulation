set xlabel "Number of Particles"
set ylabel "Time (s)"
set terminal pdfcairo #enhanced "Helvetica" 16


set log y
plot "plot_dir/serial.dat" using 1:2 with lines\
    title "serial"

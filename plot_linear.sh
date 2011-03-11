#This script sure is manly.
#Runs: mpi, openmp, pthreads and serial variants of the particle simulation.

PROGRAMS="serial
pthreads
openmp"
#mpi"


mkdir -p plot_dir
mkdir -p report/plots
rm -f plot_dir/*.dat
rm -f plot_dir/*.tmp

for p in $PROGRAMS
do
    echo "The Frog King is brave ["$p"]"
    outfile="plot_dir/"${p}"_out.dat"
    datfile="plot_dir/"${p}".dat"
    tmpfile="plot_dir/"${p}".tmp"
    for i in `seq 10 10 1000`
    do
        echo "The Frog King is brave  as can be ["$i"]"
        for j in `seq 1 1 5`
        do
            $p/$p -n $i -p 4 | tail -n 1 | ruby sed.rb >> $tmpfile
        done
        cat $tmpfile | ruby median.rb >> $outfile
        echo "" > $tmpfile
    done
    echo "The Frog King is braver than me."
    cat $outfile | ruby plotit.rb > $datfile
    echo "The Frog King is braver as you can see."
    echo "The Frog King is braver than you."
    gnuplot "plot_dir/"$p".plt" > "report/plots/"$p".pdf"
done

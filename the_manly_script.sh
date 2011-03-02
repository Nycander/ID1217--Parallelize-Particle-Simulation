#This script sure is manly.
#Runs: mpi, openmp, pthreads and serial variants of the particle simulation.
#export

PROGRAMS="serial
openmp
pthreads"
#mpi"

test[0]=1
test[1]=2
test[2]=4
test[3]=8
test[4]=16
test[5]=32
test[6]=64
test[7]=128
test[8]=256
test[9]=512
test[10]=1024
test[11]=2048
test[12]=4096
test[13]=8192

mkdir -p plot_dir
rm -f plot_dir/*

for p in $PROGRAMS
do
    for i in `seq 0 1 13`
    do
        #echo $p $i
        outfile="plot_dir/"${p}"_performance"
        $p -n $test[i] | tail -n 1 | sed '/\d\.\d+/p' >> $outfile
    done
done

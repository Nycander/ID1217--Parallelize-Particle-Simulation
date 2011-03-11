$PROGRAMS = ["serial", "pthreads", "openmp"]#, "mpi" ]
$SPEED_UP = ["pthreads", "openmp"]#, "mpi" ]
$NUMBER_OF_CORES = 4
$NUMBER_OF_REPETED_TESTS = 5
$PARTICLES_TO_SIMULATE = 2000
$RESULT = []
$one_thread

def median(floats)
    floats.sort!
    return floats[floats.length/2]
end

def main()
    $PROGRAMS.each_index do |p_i|
        p = $PROGRAMS[p_i]
        $RESULT[p_i] = []
        puts p
        1.upto($NUMBER_OF_CORES) do |c|
            time = []
            1.upto($NUMBER_OF_REPETED_TESTS) do |j|
                time[j] = `#{p}/#{p} -n #{$PARTICLES_TO_SIMULATE} -p #{c}`.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
                puts "cores: " + c.to_s + "   time: " + time[j].to_s + "   speedup: " + ($one_thread / time[j]).to_s unless c == 1
            end
            time.compact!()
            $RESULT[p_i][c] = median(time)
            $one_thread = $RESULT[p_i][1]
            puts "median:"
            puts "cores: " + c.to_s + "   time: " + $RESULT[p_i][c].to_s + "   speedup: " + ($RESULT[p_i][1] / $RESULT[p_i][c]).to_s
        end
        puts "\n"

    end
    $PROGRAMS.each_index() do |p_i|
        file = "plot_dir/" + $PROGRAMS[p_i] + "_speedup.dat"
        p = $PROGRAMS[p_i]
        File.delete(file) if File.exists?(file)
        next if p_i == 0
        1.upto($NUMBER_OF_CORES) do |c|
            puts "p_i= " + p_i.to_s  + "   c= " + c.to_s
            speedup = $RESULT[p_i][1]/$RESULT[p_i][c]
            open(file, 'a') do |f|
                f.puts c.to_s + "   " + speedup.to_s
            end
        end
        gnuplot(p)
    end
end
def gnuplot(p)
    `echo 'set xlabel "Number of Threads"
set ylabel "Speedup Factor"
set terminal pdfcairo #enhanced "Helvetica" 16
set boxwidth 0.8
set style fill solid

set xrange [ 0.5 : #{$NUMBER_OF_CORES+0.5} ]
set yrange [ 0 : #{$NUMBER_OF_CORES} ]
set mxtics 0
set mytics 2
set xtics 0,1
set ytics 0,1


plot "plot_dir/#{p}_speedup.dat" using 1:2 with boxes lt rgb "steelblue"\
    title "#{p} speedup factor"' | gnuplot > report/plots/#{p}_speedup.pdf`
end
puts ARGV
if(ARGV[0] == "-ga")
$SPEED_UP.each() { |p| gnuplot(p)}
else
    main()
end


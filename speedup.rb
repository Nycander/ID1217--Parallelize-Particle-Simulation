programs = ["serial", "pthreads", "openmp"]#, "mpi" ]
number_of_cores = 4
number_of_repeted_tests = 5
particles_to_simulate = 2000
result = []
$one_thread

def median(floats)
    floats.sort!
    return floats[floats.length/2]
end

programs.each_index do |p_i|
    p = programs[p_i]
    result[p_i] = []
    puts p
    1.upto(number_of_cores) do |c|
        time = []
        1.upto(number_of_repeted_tests) do |j|
            time[j] = `#{p}/#{p} -n #{particles_to_simulate} -p #{c}`.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
            puts "cores: " + c.to_s + "   time: " + time[j].to_s + "   speedup: " + ($one_thread / time[j]).to_s unless c == 1
        end
        time.compact!()
        result[p_i][c] = median(time)
        $one_thread = result[p_i][1]
        puts "median:"
        puts "cores: " + c.to_s + "   time: " + result[p_i][c].to_s + "   speedup: " + (result[p_i][1] / result[p_i][c]).to_s
    end
    puts "\n"

end
programs.each_index() do |p_i|
    file = "plot_dir/" + programs[p_i] + "_speedup.dat"
    p = programs[p_i]
    File.delete(file) if File.exists?(file)
    next if p_i == 0
    1.upto(number_of_cores) do |c|
        puts "p_i= " + p_i.to_s  + "   c= " + c.to_s
        speedup = result[p_i][1]/result[p_i][c]
        open(file, 'a') do |f|
            f.puts c.to_s + "   " + speedup.to_s
        end
    end
`echo 'set xlabel "Number of Threads"
set ylabel "Speedup Factor"
set terminal pdfcairo #enhanced "Helvetica" 16
set boxwidth 0.8
set style fill solid

set xrange [ 0.5 : #{number_of_cores+0.5} ]
set yrange [ 0 : #{number_of_cores+0.5} ]
set mxtics 0
set mytics 2
set xtics 0,1
set ytics 0,1


plot "plot_dir/#{p}_speedup.dat" using 1:2 with boxes lt rgb "steelblue"\
    title "#{p} speedup factor"' | gnuplot > report/plots/#{p}_speedup.pdf`
end

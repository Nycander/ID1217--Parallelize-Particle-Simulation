
(1..8).each do | i |
	times = []
	reading = []
	writing = []

	(1..10).each do | j |
		 out = `./openmp -n 4000 -p #{i}`
		 times << out.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
	end

	times.sort!

	print times[times.length/2].to_s.ljust(8)
	print "\n"
end

(100..4000).step(100) do | i |
	times = []
	reading = []
	writing = []

	(1..10).each do | j |
		 out = `./serial -n #{i}`
		 times << out.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
	end

	times.sort!

	print i
	print "\t"
	print times[times.length/2].to_s.ljust(8)
	print "\n"
end
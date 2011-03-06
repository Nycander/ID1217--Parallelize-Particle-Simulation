
(1..8).each do | i |
	times = []
	(1..5).each do | j |
		 out = `pthreads.exe -n 2000 -p #{i}`
		 times << out.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
	end

	times.sort!

	print times[times.length/2].to_s.ljust(8), "\n"
end

(1..8).each do | i |
	times = []
	reading = []
	writing = []

	(1..5).each do | j |
		 out = `pthreads.exe -n 8000 -p #{i}`
		 times << out.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
		 #reading << out.gsub(/^.*?Reading: (\d+\.\d+) s.+?/m, "\\1").to_f
		 #writing << out.gsub(/^.*?Writing: (\d+\.\d+) s.+?/m, "\\1").to_f
	end

	times.sort!
	#reading.sort!
	#writing.sort!

	print times[times.length/2].to_s.ljust(8)
	#print "\t"
	#print reading[reading.length/2].to_s.ljust(8)
	#print "\t"
	#print writing[writing.length/2].to_s.ljust(8)
	print "\n"
end
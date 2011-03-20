
(1..8).each do | i |
	times = []
	force = []
	move = []
	send = []

	(1..10).each do | j |
		 out = `mpiexec -np #{i} ./mpi -n 10000`
		 times << out.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
		 force << out.gsub(/^.*?Force: (\d+\.\d+).+?$/m, "\\1").to_f
		 move << out.gsub(/^.*?Move: (\d+\.\d+).+?$/m, "\\1").to_f
		 send << out.gsub(/^.*?Communication: (\d+\.\d+).+?$/m, "\\1").to_f
	end

	times.sort!
	force.sort!
	move.sort!
	send.sort!

	print times[times.length/2].to_s
	print ","
	print force[force.length/2].to_s
	print ","
	print move[move.length/2].to_s
	print ","
	print send[send.length/2].to_s
	print "\n"
end
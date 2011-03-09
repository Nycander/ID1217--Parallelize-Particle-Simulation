
(1..8).each do | i |
	times = []
	force = []
	move = []
	send = []
	clear = []
	receive = []

	(1..10).each do | j |
		 out = `mpiexec -np #{i} ./mpi -n 4000`
		 times << out.gsub(/^.*?simulation time = (\d+\.\d+) seconds.+?$/m, "\\1").to_f
		 force << out.gsub(/^.*?Force: (\d+\.\d+).+?$/m, "\\1").to_f
		 move << out.gsub(/^.*?Move: (\d+\.\d+).+?$/m, "\\1").to_f
		 send << out.gsub(/^.*?Send: (\d+\.\d+).+?$/m, "\\1").to_f
		 clear << out.gsub(/^.*?Clear: (\d+\.\d+).+?$/m, "\\1").to_f
		 receive << out.gsub(/^.*?Receive: (\d+\.\d+).+?$/m, "\\1").to_f
	end

	times.sort!
	force.sort!
	move.sort!
	send.sort!
	clear.sort!
	receive.sort!

	print times[times.length/2].to_s
	print ","
	print force[force.length/2].to_s
	print ","
	print move[move.length/2].to_s
	print ","
	print send[send.length/2].to_s
	print ","
	print clear[clear.length/2].to_s
	print ","
	print receive[receive.length/2].to_s
	print "\n"
end
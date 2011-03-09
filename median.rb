floats = []
ARGF.each do | a |
    next if a.strip ==''
	f = a.strip.to_f
	if f == 0.0
		$stderr.puts "FU! Invalid argument supplied. '#{a.trim}' is not a float!"
		exit
	end
	floats << f
end

floats.sort!

puts floats[floats.length/2]

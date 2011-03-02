floats = []
ARGF.each do | a |
	if a.to_f == 0.0
		$stderr.puts "FU! Invalid argument supplied. '#{a}' is not a float!"
		exit
	end
	floats << a.to_f
end

floats.sort!

puts floats[floats.length/2]
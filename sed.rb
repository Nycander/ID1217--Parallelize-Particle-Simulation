ARGF.each_with_index do |line, idx|
    print line.gsub(/^.*(\d\.\d+).*$/, "\\1")
end
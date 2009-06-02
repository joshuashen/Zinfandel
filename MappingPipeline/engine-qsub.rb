## write qsub scripts

cmd = ARGV[0..-1].join("\s")

puts '#!/bin/bash'
puts '#$ -cwd'
puts cmd



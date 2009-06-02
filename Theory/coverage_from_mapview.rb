# compute coverage from mapview file

# coverage is defined as number of reads in a given region

step = 10000 
e = step
c = 0
while line = ARGF.gets do 
  cols = line.split(/\s+/)
  pos = cols[2].to_i
  if pos < e 
    c += 1
  else  # a new bin
    puts c
    c = 1
    e = e + step
  end
end
puts c 


    

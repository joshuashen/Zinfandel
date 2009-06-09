## produce a list of random CNVs

## format: 
#  VariationID     Landmark        Chr     Start   End     VariationType  
# Variation_10631 chr1:102600719..102601411       chr1    102600719       102601411       InDel  

fa = ARGV[0]
numCNVs = ARGV[1].to_i
maxSize = 20000
minSize = 100
step = 500  

chrs = {}
name = ''
total = 0
File.new(fa,'r').each do |line|
  line.chomp!
  if line=~ /^>(\S+)/
    name = $1
    chrs[name] = 0
  else
    chrs[name] += line.size
    total += line.size
  end
end

bins = ( maxSize - minSize ) / step
cnvSize = []
eachBin =  numCNVs / bins 
i = minSize
$stderr.puts "bins: #{bins}\teachBin:  #{eachBin}"

while i < maxSize 
  a = Array.new(eachBin*2)
  a.fill(i)
  cnvSize.concat(a)
  i += step
end
cnvSize = cnvSize.sort_by {rand()}
$stderr.puts "size of cnvSize array:  #{cnvSize.size}"

cnvStart = Array.new(numCNVs)
# 0.upto(numCNVs - 1) {|i| cnvStart[i] =  rand(total)}
cnvStart.map! {|i| rand(total)}

s,e = 0,0
idKey = 0
chrs.keys.sort.each do |ref|
  $stderr.puts "#{ref}\t#{chrs[ref]}"
  e = s + chrs[ref] - 1
  
  starts = cnvStart.select {|i| i >=s and i < e  }
  
  starts.sort.each do |ss|
    ll = cnvSize.shift
    ee = ll + ss - 1
    if ee <= e
      puts "#{idKey}\t#{ref}_#{ss}_#{ee}\t#{ref}\t#{ss}\t#{ee}\tindel\t#{ll}"
    end
    idKey += 1
  end
  s = e + 1

end
  


    

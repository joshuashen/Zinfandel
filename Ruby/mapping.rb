class Mapping
  attr_accessor :chromosomes

  def initialize(mapview, ref)
    @chromosomes = Hash.new
   
    name = ''
    seqLength = 0
    File.new(ref,'r').each do |line|
      line.chomp!
      if line=~ /^>(\S+)/
        if name != ''  # the end of the last ref
          @chromosomes[name] = Chromosome.new(name, seqLength)
        end
        seqLength =  0 
        name = $1
      else
        seqLength += line.chomp.size
      end
    end
    if name != ''  # the end of the last ref
      @chromosomes[name] = Chromosome.new(name, seqLength)
    end
   
    $stderr.puts "Fasta file #{ref} loaded."
 ## mapview format example: 
#  chr18_d_9_151_a0ac8/1	chr18	9	+	143	18	...
   
    File.new(mapview, 'r').each do |line|
      cols = line.split(/\s+/)
      read,ref,s,dir,dist,flag = cols[0], cols[1],cols[2].to_i, cols[3],cols[4].to_i,cols[5]
      
 #     if !@chromosomes[ref].starts.key?(s)
      @chromosomes[ref].coverage[s-1] += 1
#      else
#        @chromosomes[ref].starts[s] += 1
#      end
      
  #    if dir == '+' && flag == '18'
  #      @insertSize[ref][s] = dist
  #    end
    end
  end

end

class Chromosome
  attr_accessor  :name, :coverage, :insertSize, :size
  
  def initialize(name, l)
    @name = name
  #  @seq = ''
    @size = l
   
   # this implementation would take a lot of memory for complex genomes
    @coverage = Array.new(@size)
    @coverage.fill(0)
    @insertSize = {}
  end
end
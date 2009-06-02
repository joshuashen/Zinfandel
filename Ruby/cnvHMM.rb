### 

## generalized hidden markov model:
## only trigger break points if  P(insertSize|normal) ** clone_coverage < 0.05 
#

# the model:
#                  +----------------+   +-----------+
#                  |                |   |           |
##    [del_w/o_breakPoints]  <--  [normal] -->  dup_w/o_breakPoints
#                                  |     |
#                 {del_w_BP}  <-- {BP}  {BP} --> {dup_w_BP}
#
## [ ] states: duration are exponentially distributed
## { } states: duration are normal or Poisson distribution, conditioned on the BreakPoint entry, ie.  duration = distance - average_insert_size  
## positive:  deletion;  negative: insertion/duplication


require File.dirname(__FILE__) + '/mapping.rb'

class CnvHMM
 	attr_accessor :avgCov
    
  def initialize(model)  
    @states = []
    @transitionMatrix = Hash.new {|h,k| h[k] = Hash.new} # two-dimensional hash
    @emissionMatrix = Hash.new {|h,k| h[k] = Hash.new}
    @maxEmission = 100  # maximum emission
    @avgCov = 1.0
    @genome = 'd'  # diploid
    @initProb = {}
    loadModel(model)

#    @chromosomes = mapping.chromosomes # class Chromosome
    # !! need better prior !!
  end

  def loadModel(model)
  # format:
  ##  genome  haploid/diploid
  # 
  # 
    numCNVs = 3000.0
    avgCNVSize = 10000.0
    chromosomes = {}
    readSize = 35.0
    depthCov = 1
    
    File.new(model, 'r').each do |line|
      if line=~ /^genome\s+(\S+)/i
        if $1 == 'haploid' or $1 == 'h'
          @states  = [:normal, :del1, :dup1]
          @genome = 'h'
        elsif $1 == 'd' or $1 == "diploid"  # diploid
          @states = [:normal, :del1, :del2, :dup1, :dup2]
        end
      elsif line=~ /^numCNVs\s+(\S+)/i
        numCNVs =   $1.to_f
      elsif line=~ /^avgCNVSize\s+(\S+)/i
        avgCNVSize  = $1.to_f
      elsif line=~ /^depthCov\s+(\S+)/i
        depthCov = $1.to_f
      elsif line=~ /^readSize\s+(\d+)/i
        readSize = $1.to_f
      end 
    end 
        
    @avgCov = depthCov / readSize
    
    @states.each do |s1|
      @initProb[s1] = Math.log(1.0/@states.size)
      
      @states.each do |s2|
        if s1 == s2 
          if s1 == :normal # stay
            @transitionMatrix[s1][s2] = Math.log(1 - (numCNVs / 3000000000.0))
          else 
            @transitionMatrix[s1][s2] = Math.log(1 - (1.0 / avgCNVSize))
          end
        else
          if s1 == :normal
            @transitionMatrix[s1][s2] = Math.log(( numCNVs / 3000000000.0 ) / (@states.size - 1.0) )
          else
            @transitionMatrix[s1][s2] = Math.log(( 1.0/avgCNVSize ) / (@states.size - 1.0))
          end
        end
      end
    end
    
    setEmissionMatrix()
  end

  def EM()  # estimate transition probabilities
    return 
  end

  def baumWelch  # EM algorithm to estimate transition probabilities
    return	  
  end
  
  ## emission = Pr(depth-coverage | copy-number) * Pr(distance | state) 
  # 
  #  add a uniform component unif:
  #   Pr(distance | CNV neutral) = unif + (1-unif) * Norm(distance | CNV-neutral) 
  #  other: 
  
  def setEmissionMatrix
    ## Poisson model
    obsArray = (0..@maxEmission).to_a
    
    lambda = @avgCov
    if @genome == 'd'
      del1 = lambda / 2
      del2 = lambda / 100.0
      dup1 = lambda * 1.5
      dup2 = lambda * 2.0
    else
      del1 = lambda / 100.0
      dup1 = lambda * 2.0
      del2 = del1  # not useful
      dup2 = dup1   # not useful
    end

    
    # log(Poisson(k,lambda)) = klog(lambda) - lambda - SUM(j)_1..k
    obsArray.each do |i|
      @emissionMatrix[:normal][i] = logPoisson(i,lambda)
      @emissionMatrix[:del1][i] = logPoisson(i,del1)
      @emissionMatrix[:del2][i] = logPoisson(i,del2)
      @emissionMatrix[:dup1][i] = logPoisson(i,dup1)
      @emissionMatrix[:dup2][i] = logPoisson(i,dup2)
    end
  end

  def emission(state, cov)
    if cov > @maxEmission 
      cov = @maxEmission
    end
    return @emissionMatrix[state][cov]  
  end
  
  def logPoisson(k,l)
 # log(Poisson(k,lambda)) = klog(lambda) - lambda - SUM(j)_1..k
    p = k * Math.log(l) - l
    1.upto(k) do |i|
      p += -1 * Math.log(i)
    end
    return p
  end
  
  def Viterbi(chr, name)
    chrSize = chr.size
    delta = Hash.new()  # store the prob 
    prevMap = Array.new(chrSize + 1)     # store the previous best state to current
    prevMap.map! {|i| i = Hash.new() }
    
    # initialization
    # pi[state] : inital prob of the state
    
    @states.each do |state|
      delta[state] = @initProb[state] + emission(state, chr[0])
      prevMap[0][state] = state;
    end

    prob = {} 
    # recursion
    1.upto(chr.size-1) do |i|  # position i
      @states.each do |state|
        @states.each do |prevState|
          prob[prevState] = delta[prevState] + @transitionMatrix[prevState][state]
        end
        ## this may be too costly: sort each time
        prevMap[i][state] = prob.keys.sort {|a,b| prob[b] <=> prob[a]}[0] 
        delta[state] = prob[prevMap[i][state]] + emission(state, chr[i])
      end
    end
 
    # termination
    lastBest = delta.keys.sort {|a,b| delta[b] <=> delta[a]}[0]
    prevMap[chr.size][:normal] = lastBest
    flag = 0
    bestPath = []
##     bestPath[chrSize -1] = delta.keys.sort {|a,b| delta[b] <=> delta[a]}[0]

    #  Path (state sequence) backtracking
    
    s = 0
    e = 0 
    cnvs = []
    lastBest = :normal
    (chr.size).downto(1) do |i|
##      prevMap.each do |y|
      state = prevMap[i][lastBest]
     
      s = i
      if state != :normal  
        if flag == 0 # enter a cnv
          flag = 1
          e = i
        elsif state != lastBest  # continue # from a CNV to a new cnv
          cnv = {}
          cnv[:s] = s
          cnv[:e] = e
          cnv[:type] = lastBest
          cnvs << cnv
          e = i
        end  
       
      elsif  flag != 0  # exit a CNV
        cnv = {} 
        cnv[:s] = s
        cnv[:e] = e
        cnv[:type] = lastBest 
        cnvs << cnv
        flag = 0
        s = i
        e = i
      end
      lastBest = state
##       bestPath[i-1] = prevMap[i][bestPath[i]]
    end
    
    cnvs.reverse.each do |cnv|
      puts "#{name}\t#{cnv[:type]}\t#{cnv[:s]}\t#{cnv[:e]}\t#{cnv[:e] - cnv[:s] + 1}"
    end
   
  end  
  
end

def help 
  $stderr.puts "Usage: ruby #{File.basename($0)} -r ref.fa -m mapview -p params  [ -o output_prefix] "
  $stderr.puts " --genome h means haploid genome; default is d, diploid genome"
  exit
end

if __FILE__ == $0 
  require 'getoptlong'

  opts = GetoptLong.new(
        ["--ref", "-r", GetoptLong::REQUIRED_ARGUMENT],
        ["--mapping", "-m",GetoptLong::REQUIRED_ARGUMENT],
        ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
        ["--params", "-p", GetoptLong::REQUIRED_ARGUMENT],
        ["--help", "-h", GetoptLong::NO_ARGUMENT]
  )
  optHash = {}
  opts.each do |opt, arg|
    optHash[opt] = arg
  end
  if optHash.key?("--help")
    help()
    exit
  end  
  if optHash.key?("--output") 
    outputPrefix = optHash["--output"]
  else 
    outputPrefix = "cnvHMM"
  end

  
  mapping = Mapping.new(optHash["--mapping"], optHash["--ref"])
  $stderr.puts "Mapping results loaded."
  
  hmm  = CnvHMM.new(optHash["--params"])
  
  mapping.chromosomes.each do |name, chr|
#   hmm.EM()
    hmm.Viterbi(chr.coverage, name)
    
  end

end




	

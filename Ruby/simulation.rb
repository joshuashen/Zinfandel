#!/usr/bin/env ruby -w

### need to make sure the cancer type CNV is a dup in the first place

class Cnv
    attr_accessor :s, :e, :chr, :indeltype, :copyNumber
    
    def initialize()
        @copyNumber = -1  # unknown copy number
        @indeltype = "indel"
    end
end

class CNVSimulation
  def initialize(ref, genome)
    @genome = genome
    @seq = {}
    @chrs = {}
    @cnvs = []
    readRef(ref)
  end

  def readRef(ref)
    chr = '';
    File.new(ref,'r').each do |line|
      if line=~/^>(\S+)/  
        chr = $1
        $stderr.puts "Loading #{chr}..."
        @chrs[chr] =[]
        @seq[chr] = ''
      else
        @seq[chr] << line.chomp
      end
    end
    $stderr.puts "#{ref} loaded."    
  end

  def simulate(prefix)  # produce the sequence of target genome in fasta format
    faOut = File.new(prefix + "_target.fa",'w')
    cnvOut = File.new(prefix + "_true_cnv.txt", 'w')
    
 
    @cnvs.each do |cnv|
      chr = cnv.chr
      @chrs[chr] << cnv
    end
 
    if @genome == 'diploid'  # simulate diploid target genome
      @chrs.keys.each do |chr|
        sec = chr + '_d'
        ncnv = @chrs[chr].size / 2
        @chrs[chr].shuffle!
        @chrs[sec] = @chrs[chr].slice!(0,ncnv)  # randomly pull half CNVs to the second chromosome   
      end
    end    
    
    @chrs.keys.sort.each do |chr|   
      target = ''
      if (chr=~ /^(\S+)\_d/)  # the second chromosome of a pair
        original = $1
        target = @seq[original] 
      else
        target = @seq[chr]
      end   
      shifts = 0
      @chrs[chr].sort {|a,b| a.s <=> b.s}.each do |cnv|

        if cnv.indeltype == 'del' 
          target[(shifts + cnv.s - 1) .. (shifts + cnv.e - 1)] = 'X' * (cnv.e - cnv.s + 1) # delete
          cnvOut.puts "#{cnv.chr}\tDEL\t#{cnv.s}\t#{cnv.e}\t#{cnv.e - cnv.s + 1}\t1"
        elsif cnv.indeltype == 'dup' or cnv.indeltype == 'cancer'

          copyNumber = cnv.copyNumber
          if copyNumber < 0  # un-defined, should be a single duplication
            copyNumber = 2
          end
          cnvOut.puts "#{cnv.chr}\tDUP\t#{cnv.s}\t#{cnv.e}\t#{cnv.e - cnv.s + 1}\t#{copyNumber + 1}"          
          inserted = target[(shifts + cnv.s - 1) .. (shifts + cnv.e - 1)] * (copyNumber - 1 )
          target[(shifts + cnv.e )..(shifts + cnv.e - 1)] = inserted
          shifts += inserted.size
        end
      end
      
      target.tr('X', '')
      faOut.puts ">#{chr}"
## max 80 bp each line
      chrSize = target.size 
      0.upto(chrSize/80) do |i|
        faOut.print(target[i*80..(i+1)*80-1] + "\n")
      end

      target = ''
    end
    
    faOut.close
    cnvOut.close
  end

  def highCopyNumber(copyNumberArray)  # randomly pick #num duplications, and make the copy number to #copyNumber
    $stderr.puts "number of cnvs: #{@cnvs.size}"
    j = 0
    copyNumberArray.each do |i|
      @cnvs[j].indeltype = "cancer"
      @cnvs[j].copyNumber = i
      j += 1
    end
  end

  def aneuploid(chr, s, e, indeltype)  # large-scale duplication or deletion: half or the entire chromosome is duplicated or deleted
    cnv = Cnv.new()
    cnv.chr, cnv.s, cnv.e, cnv.indeltype = chr, s, e, indeltype
    @cnvs << cnv
  end

  def drawfromCommon(list, num, minSize)
    cnvArray = readCommonCNVs(list, minSize)
    $stderr.puts "number of common cnvs: #{cnvArray.size}"
    @cnvs.concat(cnvArray.shuffle[0..num-1])
  end

  def readCommonCNVs(list, minSize)
    cnvArray = []
    ## format: 
#  VariationID     Landmark        Chr     Start   End     VariationType   LocusChr        LocusStart      LocusEnd        Reference       PubMedID        Method/platform Gain    Loss    TotalGainLossInv        SampleSize
# Variation_10631 chr1:102600719..102601411       chr1    102600719       102601411       InDel   chr1    102600719       102601411       Conrad et al. (2005)  16327808        Mendelianinconsistencies     
    
    File.new(list,'r').each do |line|
      cols=line.split(/\s+/)
      cnv = Cnv.new()
      cnv.chr, cnv.s, cnv.e, cnv.indeltype = cols[2], cols[3].to_i, cols[4].to_i, cols[5]
      
      ss = cnv.e - cnv.s + 1
      
      next unless ss >= minSize
      
      # make sure the chr is in the reference; 
      if  @chrs.key?(cnv.chr) and (cnv.indeltype =~ /InDel/i or cnv.indeltype =~ /CopyNumber/i ) # from common CNV file
        
        cnv.indeltype = ["dup", "del"].shuffle[0]  # randomly choose dup or del 
        cnvArray << cnv  
      end
    end
    return cnvArray
  end
end

# implement Array#shuffle ; necessary for Ruby 1.8.6 or earlier
# class Array
#  def shuffle
#    self.sort_by { rand }
#  end
# end

def help
  $stderr.puts "Usage: ruby #{File.basename($0)} -r ref.fa -d cnv_db -n num_of_CNVs [ -c cancer_amplification_number] [-m min_CNV_size] [--genome haploid/diploid] [-o prefix_of_output]"
end

# main()
if __FILE__ == $0 
  require 'getoptlong'

  opts = GetoptLong.new(
        ["--ref", "-r", GetoptLong::REQUIRED_ARGUMENT],
        ["--common", "-d", GetoptLong::REQUIRED_ARGUMENT],
        ["--number", "-n",GetoptLong::REQUIRED_ARGUMENT],
        ["--output", "-o", GetoptLong::OPTIONAL_ARGUMENT],
        ["--cancer", "-c", GetoptLong::OPTIONAL_ARGUMENT],
        ["--genome", "-g", GetoptLong::OPTIONAL_ARGUMENT],
        ["--minSize", "-m", GetoptLong::OPTIONAL_ARGUMENT],
        ["--aneuploid", "-a", GetoptLong::OPTIONAL_ARGUMENT],
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
    outputPrefix = "simulation"
  end
  
  minSize = 500 
  
  if optHash.key?("--minSize") 
    minSize = optHash["--minSize"].to_i
  end
  
  genome = "haploid"
  if optHash["--genome"] == "diploid" or optHash["--genome"] == 'd'
    genome = "diploid"
  end
  
  simu = CNVSimulation.new(optHash["--ref"], genome)
  simu.drawfromCommon(optHash["--common"], optHash["--number"].to_i, minSize)
  if optHash.key?("--cancer")
    cnumArray = optHash["--cancer"].split(',').map {|i| i.to_i}
    simu.highCopyNumber(cnumArray)  
  end
  simu.simulate(outputPrefix)
end
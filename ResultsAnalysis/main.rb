class Cnv
  attr_accessor  :type, :size, :startPos, :endPos
  def initialize(type, size, startPos, endPos)
    @type = type
    @size = size
    @startPos = startPos
    @endPos = endPos
  end
end    
    
    known = Array.new
    results = Array.new
    
    numDeletions = 50
    errorRange = 200
    
    matches = 0
    falsePositives = 0
        
    File.new(ARGV[0], 'r').each do |line|
      cols = line.split(" ");
      name, type, startPos, endPos, size = cols[0], cols[1], cols[2].to_i, cols[3].to_i, cols[4].to_i 
      
      if type.eql?("DEL")
        known << Cnv.new(type, size, startPos, endPos)        
      end
    end        
    
    File.new(ARGV[1], 'r').each do |line|
      if /^EcoliGenome/ =~ line
        cols = line.split(" ");
        name, type, startPos, endPos, size = cols[0], cols[1], cols[2].to_i, cols[3].to_i, cols[4].to_i 
        results << Cnv.new(type, size, startPos, endPos)        
      end      
    end
    
 
    #puts results.size
    #puts known.size
    
    known.each do |i|
      valid = false
      results.each do |j|
        if j.startPos >= i.startPos - errorRange and j.endPos <= i.endPos + errorRange
          valid = true;
        end                
      end
      if valid == true
        matches += 1
      end
    end
    
    results.each do |i|
      found = false
      known.each do |j|
        if i.startPos >= j.startPos - errorRange and i.endPos <= j.endPos + errorRange
          found = true
        end
      end
      if found ==  false
        falsePositives += 1
      end
    end       

    deletions = known.size
    numResults = results.size
    
    successRate = matches.to_f / deletions.to_f
    falsePositiveRate = falsePositives.to_f / numResults.to_f
    
    puts successRate
    puts falsePositiveRate
    

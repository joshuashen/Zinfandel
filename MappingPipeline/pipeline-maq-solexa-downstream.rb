## 
prefix = ARGV[0]
ref = ARGV[1]
maps = ARGV[2..-1].join("\t") 
# puts maps

# merge all .map files
system("maq mapmerge #{prefix}.map #{maps}")

# assemble
system("maq assemble #{prefix}.cns #{ref} #{prefix}.map 2> #{prefix}.cns.log")

# call SNP
system("maq cns2snp #{prefix}.cns > #{prefix}.cns.SNP")

# pileup
system("maq pileup -sP #{ref} #{prefix}.map > #{prefix}.pileup")



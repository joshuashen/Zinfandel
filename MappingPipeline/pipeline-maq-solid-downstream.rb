## 
prefix = ARGV[0]
ref = ARGV[1]
maps = ARGV[2..-1].join("\t") 
# puts maps

# merge all .map files
system("maq mapmerge #{prefix}.cs.map #{maps}")

#convert map from color space to nt space
system("maq csmap2nt #{prefix}.nt.map #{ref}  #{prefix}.cs.map")

# assemble
system("maq assemble #{prefix}.cns #{ref} #{prefix}.nt.map 2> #{prefix}.cns.log")

# call SNP
system("maq cns2snp #{prefix}.cns > #{prefix}.cns.SNP")

# pileup
system("maq pileup -sP -m 3 #{ref} #{prefix}.nt.map > #{prefix}.pileup")

# mapview
system("maq mapview -b #{prefix}.nt.map > #{prefix}.nt.map.mapview")



# 

prefix = ARGV[0]

ref = ARGV[1] # ref in csbfa
  
# fastq to bfq

system("/nfs/apollo/2/c2b2/users/saec/Software/maq/bin/maq fastq2bfq #{prefix}.read1.fastq #{prefix}.read1.bfq")
system("/nfs/apollo/2/c2b2/users/saec/Software/maq/bin/maq fastq2bfq #{prefix}.read2.fastq #{prefix}.read2.bfq")

system("/nfs/apollo/2/c2b2/users/saec/Software/maq/bin/maq map -a 5000 -c #{prefix}.map #{ref}  #{prefix}.read1.bfq #{prefix}.read2.bfq 2>> #{prefix}.maq.run")


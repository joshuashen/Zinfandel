#!/bin/sh
#$ -cwd

/usr/java/jdk1.5.0_07/bin/java -Xms4096m -Xmx8192m -classpath /ifs/scratch/c2b2/ip_lab/yg2181/Java-CNV/ Main -r /ifs/scratch/c2b2/ip_lab/yshen/Atlas-CNV/human/FixedSizeCNVs/chr18.fa -m /ifs/scratch/c2b2/ip_lab/yg2181/MapviewFiles/100-0.2Random_210K_pairs.mapview -p /ifs/scratch/c2b2/ip_lab/yg2181/MapviewFiles/params/100-0.2params.random > 100-0.2.txt

exit
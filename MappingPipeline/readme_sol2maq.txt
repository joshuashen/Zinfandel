Sol2Maq.pl
Author: Ye Liu
Last Updated: 8/28/08

How it works:
	The purpose of this script is to take Solexa export files and run them through maq. It accepts a
directory as input and runs all the paired (and only paired) Solexa reads in it, then outputting 
the intermediate files and results in a given output directory. All intermediate outputs are saved
in each lane's own subdirectory; a Run.txt file is saved in the user-specified general output directory, 
which includes all the redirected stderr outputs during the execution. The expected run time for 
7 lanes, 21.5G of Solexa reads mapped against the entire genome of C.Elegans is about 11 hours.

Since maq's mapping tools work at an optimum rate when read files are between 1 to 3 million reads
long, the script cuts each file pair into 2 million line reads and maps each smaller pair individually
before merging the maps at the end. At the same time, if the user opted to filter the reads by
chromosomes and/or mapped positions from the Solexa files, the script will fill each smaller read
file with 2 million of the filtered results (so with or without filtering, the smaller read files
will always contain 2 million or less lines). Smaller files and output files will retain the original
files' prefixes.

The smaller files are processed in pairs. First, they are converted from Solexa's format into .fastq.
The default highest quality score allowed is 40, but the user can change this as a command line option.
Each converted pair will be mapped against a reference specified by the user and the map will be
saved in the lane's subdirectory. When all the smaller files of a lane are processed, the script calls 
on maq to combine the maps into one for the lane. And when all the lanes are processed, maq will 
combine the maps from every lane and call an consensus. Outputs such as snp calls and indels will be
pulled from the consensus file and saved in the general output directory.


Requirements:
	The script expects paired Solexa reads (2 per lane), all directory paths listed with the last slash 
(/), reference files in either .bfa or .fasta, and can accept some additional options including filtering
and maximum quality scores. If a .fasta file is given, the script will check for a .bfa file with the same
name and if it exists, the .bfa file will be used. Otherwise the .fasta file will be converted and saved
at the same directory; this checking will save time needed to convert reference sequences. There is no way
to run only selected lanes in a directory; all pairs will be processed. Separate the files into different
directories if only some lanes are to be run. Pairs will be recognized as having the same prefix, namely
s_n_1 and s_n_2; other naming conventions will not work correctly. If filtering by position is used, both
a start and an end position must be specified.


Outputs:
- in each lane's output directory:
  - 2-million-read files: [name]_[1 or 2]_export.n.txt, [name]_[1 or 2]_export.n.fastq, [name]_[1 or 2]_
    export.n.bfq, where n is an integer
  - maps: part[n].map, [name].all.map
- in the user-specified general output directory:
  - full map: all.map
  - consensus: all.cns
  - snp: all.cns.snp
  - filtered snp: all.cns.filtered.snp
  - fq: all.cns.fq
  - indelpe: all.indelpe
  - run file from stderr: Run.txt


Command-line inputs (required, separated by a space): 
1. Solexa output directory
	The full path of the directory containing Solexa output files; path must end in "/".
2. Maq output directory
	The full path of the general output directory where Maq outputs will be saved. Path must end in 
  "/". Each Solexa lane will create its own subdirectory under this.
3. Reference file
	The full path and name of the reference file. Can be either .bfa, .fasta, or .txt. If *.fastq or 
  *.txt is given, the script will perform a check at its directory for *.fasta.bfa or *.txt.bfa and use 
  this file instead if it exists; this assumes that the bfa file corresponds to the fasta file with the 
  same name! Otherwise the script will convert the *.fasta(txt) file to a bfa file called *.fasta(txt).bfa.


Command-line options (optional, separated by a space):
-h (--help)
	Prints the command-line format needed to run the script.
-c (--chr)
	Filtering by chromsome. It assumes that the mapped chromosome in the Solexa files will be in the format
  "chrN.fa", so accepts a capital letter Roman numeral. For example, to get reads mapped only onto chromosome
  6, type in "-c VI" or "--chr VI". Can only be used to filter out one single chromosome.
-s (--start), -e (--end)
	Filtering by mapped position. Both a start and an end must be chosen, with the start position being
  a smaller integer than the end position. 
-q (--maxQual)
	Maximum quality accepted in the Solexa files. The default max is 40. Quality scores greater
  than the maximum number will not be converted in the file conversion step and the read containing it will not 
  be accepted as valid. 

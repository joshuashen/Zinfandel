Readme: Simple Hidden Markov Model

The program requires 3 input files:
1. Fasta file containing the chromosomes sequence of the reference 
2. Mapview file containing depth coverage information
3. Parameter file with genome type, number of CNVs, CNV size, depth coverage, and readsize 

Running the program
After compiling all the .java files, use command in the following format:
java Main.java -r reference.fasta -m mapview -p parameters

The fasta file follows the -r argument, mapview file follows -m argument, and parameters
file follows the -p argument.

Note, the arguments can be entered in any order so long as the correct file follows the 
corresponding argument. 

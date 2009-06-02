Solid.sh (Solid_LargeFiles.sh)
Author: Ye Liu
Last Updated: 8/28/08

How it works:
	The purpose of this script is to take as input a pair of ABI/SOLiD reads, along with their 
quality files, and run them through maq. Each of the SOLiD reads is broken down into 2-million reads, 
and these subfiles are processed and mapped together. In Solid.sh, all the submaps are combined at the
end and consensus files are pulled; in Solid_LargeFiles.sh, n submaps are combined together before these
combined maps are combined yet again at the end, therefore reducing the total number of files that
has to be opened simultaneously at any time. This is at the cost of slightly longer anticipated run 
time. For 30.8G of reads/qualities, the run time for Solid.sh is expected to be about 22.65 hours.

A pair of SOLiD reads, along with their quality files, are first combined using maq's solid2fastq.pl
script. The output, two compressed .gz files, are decompressed for the .fastq files. Then the script
converts the reference files as needed; if the reference files in the correct format already exist,
they will not be converted again.

The reference file names are assumed to be .[extension] added on to the 
user input's name: for example, the .csbfa file is assumed to have the name *.fasta.csbfa.
Converted reference files are saved in the same directory as the given reference files.

Following this, the pair of .fastq files are broken down into 2-million-read subfiles to fit maq's 
optimum running range, and each subfile pair is converted to binary format and mapped. 

In Solid.sh, these part maps' names are added to a list that is used to combine all maps at the end.
In Solid_LargeFiles.sh, a variable nMaps is set to the square root of the number of pair files (integer only),
and every nMaps number of part maps is combined together. At the end, all the partly combined maps 
(again, about nMaps number of them) are combined. This is done so that as SOLiD's file size increases
and the number of part maps increases, the machine will not run out of memory trying to open all the
maps simultaneously to combine them. The square root of the number of pair files is used to minimize 
memory needed by ensuring that both rounds of map-combining will open about the same number of files.

After all part maps have been combined, the script calls on maq to create a consensus file and pull out
files such as snp calls and indels. A run.txt file is also created, containing all the stderr outputs
during the run of the script.


Requirements:
	The script must be run in the same directory as the original reads and quality files!! However,
the outputs can be redirected to other directories. All directories and files, but not prefixes, in the 
command-line arguments should contain the full path name; all directories should end in a slash (/).


Note:
	The numbers in the[name].[number].map partial maps after the first round of combination for 
Solid_LargeFiles.sh will not be consecutive; do not be alarmed if some numbers seem missing. This is due to
the implementation of the codes. 

Outputs:
- 2-million-read subfiles in the format of [name].read[1 or 2].[n].fastq
- corresponding .bfq files for each .fastq subfile
- partial maps in the format of [name].part[n].map
- Solid_LargeFiles.sh only: partial maps in the format of [name].[n].map, from the second round of combining
- full map: [name].all.map
- consensus file: [name].cns.cns
- snps: [name].cns.snp
- filtered snp: [name].cns.filtered.snp
- fq: [name].cns.fq
- indelpe: [name].indelpe
- run file from stderr: [name].run.txt
- whole read files: [name].read[1 or 2].fastq 
- reference files, as needed

All outputs are in the user-specified output directory, except for the last three in the above list. The reference
files will be saved where the given reference was, and the other two, in the directory containing the original 
reads and quality files.

Command-line inputs (required, separated by a space):
1. Prefix name
	Prefix of the SOLiD/quality files; it is everything before "F3" and "R3", including the underscore at 
  the end. Does not need to contain the full path name.
2. Output name
	Prefix of the output files; can be different from input prefix. Does not need to contain full path.
3. Reference file
	Must contain the full path name, and should be in .fasta or .txt format.
4. Output directory
	Should contain the full path name and end in a slash (/).

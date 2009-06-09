package cnv_simulation;
import cnv_HMM.HMM;
import java.io.*;
import java.util.*;

public class GenomeSegments {
    public static enum chrID {
	chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10,
	chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20,
	chr21, chr22, chrX, chrY, chrM
		};

    public static int[] chr_length = {
	247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
	158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
	114142980, 106368585, 100338915, 88827254, 78774742, 76120000,
	63811651, 62440000, 46944323, 49691432, 154913754, 57772954, 16571
   };


    static class SNP implements Comparable<SNP> {
	chrID chromosome;
	int location;
	int relative_location = -1;
	boolean heterozygous;
	char new_base;

	public int compareTo(SNP o) {
	    int chr_compare = chromosome.compareTo(o.chromosome);
	    return chr_compare == 0 ? location - o.location : chr_compare;
	}

	SNP(chrID c, int l) { chromosome = c; location = l; }

	SNP(String line) {
	    String[] tokens = line.split("\\t");
	    chromosome = getchrID(tokens[0]);
	    location = Integer.parseInt(tokens[4]);
	    heterozygous = tokens[2].startsWith("hetero");
	    new_base = tokens[7].charAt(2);
	}

	void apply(char[] seq) {
	    if(seq != null)
		seq[relative_location] = new_base;
	}
    }

    public static class Segment implements Comparable<Segment> {
	File file;
	chrID chromosome;
	int start = -1;
	int length;
	boolean has_hetero = false;
	int num_reads = 0;
	NavigableSet<CNV> variations = null;

	Segment() {}

	public Segment(chrID c, int s) { chromosome = c; start = s; }

	boolean addVariation(CNV variation) {
	    if(variation.heterozygous)
		has_hetero = true;
	    if(variations == null)
		variations = new TreeSet<CNV>(CNV.ByStart);
	    CNV previous = variations.ceiling(variation);
	    int min_start = previous == null ? 0 : previous.relative_start+previous.size;
	    CNV next = variations.floor(variation);
	    int max_start = (next == null ? length : next.relative_start) - variation.size;
	    if(max_start < min_start)
		return false;
	    variation.relative_start = Math.max(min_start, Math.min(max_start, variation.relative_start));
	    if(start >= 0)
		variation.absolute_start = variation.relative_start + start;
	    variations.add(variation);
	    return true;
	}

	public int compareTo(chrID chromosome, int loc)
	{
	    if(chromosome != this.chromosome)
		return this.chromosome.compareTo(chromosome);
	    if(loc < start)
		return 1;
	    if(loc >= start+length)
		return -1;
	    return 0;
	}

	public void set_carpet_read() {
	    num_reads = length - sequence_read.mean_block_size;
	}

	public int compareTo(Segment o) {
	    int chr_compare = chromosome.compareTo(o.chromosome);
	    int result = chr_compare == 0 ? start - o.start : chr_compare;
	    return result;
	}

    Collection<CNV> run_HMM(PrintStream out)
    throws Exception
    {
	long distance_between_duplications = Simulation.var_sizes[0] == null
		? Integer.MAX_VALUE :
		(Simulation.numbered.total_length + Simulation.XY.total_length) / Simulation.var_sizes[0].length;
	long distance_between_deletions = Simulation.var_sizes[1] == null
		? Integer.MAX_VALUE :
		(Simulation.numbered.total_length + Simulation.XY.total_length) / Simulation.var_sizes[1].length;
	HMM hmm = null; // new HMM(chromosome, distance_between_duplications, distance_between_deletions, Simulation.mean_var_size, Simulation.homozygous_probability, 1, null, null, null);
	CNV current_cnv = null;
	for(int pos=0; pos <= length; pos++){
	    int state = pos == length ? hmm.normal_state : hmm.make_transition();
	    if(state != hmm.normal_state && current_cnv == null){
		current_cnv = new CNV(
			state > hmm.normal_state ? CNV.CNV_type.DUP : CNV.CNV_type.DEL,
			chromosome, length, 1, state==1 || state==3);
		current_cnv.relative_start = pos;
	    } else if(state == hmm.normal_state && current_cnv != null){
		current_cnv.size = pos - current_cnv.relative_start;
		addVariation(current_cnv);
		current_cnv = null;
	    }
	}
	return variations == null ? new Vector<CNV>() : variations;
    }

    void do_simulation()
    throws IOException
    {
	if(num_reads > 0)
	{
	    char[] reference_sequence = new char[length];
	    BufferedReader infile = new BufferedReader(new FileReader(file));
	    start = Integer.parseInt(infile.readLine().split("\\s+")[1]);
	    int index = 0;
	    for(String next = infile.readLine(); next != null; index += next.length(), next = infile.readLine())
		next.toUpperCase().getChars(0, next.length(), reference_sequence, index);
	    String altered_sequence = variations == null ? null : new String(reference_sequence);
	    String altered_sequence2 = has_hetero ? altered_sequence : null;
	    if(Simulation.SNPs != null){
		List<SNP> SNPs = Simulation.SNPs.get(chromosome.ordinal());
		int start_index = Collections.binarySearch(SNPs, new SNP(chromosome, start));
		if(start_index < 0)
		    start_index = - start_index - 1;
		int end_index = Collections.binarySearch(SNPs, new SNP(chromosome, start+length));
		if(end_index < 0)
		    end_index = - end_index - 1;
		if(end_index > start_index){
		  char[] seq1 = reference_sequence, seq2 = null;
		  if(altered_sequence != null){
		    seq1 = new char[altered_sequence.length()];
		    altered_sequence.getChars(0, seq1.length, seq1, 0);
		  }
		  if(has_hetero){
		    seq2 = new char[altered_sequence2.length()];
		    altered_sequence2.getChars(0, seq2.length, seq2, 0);
		  }
		  for(int s = start_index; s < end_index; s++){
		    SNP snp = SNPs.get(s);
		    if(snp.heterozygous && !has_hetero){
			has_hetero = true;
			seq2 = seq1.clone();
		    }
		    snp.relative_location = snp.location - start;
		    if(snp.heterozygous)
			if(Simulation.generator.nextFloat() < .5)
			    snp.apply(seq1);
			else
			    snp.apply(seq2);
		    else {
			snp.apply(seq1);
			snp.apply(seq2);
		    }
		  }
		  if(seq1 != reference_sequence)
		    altered_sequence = new String(seq1);
		  if(seq2 != null)
		    altered_sequence2 = new String(seq2);
		}
	    }
	    if(variations != null)
	    for(Iterator<CNV> iter = variations.iterator(); iter.hasNext(); ){
		CNV variation = iter.next();
		variation.absolute_start = start + variation.relative_start;
		if(variation.heterozygous)
		    if(Simulation.generator.nextFloat() < .5)
			altered_sequence = variation.apply(altered_sequence);
		    else
			altered_sequence2 = variation.apply(altered_sequence2);
		else {
		    altered_sequence = variation.apply(altered_sequence);
		    altered_sequence2 = variation.apply(altered_sequence2);
		}
	    }
	    long file_num = Simulation.reads_per_file > 0
			? Simulation.num_reads_done / Simulation.reads_per_file + 1 : -1L;
	    String file_name = file_num > 0 ? Simulation.output_name + '.' + file_num : Simulation.output_name;
	    PrintStream reads_file_F = new PrintStream(new FileOutputStream(file_name + ".F", Simulation.num_reads_done > 0));
	    PrintStream reads_file_R = new PrintStream(new FileOutputStream(file_name + ".R", Simulation.num_reads_done > 0));
	    for(int i = 0; i < num_reads; i++){
		if(file_num > 0 && Simulation.num_reads_done++ / Simulation.reads_per_file >= file_num){
		    reads_file_F.close();
		    reads_file_R.close();
		    file_name = Simulation.output_name + '.' + (++file_num);
		    reads_file_F = new PrintStream(new FileOutputStream(file_name + ".F", false));
		    reads_file_R = new PrintStream(new FileOutputStream(file_name + ".R", false));
		}
		sequence_read read = sequence_read.perform_read
		  (reference_sequence,
			has_hetero && Simulation.generator.nextFloat() < .5
			? altered_sequence2 : altered_sequence,
			i);
		String read_name =
			(Simulation.fastq_format ? "@" : ">")
			+ "read" + chromosome
			+ '_' + (start+read.start) + '_'
			+ (start+read.start+read.block_size-read.read_length);
		reads_file_F.println(read_name + "/1");
		reads_file_F.println(read.read[0]);
		reads_file_R.println(read_name + "/2");
		reads_file_R.println(read.read[1]);
		if(Simulation.fastq_format){
		    reads_file_F.println("+");
		    reads_file_R.println("+");
		    reads_file_F.println(Simulation.quality_string);
		    reads_file_R.println(Simulation.quality_string);
		}
	    }
	    reads_file_F.close();
	    reads_file_R.close();
	}
    }

    }

    public Segment[] segments;
    long total_length;
    Map<chrID, Integer> length_per_chromosome;
    static final long genome_length = 2974076126L;
    boolean sorted = false;

    private static class XY_segment implements FilenameFilter {
	public boolean accept(File dir, String name) {
	  try {
	    return (name.startsWith("chrX.segment") || name.startsWith("chrY.segment"))
			&& Integer.parseInt(name.substring(12)) > 0;
	  } catch (NumberFormatException e) { return false; }
	}
    }

    private static class numbered_segment implements FilenameFilter {
	public boolean accept(File dir, String name) {
	  try {
	    int segindex = name.indexOf(".segment");
	    return segindex > 3 && name.startsWith("chr") &&
		Integer.parseInt(name.substring(3,segindex)) > 0 &&
		Integer.parseInt(name.substring(segindex+8)) > 0;
	  } catch (NumberFormatException e) { return false; }
	}
    }

    public GenomeSegments(String dirname, boolean XY){
	File dir = new File(dirname);
	String[] file_names = dir.list(XY ? new XY_segment() : new numbered_segment());
	List<Segment> seglist = new Vector<Segment>();
	total_length = 0;
	length_per_chromosome = new HashMap<chrID, Integer>();
	for(int i = 0; i < file_names.length; i++){
	    Segment seg = new Segment();
	    seg.chromosome = getchrID(file_names[i].substring(3, file_names[i].indexOf('.')));
	    if(Simulation.chromosomes != null && !Simulation.chromosomes.contains(seg.chromosome))
		continue;
	    seg.file = new File(dir, file_names[i]);
	    int num_lines = (int)(seg.file.length() - 21) / 51 + 1;
	    total_length += seg.length = (int)seg.file.length() - 20 - num_lines;
	    Integer previous_chromosome_length = length_per_chromosome.get(seg.chromosome);
	    if(previous_chromosome_length == null)
		previous_chromosome_length = 0;
	    length_per_chromosome.put(seg.chromosome,
				previous_chromosome_length.intValue() +
				seg.length);
	    seglist.add(seg);
	}
	segments = new Segment[seglist.size()];
	seglist.toArray(segments);
    }

    public static chrID getchrID(String chrname)
    {
	return chrID.valueOf(chrID.class, "chr" + chrname);
    }

    public static boolean isNumbered(GenomeSegments.chrID chr) {
	return chr.compareTo(GenomeSegments.chrID.chrX) < 0;
    }

    public void sort()
    throws IOException
    {
	if(sorted)
	    return;
	for(int seg = 0; seg < segments.length; seg++){
	    BufferedReader infile = new BufferedReader(new FileReader(segments[seg].file));
	    segments[seg].start = Integer.parseInt(infile.readLine().split("\\s+")[1]);
	    infile.close();
	}
	Arrays.sort(segments);
	sorted = true;
    }

    Collection<CNV> run_HMMs(PrintStream out)
    throws Exception
    {
	sort();
	Collection<CNV> result = new Vector<CNV>();
	if(segments.length > 0){
	chrID current_chr = segments[0].chromosome;
	int current_pos = 1;
	for(int segment = 0; segment < segments.length; segment++){
	  if(segments[segment].chromosome != current_chr){
	    if(out != null && chr_length[current_chr.ordinal()] >= current_pos)
		out.println(current_chr.toString() + ':' + current_pos
			+ '-' + chr_length[current_chr.ordinal()] + " :");
	    current_chr = segments[segment].chromosome;
	    current_pos = 1;
	  }
	  if(out != null && segments[segment].start > current_pos)
		out.println(current_chr.toString() + ':' + current_pos
			+ '-' + (segments[segment].start-1) + " :");
	  result.addAll(segments[segment].run_HMM(out));
	  current_pos = segments[segment].start + segments[segment].length;
	}
	if(out != null && chr_length[current_chr.ordinal()] >= current_pos)
		out.println(current_chr.toString() + ':' + current_pos
			+ '-' + chr_length[current_chr.ordinal()] + " :");
	}
	return result;
    }
}

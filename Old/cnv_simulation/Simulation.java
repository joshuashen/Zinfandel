package cnv_simulation;
import java.io.*;
import java.util.*;

public class Simulation {
    static GenomeSegments numbered;
    static GenomeSegments XY;
    public static Random generator = new Random();
    static long num_reads_done = 0;
    static int reads_per_file = -1;
    static String output_name = "read.results";
    public static Collection<GenomeSegments.chrID> chromosomes = null;
    static boolean fastq_format = false;
    static String quality_string = "";
    static List< List<GenomeSegments.SNP> > SNPs = null;
    static boolean use_HMM = false;
    static boolean HMM_dummy = false;
    static float homozygous_probability = (float).1;
    static int mean_var_size = 80000;
    static double coverage = 0;
    static int[][] var_sizes = {null, null, null, null};

    static GenomeSegments.Segment random_segment(boolean include_XY, int min_length) {
	long total_length = numbered.total_length;
	if(include_XY)
	    total_length = total_length * 2 + XY.total_length;
	for(;;){
	    long location = Math.abs(generator.nextLong()) % total_length;
	    GenomeSegments segments = location / numbered.total_length > 1 ? XY : numbered;
	    int index = -1;
	    for(long i = location % numbered.total_length; i >= 0;
		i -= segments.segments[++index].length);
	    if(segments.segments[index].length >= min_length)
		return segments.segments[index];
	}
    }

    static GenomeSegments.Segment segmentAt(int index) {
	return index < numbered.segments.length
		? numbered.segments[index]
		: XY.segments[index - numbered.segments.length];
    }

    static int[] get_var_sizes(String[] tokens) {
	int[] result;
	if(tokens[1].equalsIgnoreCase("sizes")){
	  result = new int[tokens.length-2];
	  for(int i=2; i < tokens.length; i++)
	    result[i-2] = Integer.parseInt(tokens[i]);
	} else {
	  result = new int[Integer.parseInt(tokens[1])];
	  Arrays.fill(result, 0);
	}
	return result;
    }

    public static void main(String[] args) {
      try {
	if(args.length != 2){
	    System.err.println("usage: Simulation param-file sequence-directory");
	    System.exit(0);
	}

	int min_var_size = 100, max_var_size = 1000000, min_big_vars = 0;
	int quality = 20;
	int num_common_CNVs = 0;
	List<CNV> common_CNVs = null;
	List<String> common_CNV_files = null;
	String SNPs_file = null;

	BufferedReader paramfile = new BufferedReader(new FileReader(args[0]));
	for(String paramline = paramfile.readLine(); paramline != null; paramline = paramfile.readLine())
	{
	    String[] tokens = paramline.split("\\s+");
	    if(tokens[0].equalsIgnoreCase("amplifications") || tokens[0].equalsIgnoreCase("duplications"))
		var_sizes[CNV.CNV_type.DUP.ordinal()] = get_var_sizes(tokens);
	    else if(tokens[0].equalsIgnoreCase("deletions"))
		var_sizes[CNV.CNV_type.DEL.ordinal()] = get_var_sizes(tokens);
	    else if(tokens[0].equalsIgnoreCase("inversions"))
		var_sizes[CNV.CNV_type.INV.ordinal()] = get_var_sizes(tokens);
	    else if(tokens[0].equalsIgnoreCase("commonVariations")){
		num_common_CNVs = Integer.parseInt(tokens[1]);
		common_CNV_files = Arrays.asList(tokens).subList(2,tokens.length);
	    } else if(tokens[0].equalsIgnoreCase("homozygouspercentage"))
		homozygous_probability = Float.parseFloat(tokens[1]) / 100;
	    else if(tokens[0].equalsIgnoreCase("variationsizerange")){
		min_var_size = Integer.parseInt(tokens[1]);
		max_var_size = Integer.parseInt(tokens[2]);
	    } else if(tokens[0].equalsIgnoreCase("meanvariationsize"))
		mean_var_size = Integer.parseInt(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("minimumlargevariations"))
		min_big_vars = Integer.parseInt(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("coverage")){
		if(tokens[1].equalsIgnoreCase("carpet"))
		  sequence_read.carpet_read = true;
		else
		  coverage = Double.parseDouble(tokens[1]);
	    } else if(tokens[0].equalsIgnoreCase("readlength"))
		sequence_read.read_length = Integer.parseInt(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("readblockmean") ||
			tokens[0].equalsIgnoreCase("readblocksize"))
		sequence_read.mean_block_size = Integer.parseInt(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("readblockstd"))
		sequence_read.std_block_size = Float.parseFloat(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("readblockdistribution"))
		sequence_read.readSizeDistribution(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("SNPs"))
		SNPs_file = tokens[1];
	    else if(tokens[0].equalsIgnoreCase("readsperfile"))
		reads_per_file = Integer.parseInt(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("outputname"))
		output_name = tokens[1] + ".read.results";
	    else if(tokens[0].equalsIgnoreCase("quality"))
		quality = Integer.parseInt(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("chromosome")){
		chromosomes = new Vector<GenomeSegments.chrID>();
		for(int i=1; i < tokens.length; i++)
		    chromosomes.add(GenomeSegments.getchrID(tokens[i]));
	    } else if(tokens[0].equalsIgnoreCase("HMM")){
		use_HMM = true;
		if(tokens.length > 1 && tokens[1].equalsIgnoreCase("WithDummy"))
		    HMM_dummy = true;
	    } else if(tokens[0].equalsIgnoreCase("oppositestrands")){
		if(tokens[1].toLowerCase().charAt(0) == 'y')
		    sequence_read.opposite_strands = true;
	    } else if(tokens[0].equalsIgnoreCase("outputformat")){
		if(tokens[1].toLowerCase().equalsIgnoreCase("fastq"))
		    fastq_format = true;
	    } else {
		System.err.println("Unknown parameter: " + tokens[0]);
		System.exit(0);
	    }
	}
	if(coverage <= 0 && !sequence_read.carpet_read){
	    System.err.println("coverage not specified");
	    System.exit(0);
	}
	if(fastq_format){
	    char q = (char)(33 + quality);
	    for(int i=0; i < sequence_read.read_length; i++)
		quality_string += q;
	}

	sequence_read.error_probability = Math.pow(10, -quality/10.);

	numbered = new GenomeSegments(args[1], false);
	XY = new GenomeSegments(args[1], true);

	if(SNPs_file != null){
	    SNPs = new Vector< List<GenomeSegments.SNP> >();
	    for(int i=0; i < GenomeSegments.chrID.values().length; i++)
		SNPs.add(new Vector<GenomeSegments.SNP>());
	    BufferedReader infile = new BufferedReader(new FileReader(SNPs_file));
	    for(String next = infile.readLine(); next != null; next = infile.readLine())
	    {
		GenomeSegments.SNP snp = new GenomeSegments.SNP(next);
		SNPs.get(snp.chromosome.ordinal()).add(snp);
	    }
	}

	int log_sizerange = (int)Math.ceil(Math.log(max_var_size / min_var_size) / Math.log(2));
	List<CNV> all_variations = new Vector<CNV>();
	CNV variation;

	if(num_common_CNVs > 0){
	    if(use_HMM){
		System.err.println("Can't use both HMM and common CNVs");
		System.exit(0);
	    }
	    numbered.sort();
	    XY.sort();
	    common_CNVs = new Vector<CNV>();
	    for(Iterator<String> iter = common_CNV_files.iterator(); iter.hasNext(); )
	    {
		BufferedReader cnv_file = new BufferedReader(new FileReader(iter.next()));
		cnv_file.readLine();
		for(String next = cnv_file.readLine(); next != null; next = cnv_file.readLine())
		{
		    variation = new CNV(next);
		    if((chromosomes == null || chromosomes.contains(variation.chromosome))
				&& variation.size > min_var_size)
			common_CNVs.add(variation);
		}
	    }
	}

	NavigableSet<CNV> variations = new TreeSet<CNV>(CNV.BySize);
	if(!use_HMM)
	for(CNV.CNV_type type : CNV.CNV_type.values())
	if(var_sizes[type.ordinal()] != null){
	  variations.clear();
	  boolean random_sizes = false;
	  for(int i=0; i < var_sizes[type.ordinal()].length; i++){    
	    int size = var_sizes[type.ordinal()][i];
	    boolean heterozygous = true;;
	    if(size == 0){
	      double size_rank = Math.floor(Math.pow(2,
		generator.nextFloat() * Math.log(log_sizerange+1) / Math.log(2)));
	      size = min_var_size * (int)Math.pow(2, size_rank-1);
	      random_sizes = true;
	      heterozygous = generator.nextFloat() > homozygous_probability;
	    } else if(size < 0){
	      heterozygous = false;
	      size = -size;
	    }
	    GenomeSegments.Segment segment = random_segment(heterozygous, size);
	    variation = new CNV(type, segment.chromosome, segment.length,
				size, heterozygous);
	    if(!segment.addVariation(variation)){
		i--;
		continue;
	    }
	    variations.add(variation);
	    all_variations.add(variation);
	  }
	  if(random_sizes){
	    int num_big = 0;
	    for(Iterator<CNV> iter = variations.iterator(); iter.hasNext() && ((variation = iter.next()).size > max_var_size/2 || num_big < min_big_vars); num_big++)
		variation.size = max_var_size;
	  }
	}

	variations = new TreeSet<CNV>(CNV.ByStart);
	for(int i=0; i < num_common_CNVs; i++){
	    variation = common_CNVs.get(generator.nextInt(common_CNVs.size()));
	    CNV previous_variation = variations.ceiling(variation);
	    if(variation.equals(previous_variation))
		previous_variation.heterozygous = false;
	    else {
		GenomeSegments.Segment segment = variation.find_segment();
		if(segment == null || !segment.addVariation(variation)){
		    i--;
		    continue;
		}
		variations.add(variation);
		all_variations.add(variation);
	    }
	}

	if(use_HMM){
	    PrintStream coverage_file =
		HMM_dummy ? new PrintStream(new FileOutputStream(
			output_name.substring(0,output_name.length()-12)
			+ "dummy.coverage")) : null;
	    all_variations.addAll(numbered.run_HMMs(coverage_file));
	    all_variations.addAll(XY.run_HMMs(coverage_file));
	}

	long num_reads = Math.round(coverage * (2*numbered.total_length + XY.total_length) /
				(sequence_read.read_length * 2));


	if(sequence_read.carpet_read){
	    for(int i=0; i < numbered.segments.length; i++)
		numbered.segments[i].set_carpet_read();
	    for(int i=0; i < XY.segments.length; i++)
		XY.segments[i].set_carpet_read();
	} else for(long i=0; i < num_reads; i++){
	    GenomeSegments.Segment segment = random_segment(true, 0);
	    segment.num_reads++;
	}

	List<Integer> segment_order = new Vector<Integer>();
	for(int i=0; i < numbered.segments.length + XY.segments.length; i++)
	    segment_order.add(i);
	Collections.shuffle(segment_order);

	for(Iterator<Integer> iter = segment_order.iterator(); iter.hasNext(); )
	    segmentAt(iter.next().intValue()).do_simulation();

	Collections.sort(all_variations, Collections.reverseOrder(CNV.ByStart));
	if(!all_variations.isEmpty()){
	PrintStream var_file = new PrintStream(new FileOutputStream(
		output_name.substring(0,output_name.length()-12)
			+ "variations"));
	for(Iterator<CNV> iter = all_variations.iterator(); iter.hasNext(); )
	    var_file.println(iter.next().toString());
	}
      } catch (Exception e) {
	e.printStackTrace();
      }
    }
}

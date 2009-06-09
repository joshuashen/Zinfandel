package cnv_HMM;
import cnv_simulation.*;
import java.io.*;
import java.util.*;

public class findLongInserts {
   
    static class insert implements Comparable<insert> {
	int start;
	int size;
	String name;
	insert(int st, int si, String n)
	{ start = st; size = si; name = n; }
	public int compareTo(insert o)
	{
	    int result = start + size - o.start - o.size;
	    if(result == 0)
		result = start - o.start;
	    if(result == 0)
		return name.compareTo(o.name);
	    return result;
	}
    }

    public static void main(String[] args) {
      try {
	GenomeSegments.chrID chromosome = null;
	File params = null, segments = null, mapping = null,
		error_distr = null, CNV_file = null;
	int orientation = 1;
	double logqvalthreshold = -3;
	boolean verbose = false;

	for(int i=0; i < args.length; i++){
	    String[] tokens = args[i].split("[=,]");
	    if(tokens.length == 1){
	      if(tokens[0].equalsIgnoreCase("verbose"))
		verbose = true;
	      else if(params == null)
		params = new File(tokens[0]);
	    } else if(tokens[0].equalsIgnoreCase("params"))
		params = new File(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("mapping"))
		mapping = new File(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("segments"))
		segments = new File(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("qvalthreshold"))
		logqvalthreshold = Math.log(Double.parseDouble(tokens[1]));
	    else if(tokens[0].equalsIgnoreCase("badmappingdistribution"))
		error_distr = new File(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("CNVs"))
		CNV_file = new File(tokens[1]);
	}

	if(params != null){
	  BufferedReader paramfile = new BufferedReader(new FileReader(params));
	  for(String paramline = paramfile.readLine(); paramline != null; paramline = paramfile.readLine())
	  {
	    String[] tokens = paramline.split("\\s+");
	    if(tokens[0].equalsIgnoreCase("chromosome") && tokens.length == 2)
		chromosome = GenomeSegments.getchrID(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("oppositestrands") && tokens[1].toLowerCase().charAt(0) == 'y')
		orientation = 2;
	    else if(tokens[0].equalsIgnoreCase("readblockdistribution"))
		sequence_read.readSizeDistribution(tokens[1]);
	    else if(tokens[0].equalsIgnoreCase("readlength"))
		sequence_read.read_length = Integer.parseInt(tokens[1]);
	  }
	  paramfile.close();
	}

	if(chromosome == null){
	    System.err.println("must specify one chromosome");
	    System.exit(1);
	} else if(sequence_read.block_size_distribution == null){
	    System.err.println("must specify insert size distribution");
	    System.exit(1);
	} else if(error_distr == null){
	    System.err.println("must specify insert size distribution");
	    System.exit(1);
	}

	double[] cumulative_block_size_distribution = new double[sequence_read.block_size_distribution.length];
	long total_count = 0, running_count = 0;
	for(int i=0; i < sequence_read.block_size_distribution.length; i++)
	    total_count += (sequence_read.min_block_size + i) * sequence_read.block_size_distribution[i];
	for(int i=0; i < sequence_read.block_size_distribution.length; i++){
	    running_count += (sequence_read.min_block_size + i) * sequence_read.block_size_distribution[i];
	    cumulative_block_size_distribution[i] = (double)running_count / total_count;
	}

	int max_block_size = sequence_read.min_block_size +
				sequence_read.block_size_distribution.length - 1;

	Map<Integer,Integer> error_sizes = new TreeMap<Integer,Integer>();
	BufferedReader error_distr_input = new BufferedReader(new FileReader(error_distr));
	total_count = running_count = 0;
	double error_prob = Double.parseDouble(error_distr_input.readLine());
	for(String next = error_distr_input.readLine(); next != null;
		   next = error_distr_input.readLine())
	{
	    int blank = next.indexOf(' '),
	        size = Integer.parseInt(next.substring(0,blank)),
		count = Integer.parseInt(next.substring(blank+1)) * size;
	    error_sizes.put(size, count);
	    total_count += count;
	}

	NavigableMap<Integer,Double> cumulative_error_distr = new TreeMap<Integer,Double>();
	for(Iterator< Map.Entry<Integer,Integer> > iter = error_sizes.entrySet().iterator(); iter.hasNext(); )
	{
	    Map.Entry<Integer,Integer> entry = iter.next();
	    running_count += entry.getValue().intValue();
	    cumulative_error_distr.put(entry.getKey(), (double)running_count / total_count);
	}

	boolean[] sequenced = new boolean[GenomeSegments.chr_length[chromosome.ordinal()]+2];
	Arrays.fill(sequenced, true);
	sequenced[0] = sequenced[sequenced.length-1] = false;
	int num_pos = sequenced.length - 2;

	if(segments != null){
	    GenomeSegments sequenced_segments = null;
	    Simulation.chromosomes = new Vector<GenomeSegments.chrID>();
	    Simulation.chromosomes.add(chromosome);
	    sequenced_segments = new GenomeSegments(segments.getPath(), !GenomeSegments.isNumbered(chromosome));
	    sequenced_segments.sort();
	    int current_segment = 0;
	    for(int next_pos = 1; next_pos < sequenced.length; next_pos++)
	    {
	      while(current_segment < sequenced_segments.segments.length-1 &&
		sequenced_segments.segments[current_segment].compareTo(chromosome, next_pos)<0)
		current_segment++;
	      if(sequenced_segments.segments[current_segment].compareTo(chromosome, next_pos) != 0)
	      {
		sequenced[next_pos] = false;
		num_pos--;
	      }
	    }
	}

	int[] pos_inserts = new int[sequenced.length];
	int[] quartile = new int[sequenced.length];
	double[] logpvalue = new double[sequenced.length];
	Arrays.fill(pos_inserts, 0);
	Arrays.fill(quartile, 0);
	Arrays.fill(logpvalue, 0);

	int current_pos = 0;
	LinkedList<insert> inserts = new LinkedList<insert>();
	int[] sizes = new int[1000];

	BufferedReader infile = new BufferedReader(mapping == null
				? new InputStreamReader(System.in)
				: new FileReader(mapping));

	boolean finished = false;
	while(!finished){
	  String next = infile.readLine(); 
	  String[] tokens = next == null ? null : next.split("\\t");
	  if(current_pos == 0){
	    if(tokens == null){
	      System.err.println("no mapping data for " + chromosome);
	      System.exit(1);
	    }
	    if(!tokens[1].equals(chromosome.name()))
	      continue;
	  }
	  int flag = tokens == null ? orientation : Integer.parseInt(tokens[5]);
	  boolean forward = tokens == null || tokens[3].charAt(0) == '+';
	  if(flag != orientation && flag != orientation + 16 || !forward)
	    continue;
	  finished = tokens == null || !tokens[1].equals(chromosome.name());
	  int pos = finished ? Integer.MAX_VALUE : Integer.parseInt(tokens[2]);
	  int size = finished ? 0 : Integer.parseInt(tokens[4]);
	  for(;current_pos < pos && !inserts.isEmpty(); current_pos++){
	    int num_inserts = 0;
	    for(Iterator<insert> iter = inserts.iterator(); iter.hasNext(); ){
		insert ins = iter.next();
		sizes[num_inserts++] = ins.size;
		if(ins.start + ins.size == current_pos+1)
		    iter.remove();
	    }
	    if(!sequenced[current_pos])
		continue;
	    Arrays.sort(sizes, 0, num_inserts);
	    int quartile_index = num_inserts*3/4;
	    quartile[current_pos] = sizes[quartile_index];
	    pos_inserts[current_pos] = num_inserts;
	    if(quartile[current_pos] > sequence_read.min_block_size){
		double prob_cond_good =
			quartile[current_pos] > max_block_size ? 1. :
			cumulative_block_size_distribution[quartile[current_pos] - sequence_read.min_block_size -1];
		double prob_cond_bad = cumulative_error_distr.lowerEntry(quartile[current_pos]).getValue().doubleValue();
		double prob = error_prob * prob_cond_bad + (1 - error_prob) * prob_cond_good;
		logpvalue[current_pos] = logbinocdf.compute(quartile_index, num_inserts, prob);
		/*if(quartile[current_pos] > 8000000 && logpvalue[current_pos] <= 0)
		    System.err.println("quartile " + quartile[current_pos] + " out of "
				+ num_inserts + " prob good " + prob_cond_good
				+ " bad " + prob_cond_bad + " creates success prob "
				+ prob + " logpvalue " + logpvalue[current_pos]);
		*/
	    }
	  }
	  if(size > sequence_read.read_length)
	    inserts.add(new insert(current_pos = pos, size, tokens[0]));
	}

	// System.err.println("starting sort of length " + logpvalue.length + ' ' + (new java.util.Date()).toString());
	double[] logpval_set = logpvalue.clone();
	Arrays.sort(logpval_set);
	// System.err.println("ending sort " + (new java.util.Date()).toString());
	double log_num_pos = Math.log(num_pos);
	TreeMap<Double,Double> logqval_map = new TreeMap<Double,Double>();
	for(int i=1; i <= logpval_set.length; i++)
	if(i == logpval_set.length || logpval_set[i] != logpval_set[i-1])
	    logqval_map.put(logpval_set[i-1], logpval_set[i-1] + log_num_pos - Math.log(i));
	// System.err.println("qval map computed " + (new java.util.Date()).toString());
	double[] logqvalue = logpval_set;
	for(int i=0; i < logpvalue.length; i++)
	    logqvalue[i] = logqval_map.get(logpvalue[i]).doubleValue();

	List<CNV> inferred_CNVs = new Vector<CNV>();
	CNV cnv;
	for(int index = 0;;index++){
	    for(;index < logqvalue.length && logqvalue[index]>logqvalthreshold; index++);
	    if(index >= logqvalue.length)
		break;
	    int start = index++;
	    int[] num_insert_range = {pos_inserts[start], pos_inserts[start]};
	    int[] quartile_range = {quartile[start], quartile[start]};
	    double[] logpval_range = {logpvalue[start], logpvalue[start]};
	    double[] logqval_range = {logqvalue[start], logqvalue[start]};
	    for(int count = 1; index < logqvalue.length && (float)count/(index-start)>2./3; index++)
	    if(logqvalue[index] <= logqvalthreshold){
		count++;
		num_insert_range[0] = Math.min(num_insert_range[0], pos_inserts[index]);
		quartile_range[0] = Math.min(quartile_range[0], quartile[index]);
		logpval_range[0] = Math.min(logpval_range[0], logpvalue[index]);
		logqval_range[0] = Math.min(logqval_range[0], logqvalue[index]);
		num_insert_range[1] = Math.max(num_insert_range[1], pos_inserts[index]);
		quartile_range[1] = Math.max(quartile_range[1], quartile[index]);
		logpval_range[1] = Math.max(logpval_range[1], logpvalue[index]);
		logqval_range[1] = Math.max(logqval_range[1], logqvalue[index]);
	    }
	    while(logqvalue[--index] > logqvalthreshold);
	    
	    cnv = new CNV(CNV.CNV_type.DEL, chromosome, index-start+2, index-start+1, true);
	    cnv.absolute_start = start;
	    inferred_CNVs.add(cnv);
	    if(verbose)
		System.out.println(String.valueOf(start) + '-' + index + " size " + (index-start+1) + " quartile " + quartile_range[0] + '-' + quartile_range[1] + " out of " + num_insert_range[0] + '-' + num_insert_range[1] + " logpvalue " + logpval_range[0] + '-' + logpval_range[1] + " logqvalue " + logqval_range[0] + '-' + logqval_range[1]);
	}

	if(CNV_file != null){
	    boolean[] non_deleted = sequenced;
	    int num_non_deleted = num_pos, false_positives = 0;
	    List<CNV> CNVs = new Vector<CNV>();
	    BufferedReader CNVfile = new BufferedReader(new FileReader(CNV_file));
	    for(String next = CNVfile.readLine(); next != null; next = CNVfile.readLine())
		if((cnv = new CNV(next)).chromosome == chromosome && cnv.type == CNV.CNV_type.DEL)
		{
		    CNVs.add(cnv);
		    for(int i=Math.max(0,cnv.absolute_start-max_block_size); i < Math.min(non_deleted.length, cnv.absolute_start+cnv.size+max_block_size); i++)
		    {
			non_deleted[i] = false;
			num_non_deleted--;
		    }
		}
	    CNVfile.close();
	    for(int i=0; i < logqvalue.length; i++)
		if(non_deleted[i] && logqvalue[i] <= logqvalthreshold)
		    false_positives++;
	    if(verbose)
		System.out.println(String.valueOf(false_positives) + " false positives, " + (float)false_positives/num_non_deleted*100 + '%');
	    for(Iterator<CNV> iter = CNVs.iterator(); iter.hasNext(); ){
		cnv = iter.next();
		int count = 0;
		for(int i=0; i < cnv.size; i++)
		  if(logqvalue[cnv.absolute_start+i] <= logqvalthreshold)
		    count++;
		Arrays.sort(logpvalue, cnv.absolute_start, cnv.absolute_start+cnv.size);
		Arrays.sort(logqvalue, cnv.absolute_start, cnv.absolute_start+cnv.size);
		if(verbose)
		     System.out.println(cnv.toString() + " size " + cnv.size + ": " + count + " true positives, " + (float)count/cnv.size*100
+ "%, range " + logpvalue[cnv.absolute_start] + " - " + logpvalue[cnv.absolute_start+cnv.size-1] + " 90% at " + logpvalue[cnv.absolute_start+cnv.size*9/10] + " 95% at " + logpvalue[cnv.absolute_start+cnv.size*19/20]
+ ", qvalue range " + logqvalue[cnv.absolute_start] + " - " + logqvalue[cnv.absolute_start+cnv.size-1] + " 90% at " + logqvalue[cnv.absolute_start+cnv.size*9/10] + " 95% at " + logqvalue[cnv.absolute_start+cnv.size*19/20]);
	    }
	    runViterbi.compare_CNV_lists(CNVs, inferred_CNVs);
	}
      } catch (Exception e) {
	e.printStackTrace();
      }
    }

}

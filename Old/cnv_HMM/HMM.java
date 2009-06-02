package cnv_HMM;
import cnv_simulation.*;
import java.io.*;
import java.util.*;
import org.apache.commons.math.distribution.PoissonDistributionImpl;

public class HMM {
    private GenomeSegments.chrID chr;
    private double[][] transition;
    private double[][][][] emission;
    private double[][] N_emission;
    int[][] expectation = null;
    private int current_pos;
    public int normal_state;
    private int num_states;
    private double best_path_log_probability[][];
    private int best_path_CNV_start[];
    private Vector<CNV> best_path_CNVs;
    private GenomeSegments segments;
    private int current_state;
    private static int max_emission = 50;
    static double distribution_threshold = .0001;
    static double EM_change_threshold = .01;
    static double EM_max_iterations = 10;
    static boolean infer_emission_distribution = false;
    static boolean infer_transition_distribution = false;
    static boolean uniform_start = false;
    int[][] data;
    int[] read_lengths;
    boolean[] sequenced;
    float[] average_cov;

    static boolean debug = false;

    public HMM(GenomeSegments.chrID chr, long space_between_duplications, long space_between_deletions, int average_CNV_length, double homozygous_probability, int[] read_lengths, File[][] carpet, GenomeSegments segments, int[][] data)
    throws Exception
    {
	this.chr = chr;
	this.segments = segments;
	this.data = data;
	this.read_lengths = read_lengths;
	current_pos = 1;
	sequenced = new boolean[GenomeSegments.chr_length[chr.ordinal()]+2];
	Arrays.fill(sequenced, true);
	sequenced[0] = false;
	if(segments != null){
	  int current_segment = Arrays.binarySearch(segments.segments,
		new GenomeSegments.Segment(chr, 1));
	  if(current_segment < 0)
	    current_segment = - current_segment - 1;
	  for(int next_pos = 1; next_pos < sequenced.length; next_pos++)
	  {
	    while(current_segment < segments.segments.length-1 &&
		segments.segments[current_segment].compareTo(chr, next_pos)<0)
		current_segment++;
	    if(segments.segments[current_segment].compareTo(chr, next_pos) != 0)
		sequenced[next_pos] = false;
	  }
	}
	int max_expectation = 1;
	int[] normal_expectation = new int[read_lengths.length];
	Arrays.fill(normal_expectation, 1);
	if(carpet != null){
	  expectation = new int[carpet.length][];
	  Arrays.fill(expectation, null);
	  int default_read_length = -1;
	  while(carpet[++default_read_length] == null)
	    expectation[default_read_length] = null;
	  for(int i = default_read_length; expectation[i] == null; i = (i+1) % carpet.length)
	  if(carpet[i] == null){
	    expectation[i] = expectation[default_read_length];
	    normal_expectation[i] = normal_expectation[default_read_length];
	  } else {
	    expectation[i] = runViterbi.readPileup(carpet[i], chr, read_lengths[i]);
	    max_expectation = Math.max(max_expectation, expectation[i][0]);
	    normal_expectation[i] = carpet[i].length*2;
	  }
	}
	current_state = normal_state = GenomeSegments.isNumbered(chr) ? 2 : 1;
	int distribution_max = 0;
	average_cov = new float[read_lengths.length];
	for(int i=0; i < read_lengths.length; i++)
	if(data[i] != null){
	  long total_sum = 0;
	  int total_count = 0;
	  for(int next_pos = 1; next_pos < data[i].length; next_pos++)
	    if(sequenced[next_pos]){
		total_sum += data[i][next_pos];
		total_count++;
	    }
	  average_cov[i] = (float)total_sum/total_count;
	  System.err.println("read length " + read_lengths[i] + " average coverage " + average_cov[i]/normal_state + " * " + read_lengths[i] + " = " + average_cov[i]*read_lengths[i]/normal_state);
	}

	num_states = normal_state * 2 + 1;
	best_path_log_probability = new double[2][];
	best_path_log_probability[0] = new double[num_states];
	best_path_log_probability[1] = new double[num_states];
	Arrays.fill(best_path_log_probability[0], Double.NEGATIVE_INFINITY);
	best_path_log_probability[0][normal_state] = 0;
	best_path_CNVs = new Vector<CNV>();
	best_path_CNV_start = new int[num_states];
	transition = new double[num_states][];
	N_emission = new double[num_states][];
	emission = new double[read_lengths.length][][][];
	for(int rl=0; rl < read_lengths.length; rl++)
	if(data[rl] != null){
	  emission[rl] = new double[max_expectation+1][][];
	  for(int expec=0; expec <= max_expectation; expec++){
	    emission[rl][expec] = new double[num_states][];
	    for(int i=0; i < num_states; i++){
	      emission[rl][expec][i] = new double[max_emission];
	      float weight_factor = i * (expec == 0 ? (float).2 : (float)expec) /
                        (normal_state * normal_expectation[rl]);
	      if(i == 0){
		Arrays.fill(emission[rl][expec][i], Double.NEGATIVE_INFINITY);
		emission[rl][expec][i][0] = 0;
	      } else {
	        PoissonDistributionImpl distribution = new PoissonDistributionImpl(
			average_cov[rl] * weight_factor);
	        for(int j=0; j < max_emission; j++)
		  emission[rl][expec][i][j] = Math.log(distribution.probability(j));
	      }
	    }
	  }
	}
	for(int i=0; i < num_states; i++){
	    N_emission[i] = new double[max_emission];
	    Arrays.fill(N_emission[i], i == normal_state ? Math.log(1./max_emission) : Double.NEGATIVE_INFINITY);
	    transition[i] = new double[num_states];
	    if(i == normal_state){
		transition[i][i] = Math.log(1. - 1. / space_between_duplications - 1. / space_between_deletions);
		if(i == 1){
		    transition[i][0] = Math.log(1. / space_between_deletions);
		    transition[i][2] = Math.log(1. / space_between_duplications);
		} else {
		    transition[i][0] = Math.log(homozygous_probability / space_between_deletions);
		    transition[i][4] = Math.log(homozygous_probability / space_between_duplications);
		    transition[i][1] = Math.log((1 - homozygous_probability) / space_between_deletions);
		    transition[i][3] = Math.log((1 - homozygous_probability) / space_between_duplications);
		}
	    } else {
		Arrays.fill(transition[i], Double.NEGATIVE_INFINITY);
		double continue_probability = Math.pow(.5, 1. / average_CNV_length);
		transition[i][i] = Math.log(continue_probability);
		transition[i][normal_state] = Math.log(1 - continue_probability);
	    }
	}

	/*if(infer_transition_distribution){
	  double EM_change = 1;
	  int num_iterations = -1;
	  if(uniform_start){
	    for(int i=0; i < num_states; i++)
	    if(i == normal_state)
		Arrays.fill(transition[i], 1. / num_states);
	    else {
		Arrays.fill(transition[i], 0);
		transition[i][i] = transition[i][normal_state] = .5;
	    }
	  }
	  while(++num_iterations < EM_max_iterations && EM_change > EM_change_threshold)
	  {
	  if(debug){
	    System.err.println("At " + (new java.util.Date()).toString() + " transition probabilities before iteration "
		+ num_iterations + ", change " + EM_change + ':');
	    for(int i=0; i < num_states; i++)
	      System.err.println("For " + i + ": " + Arrays.toString(transition[i]));
	  }
	  EM_change = EM_iteration();
	  }
	}
	*/

	if(debug){
	  System.err.println("Starting HMM run at " + (new java.util.Date()).toString());
	  /* System.err.println(chr.toString() + " emission probabilities:");
	  for(int expec=0; expec < emission.length; expec++)
	  if(emission[expec] != null){
	   System.err.println("for expectation " + expec + ':');
	   for(int i=0; i < num_states; i++)
	    System.err.println("For " + i + ": " + Arrays.toString(emission[expec][i]));
	  }
	  */
	  System.err.println("transition probabilities:");
	  for(int i=0; i < num_states; i++)
	    System.err.println("For " + i + ": " + Arrays.toString(transition[i]));
	}
	for(int next_pos = 1; next_pos < data[0].length; next_pos++)
	    observed_emission();
    }

    private void observed_emission()
    throws Exception
    {
      int state = -30, best_predecessor = -30;
      int index = current_pos % 2;
      int e = -1;
      try {
        Arrays.fill(best_path_log_probability[index], Double.NEGATIVE_INFINITY);
	boolean possible = false;
	for(int s = -1; s < num_states; s++){
	    if(s == normal_state)
		continue;
	    state = s < 0 ? normal_state : s;
	    best_predecessor = normal_state;
	    for(int predecessor = 0;  predecessor < num_states; predecessor++){
		double log_prob = best_path_log_probability[1 - index][predecessor]
			+ transition[predecessor][state];
		for(int i=0; i < data.length; i++)
		if(data[i] != null){
		    double[][] emission_to_use =
		    sequenced[current_pos]
			? emission[i][expectation == null ? 1 : expectation[i][current_pos]]
			: N_emission;
		    log_prob += emission_to_use[state][e = data[i][current_pos]];
		}
		if(log_prob > best_path_log_probability[index][state]){
		    best_path_log_probability[index][state] = log_prob;
		    best_predecessor = predecessor;
		    possible = true;
		}
	    }
	    if(state == normal_state && best_predecessor != normal_state){
		while(!best_path_CNVs.isEmpty() &&
			best_path_CNVs.lastElement().absolute_start + best_path_CNVs.lastElement().size
				>= best_path_CNV_start[best_predecessor])
		    best_path_CNVs.remove(best_path_CNVs.size() - 1);
		best_path_CNVs.add(new CNV(
			best_predecessor > normal_state ? CNV.CNV_type.DUP : CNV.CNV_type.DEL,
			chr, 2, 1, best_predecessor==1 || best_predecessor==3));
		best_path_CNVs.lastElement().absolute_start = best_path_CNV_start[best_predecessor];
		best_path_CNVs.lastElement().size = current_pos - best_path_CNV_start[best_predecessor];
	    } else if(state != normal_state && best_predecessor == normal_state)
		best_path_CNV_start[state] = current_pos;
	}
	if(!possible)
	    throw new Exception("Impossible emission");
	current_pos++;
	/*
	if(debug && (current_pos%1000000000 == 0 || current_pos >= 15400856 && current_pos <= 15410000)){
	    System.err.println((new java.util.Date()).toString() + " at position " + current_pos + " emission " + e + " best probs " + Arrays.toString(best_path_log_probability[index]) 
		+ " starting points " + Arrays.toString(best_path_CNV_start));
	}
	*/
      } catch (Exception ex) {
	System.err.println("at position " + current_pos + ", state " + state
		+ ", predecessor " + best_predecessor);
	System.err.println("best probs " + Arrays.toString(best_path_log_probability[index]));
	throw ex;
      }
    }

    Collection<CNV> optimalSolution() {
	return best_path_CNVs;
    }

    public int make_transition() {
	double p = Simulation.generator.nextDouble();
	int next_state = -1;
	while(p >= 0)
	    p -= transition[current_state][next_state];
	return current_state = next_state;
    }

    public int make_emission() {
	double p = Simulation.generator.nextDouble();
	int result = -1;
	while(p >= 0)
	    p -= emission[0][1][current_state][++result];
	return result;
    }

    /*
    private double EM_iteration() 
    throws IOException {
	File forward_prob_file = File.createTempFile("probs",null, new File("."));
	double[][] probs = new double[2][];
	double[] normal_probs = new double[num_states];
	Arrays.fill(normal_probs, 0);
	normal_probs[normal_state] = 1;
	probs[0] = new double[num_states];
	probs[1] = new double[num_states];
	double[][] reverse_transition = new double[num_states][];
	for(int j = 0; j < num_states; j++){
	    double total_prob = 0;
	    for(int k = 0; k < num_states; k++)
		total_prob += transition[k][j];
	    reverse_transition[j] = new double[num_states];
	    if(total_prob == 0)
		Arrays.fill(reverse_transition[j], 1. / num_states);
	    else for(int k = 0; k < num_states; k++)
		reverse_transition[j][k] = transition[k][j] / total_prob;
	}
	if(debug){
	    System.err.println("reverse transition probabilities:");
	    for(int i=0; i < num_states; i++)
		System.err.println("For " + i + ": " + Arrays.toString(reverse_transition[i]));
	}
	DataOutputStream output = new DataOutputStream(new FileOutputStream(forward_prob_file));
	for(int i = 0, index = 0; i < data.length; i++, index = 1 - index)
	{
	  if(sequenced[i]){
	    Arrays.fill(probs[index], 0);
	    double total_prob = 0;
	    for(int j = 0; j < num_states; j++)
	    for(int k=0; k < num_states; k++)
		probs[index][j] += probs[1 - index][k] * transition[k][j];
	    for(int j=0; j < num_states; j++)
		total_prob += probs[index][j] *= emission[expectation[i]][j][data[i]];
	    for(int j=0; j < num_states; j++)
		probs[index][j] /= total_prob;
	  } else
	    System.arraycopy(normal_probs, 0, probs[index], 0, num_states);
	  for(int j=0; j < num_states; j++)
	    output.writeDouble(probs[index][j]);
	}
	output.close();

	double[] total_per_state = new double[num_states];
	Arrays.fill(total_per_state, 0);
	double[][] new_transition = new double[num_states][];
	for(int j=0; j < num_states; j++)
	    new_transition[j] = total_per_state.clone();

	RandomAccessFile forward_input = new RandomAccessFile(forward_prob_file, "r");
	int location_size = num_states * Double.SIZE / Byte.SIZE;
	long forward_offset = (long)(data.length - 1) * location_size;
	for(int i = data.length, index = 0; i > 0;
		i--, index = 1 - index, forward_offset -= location_size)
	{
	  if(sequenced[i]){
	    Arrays.fill(probs[index], 0);
	    double total_prob = 0;
	    for(int j = 0; j < num_states; j++)
	    for(int k=0; k < num_states; k++)
		probs[index][j] += probs[1 - index][k] * reverse_transition[k][j];
	    for(int j=0; j < num_states; j++)
		total_prob += probs[index][j] *= emission[expectation[i]][j][data[i]];
	    for(int j=0; j < num_states; j++)
		probs[index][j] /= total_prob;
	    if(sequenced[i-1]){
	      forward_input.seek(forward_offset);
	      for(int j = 0; j < num_states; j++){
	       double forward_prob = forward_input.readDouble();
	       for(int k = 0; k < num_states; k++)
	       if(j == normal_state || k == normal_state || j == k)
	       {
		double kj_prob = probs[index][k] * forward_prob;
		new_transition[j][k] += kj_prob;
		total_per_state[j] += kj_prob;
	       }
	      }
	    }
	  } else
	    System.arraycopy(normal_probs, 0, probs[index], 0, num_states);
	}
	forward_input.close();
	forward_prob_file.delete();

	double highest_change = 0;
	for(int j = 0; j < num_states; j++)
	for(int k = 0; k < num_states; k++){
	    new_transition[j][k] /= total_per_state[j];
	    if(new_transition[j][k] != 0 || transition[j][k] != 0)
	      if(transition[j][k] == 0 || new_transition[j][k] == 0)
		highest_change = 1;
	      else
		highest_change = Math.max(highest_change,
			Math.abs(new_transition[j][k] - transition[j][k])
				/ transition[j][k]);
	}
	transition = new_transition;
	return highest_change;
    }
    */
}

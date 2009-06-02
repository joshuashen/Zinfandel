package cnv_simulation;
import java.io.*;
import java.util.*;

public class sequence_read {

    String[] read = new String[2];
    int start;
    int block_size;

    public String toString() {
	return String.valueOf(start) + ' ' + block_size + ' '
		+ read[0] + " ----- " + read[1];
    }

    static boolean carpet_read = false;
    public static int read_length = 35;
    static int mean_block_size = 4000;
    static float std_block_size = 0;
    static boolean opposite_strands = false;
    static double error_probability;
    public static int[] block_size_distribution = null;
    public static int min_block_size = -1;
    static int block_size_sample_size = 0;

    static char[] letters = {'A', 'C', 'G', 'M', 'N', 'R', 'T'};
    static char[] opposite = {'T', 'G', 'C', 'R', 'N', 'M', 'A'};
    static String [] alternates = {"CGT", "AGT", "ACT", null, null, null, "ACG"};

    public static void readSizeDistribution(String file)
    throws IOException
    {
	block_size_distribution = new int[5000];
	Arrays.fill(block_size_distribution, 0);
	BufferedReader infile = new BufferedReader(new FileReader(file));
	for(String next = infile.readLine(); next != null; next = infile.readLine()){
	    String[] tokens = next.split("\\t");
	    if(min_block_size < 0)
		min_block_size = Integer.parseInt(tokens[0]);
	    int val = Integer.parseInt(tokens[1]);
	    block_size_distribution[Integer.parseInt(tokens[0]) - min_block_size] = val;
	    block_size_sample_size += val;
	}
    }

    static sequence_read perform_read(char[] ref_sequence, String sequence, int num_read)
    {
	sequence_read result = new sequence_read();
	if(block_size_distribution == null){
	  result.block_size = mean_block_size;
	  if(std_block_size > 0)
	    result.block_size += Math.round(Simulation.generator.nextGaussian() * std_block_size);
	} else {
	  int size_index = -1;
	  for(int num = Simulation.generator.nextInt(block_size_sample_size); num >= 0; num -= block_size_distribution[++size_index]);
	  result.block_size = min_block_size + size_index;
	}
	result.start = carpet_read ? num_read :
		Simulation.generator.nextInt(
		(sequence == null ? ref_sequence.length : sequence.length()) - result.block_size);
	if(sequence == null){
	    result.read[0] = new String(ref_sequence, result.start, read_length);
	    result.read[1] = new String(ref_sequence, result.start + result.block_size - read_length, read_length);
	} else {
	    result.read[0] = sequence.substring(result.start, result.start+read_length);
	    result.read[1] = sequence.substring(result.start + result.block_size - read_length, result.start + result.block_size);
	}
	try {
	if(opposite_strands)
	    result.read[1] = complement(result.read[1]);
		
	for(int i = 0; i < 2; i++)
	for(int s = 0; s < read_length; s++){
	    double r = Simulation.generator.nextDouble();
	    if(r >= error_probability)
		continue;
	    String alt = alternates[Arrays.binarySearch(letters, result.read[i].charAt(s))];
	    if(alt == null)
		continue;
	    result.read[i] = result.read[i].substring(0, s)
		+ alt.charAt((int)Math.floor(r*3/error_probability))
		+ result.read[i].substring(s+1);
	}
	return result;
	} catch (RuntimeException e) {
	    System.err.println("While processing " + result.toString());
	    throw e;
	}
    }

    static String complement(String source)
    {
	String result = "";
	for(int s=source.length()-1; s >= 0; s--)
	  result += opposite[Arrays.binarySearch(letters, source.charAt(s))];
	return result;
    }
}

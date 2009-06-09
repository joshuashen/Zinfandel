package cnv_HMM;
import cnv_simulation.*;
import java.io.*;
import java.util.*;

public class runViterbi {
   
  static class DividedCNV extends Vector<CNV> {
	  DividedCNV(CNV c) { add(c); }

	  public String toString() {
	    return size() == 1 ? get(0).toString() : super.toString();
	  }

	  float overlap(DividedCNV other)
	  {
	    float overlap1 = 0, overlap2 = 0;
	    for(Iterator<CNV> iter = iterator(); iter.hasNext(); )
		    overlap1 += other.get(0).overlap(iter.next());
	    for(Iterator<CNV> iter = other.iterator(); iter.hasNext(); )
		    overlap2 += get(0).overlap(iter.next());
	    return Math.max(overlap1, overlap2);
	  }
  }

  static int report_count_threshold = Integer.MAX_VALUE;
  static Map<GenomeSegments.chrID, HMM> hmm_map = new HashMap<GenomeSegments.chrID, HMM>();
    
  static int default_read_length = 35;
  static float overlap_threshold = (float).5;

  public static void main(String[] args) {
    try {
	    int num_duplications = 0;
	    int num_deletions = 0;
	    int average_CNV_length = 80000;
    	double homozygous_probability = .000045;
  	  int default_carpet_read_length = 35;
    	Collection<GenomeSegments.chrID> chromosomes = null;
    	File params = null, segments = null, actual_file = null;
	    File[][] carpet = null, pileup = null;
	    
 	    Map<Integer, List<File> > carpet_map = null, pileup_map = null;
// pileup_map is a holder for depth-coverage data, specifically, the program reads a file similar to mapview format from Maq, where each read is a line, and the start position of the mapping is specified at the third column.
// mapview format example:
//   chr18_d_9_151_a0ac8/1	chr18	  9	  +	  143	18	... 
//    readName           refName    startPositon  direction distance flag ....


//parse argument 
    	for(int i=0; i < args.length; i++){
	      String[] tokens = args[i].split("[=,]");
	      if(tokens.length == 1) {
          if(params == null)
		        params = new File(tokens[0]);
	      } else if(tokens[0].equalsIgnoreCase("params"))
		      params = new File(tokens[1]);

    // get mapping file (in mapview-like format)
	      else if(tokens[0].equalsIgnoreCase("mapview") || tokens[0].equalsIgnoreCase("mapping")) {
  		    pileup_map = new TreeMap<Integer, List<File> >();
		      int current_read_size = default_read_length;
		      for(int j=1; j < tokens.length; j++)
		      try {
		        current_read_size = Integer.parseInt(tokens[j]);
		      } catch(NumberFormatException e) {
		        List<File> file_list = pileup_map.get(current_read_size);
		        if(file_list == null) {
			        file_list = new Vector<File>();
			        pileup_map.put(current_read_size, file_list);
		        }
		        file_list.add(new File(tokens[j]));
		      }
	      }  

    // reference segments	  
 	      else if(tokens[0].equalsIgnoreCase("segments"))
		      segments = new File(tokens[1]);
		      
     // true CNV list, useful in simulation where the CNVs are all known		      
	      else if(tokens[0].equalsIgnoreCase("actual"))
      		actual_file = new File(tokens[1]);
      		
      // carpet mapping, useful for background normalization       		
	      else if(tokens[0].equalsIgnoreCase("carpet")){
      		carpet_map = new TreeMap<Integer, List<File> >();
		      int current_read_size = default_carpet_read_length;
		      for(int j=1; j < tokens.length; j++)
		        try {
		          current_read_size = Integer.parseInt(tokens[j]);
		        } catch(NumberFormatException e) {
		        List<File> file_list = carpet_map.get(current_read_size);
		        if(file_list == null){
			        file_list = new Vector<File>();
			        carpet_map.put(current_read_size, file_list);
		        }
		        file_list.add(new File(tokens[j]));
		      }
	      }
	    }

	    if(pileup_map == null){
	      System.err.println("mapping input not specified");
	      System.exit(1);
	    }

// initialize file handles? 
	    Set<Integer> size_set = new TreeSet<Integer>(pileup_map.keySet());
	    if(carpet_map != null)
	      size_set.addAll(carpet_map.keySet());

	    int[] sizes = new int[size_set.size()];
	    
	    // file handles for mapping output? 
	    pileup = new File[size_set.size()][];
	    if(carpet_map != null)
	      carpet = new File[size_set.size()][];

	    int index = 0;
	    for(Iterator<Integer> iter = size_set.iterator(); iter.hasNext(); index++) {
	      sizes[index] = iter.next().intValue();
	      List<File> pileup_files = pileup_map.get(sizes[index]);
	      if(pileup_files == null)
		      pileup[index] = null;
	      else {
		      pileup[index] = new File[pileup_files.size()];
		      pileup_files.toArray(pileup[index]);
	      }
	      
	      if(carpet_map != null){
	        List<File> carpet_files = carpet_map.get(sizes[index]);
	        if(carpet_files == null)
		        carpet[index] = null;
	        else {
		        carpet[index] = new File[carpet_files.size()];
		        carpet_files.toArray(carpet[index]);
	        }
	      }
	    }

// parse params 
	    if(params != null){
	      BufferedReader paramfile = new BufferedReader(new FileReader(params));
	      for(String paramline = paramfile.readLine(); paramline != null; paramline = paramfile.readLine()) {
  	      String[] tokens = paramline.split("\\s+");
	        if(tokens[0].equalsIgnoreCase("amplifications") || tokens[0].equalsIgnoreCase("duplications"))
		        num_duplications = Integer.parseInt(tokens[1]);
	        else if(tokens[0].equalsIgnoreCase("deletions"))
		        num_deletions = Integer.parseInt(tokens[1]);
	        else if(tokens[0].equalsIgnoreCase("averageCNVLength"))
		        average_CNV_length = Integer.parseInt(tokens[1]);
	        else if(tokens[0].equalsIgnoreCase("homozygouspercentage"))
		        homozygous_probability = Double.parseDouble(tokens[1]) / 100;
	        else if(tokens[0].equalsIgnoreCase("chromosome")){
		        if(tokens[1].equalsIgnoreCase("ALL"))
		          chromosomes = Arrays.asList(GenomeSegments.chrID.values());
		        else {
		          chromosomes = new Vector<GenomeSegments.chrID>();
		          for(int i=1; i < tokens.length; i++)
			          chromosomes.add(GenomeSegments.getchrID(tokens[i]));
		        }
	        } 
	        else if(tokens[0].equalsIgnoreCase("readlength"))
		        default_read_length = Integer.parseInt(tokens[1]);
	        else if(tokens[0].equalsIgnoreCase("carpetreadlength"))
		        default_carpet_read_length = Integer.parseInt(tokens[1]);
	    /*else if(tokens[0].equalsIgnoreCase("infertransitiondistribution")){
		HMM.infer_transition_distribution = true;
		HMM.uniform_start = tokens.length > 1;
	    }*/
	        else if(tokens[0].equalsIgnoreCase("distributionthreshold"))
		        HMM.distribution_threshold = Double.parseDouble(tokens[1]);
	        else if(tokens[0].equalsIgnoreCase("overlapthreshold"))
		        overlap_threshold = Float.parseFloat(tokens[1]) / 100;
	        else if(tokens[0].equalsIgnoreCase("reportcountthreshold"))
		        report_count_threshold = Integer.parseInt(tokens[1]);
	        else if(tokens[0].equalsIgnoreCase("debug"))
		        HMM.debug = true;
	        else {
		        System.err.println("Unknown parameter: " + tokens[0]);
		        System.exit(0);
	        }
	      }
	      paramfile.close();
	    }

// if no chromosomes specified, read the mapping file to figure out
	    if(chromosomes == null){
	      chromosomes = new HashSet<GenomeSegments.chrID>();
	      BufferedReader pileupfile = new BufferedReader(new FileReader(pileup_map.values().iterator().next().get(0)));
    	  for(String next = pileupfile.readLine(); next != null; next = pileupfile.readLine())
    	    chromosomes.add(GenomeSegments.chrID.valueOf(GenomeSegments.chrID.class, next.split("\\s")[1]));
    	  pileupfile.close();
	    }

// set average indel size for estimation of transition probabilities
	    long total_length = 0;
	    for(Iterator<GenomeSegments.chrID> iter = chromosomes.iterator(); iter.hasNext(); )
	      total_length += GenomeSegments.chr_length[iter.next().ordinal()];

    	long space_between_duplications = 15000000, space_between_deletions = 15000000;
	    if(num_duplications + num_deletions > 0){
	      space_between_duplications = num_duplications == 0 ? Integer.MAX_VALUE : total_length/num_duplications;
	      space_between_deletions = num_deletions == 0 ? Integer.MAX_VALUE : total_length/num_deletions;
	    }

	    GenomeSegments numbered_segments = null, XY_segments = null;
	    if(segments != null){
	      Simulation.chromosomes = chromosomes;
	      numbered_segments = new GenomeSegments(segments.getPath(), false);
	      numbered_segments.sort();
	      XY_segments = new GenomeSegments(segments.getPath(), true);
	      XY_segments.sort();
	    }

// reference CNVs
    	List<CNV> actual = null;
    	if(actual_file != null){
	      actual = new Vector<CNV>();
	      BufferedReader actfile = new BufferedReader(new FileReader(actual_file));
	      for(String next = actfile.readLine(); next != null; next = actfile.readLine())
		      actual.add(new CNV(next));
	    }

// Decode via HMM  
	    List<CNV> optimal_solution = new Vector<CNV>();


/*
Mapping should be a separate class
*/

    	for(Iterator<GenomeSegments.chrID> iter = chromosomes.iterator(); iter.hasNext(); )	{
	      GenomeSegments.chrID current_chr = iter.next();
	      int[][] pileup_data = new int[pileup.length][];
	      for(int i=0; i < pileup.length; i++)
	      // read mapping file 
		      pileup_data[i] = readMapview(pileup[i], current_chr);
	      HMM hmm = new HMM(current_chr, space_between_duplications, space_between_deletions, average_CNV_length, homozygous_probability, sizes, carpet,		GenomeSegments.isNumbered(current_chr) ? numbered_segments : XY_segments, pileup_data);
	      optimal_solution.addAll(hmm.optimalSolution());
	      if(report_count_threshold < GenomeSegments.chr_length[current_chr.ordinal()])
		      hmm_map.put(current_chr, hmm);
	      }
	      


	
	    if(actual == null)
	      for(Iterator<CNV> iter = optimal_solution.iterator(); iter.hasNext(); )
	      System.out.println(iter.next().toString());
	    else
	      compare_CNV_lists(actual, optimal_solution);
    } 
      
    catch (Exception e) {
	    e.printStackTrace();
    }
  }

  static void compare_CNV_lists(List<CNV> actual, List<CNV> optimal_solution) {
	  ListIterator<CNV> actual_iter = actual.listIterator(),
			 optimal_iter = optimal_solution.listIterator();
	  CNV current_actual = actual_iter.hasNext() ? actual_iter.next() : null;
	  CNV current_optimal = optimal_iter.hasNext() ? optimal_iter.next() : null;
	  while(current_actual != null || current_optimal != null){
	    if(current_actual != null && current_optimal != null &&
			current_actual.overlap(current_optimal) > 0)
	    {
	      DividedCNV divided_actual = new DividedCNV(current_actual), divided_optimal = new DividedCNV(current_optimal);
	      for(int index = actual_iter.nextIndex(); index < actual.size() && actual.get(index).overlap(current_optimal) > 0; index++)
		divided_actual.add(actual.get(index));
	      for(int index = optimal_iter.nextIndex(); index < optimal_solution.size() && optimal_solution.get(index).overlap(current_actual) > 0; index++)
		divided_optimal.add(optimal_solution.get(index));
	      float overlap = divided_optimal.overlap(divided_actual);
	      if(overlap > overlap_threshold){
	        System.out.println(divided_actual.toString() + " found " +
				    divided_optimal.toString() + " overlap " +
				    Math.round(overlap*100) + '%');
		for(int i=0; i < divided_actual.size(); i++)
		    current_actual = actual_iter.hasNext() ? actual_iter.next() : null;
		for(int i=0; i < divided_optimal.size(); i++)
		    current_optimal = optimal_iter.hasNext() ? optimal_iter.next() : null;
		continue;
	      }
	    }
	    if(current_actual == null ||
			current_optimal != null && CNV.ByStart.compare(current_optimal, current_actual) > 0)
	    {
		System.out.println("false positive: " + report(current_optimal));
		current_optimal = optimal_iter.hasNext() ? optimal_iter.next() : null;
	    } else {
		System.out.println("undetected: " + report(current_actual));
		current_actual = actual_iter.hasNext() ? actual_iter.next() : null;
	    }
	  }
  }


// parse mapping file ;  
// TODO:  
  static int[] readMapview(File[] mapview, GenomeSegments.chrID chr) throws Exception {
    if (mapview == null)
      return null;
        
    int[] result = new int[GenomeSegments.chr_length[chr.ordinal()]+1];
    Arrays.fill(result, 0);
      
    for(int i=0; i < mapview.length; i++) {
      BufferedReader mapviewFile = new BufferedReader(new FileReader(mapview[i]));
        
      int current_pos = 1;
      for (String next = mapviewFile.readLine(); next != null; next = mapviewFile.readLine()) {
        String[] cols = next.split("\\s");
       // cols[0] : read name
       // cols[1] : ref name
       // cols[2] : start position
       // cols[3] : direction
       // cols[4] : distance 
       // cols[5] : maq flag; 18 mean good pairs. Should be pre-processed
        if(GenomeSegments.chrID.valueOf(GenomeSegments.chrID.class, cols[1]) != chr)
          continue;
       
        int start = Integer.parseInt(cols[2]);
        result[start] += 1;
        result[0] = Math.max(result[0], result[start]);
       
      }
  
    }
    return result;
  }


   static int[] readPileup(File[] pileup, GenomeSegments.chrID chr, int read_length)
Ê Ê throws Exception
Ê Ê {
Ê Ê Ê Ê Ê if(pileup == null)
Ê Ê Ê Ê Ê Ê return null;
Ê Ê Ê Ê Ê int[] result = new int[GenomeSegments.chr_length[chr.ordinal()]+1];
Ê Ê Ê Ê Ê Arrays.fill(result, 0);
Ê Ê Ê Ê Ê for(int i = 0; i < pileup.length; i++){
Ê Ê Ê Ê Ê Ê BufferedReader pileupfile = new BufferedReader(new FileReader(pileup[i]));
Ê Ê Ê Ê Ê Ê int current_pos = 1;
Ê Ê Ê Ê Ê Ê for(String next = pileupfile.readLine(); next != null; next = pileupfile.readLine())
Ê Ê Ê Ê Ê Ê {
Ê Ê Ê Ê Ê Ê Ê Ê String[] tokens = next.split("\\t");
Ê Ê Ê Ê Ê Ê Ê Ê boolean compact_format = tokens.length < 5 ||
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê tokens[4].charAt(0) != '@';
Ê Ê Ê Ê Ê Ê Ê Ê if(compact_format)
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê tokens = next.split("[ :]+");
Ê Ê Ê Ê Ê Ê Ê Ê if(GenomeSegments.chrID.valueOf(GenomeSegments.chrID.class, tokens[0]) != chr)
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê continue;
Ê Ê Ê Ê Ê Ê Ê Ê int num_pos = 1;
Ê Ê Ê Ê Ê Ê Ê Ê int start = -1;
Ê Ê Ê Ê Ê Ê Ê Ê if(tokens[1].indexOf('-') > 0){
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê String[] pos = tokens[1].split("-");
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê start = Integer.parseInt(pos[0]);
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê num_pos = Integer.parseInt(pos[1]) - start + 1;
Ê Ê Ê Ê Ê Ê Ê Ê } else
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê start = Integer.parseInt(tokens[1]);
Ê Ê Ê Ê Ê Ê Ê Ê if(start != current_pos)
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê throw new Exception("bad pileup line: " + (current_pos-1) + " followed by " + start);
Ê Ê Ê Ê Ê Ê Ê Ê int section_coverage = 0;
Ê Ê Ê Ê Ê Ê Ê Ê if(compact_format){
Ê Ê Ê Ê Ê Ê Ê Ê Ê for(int j=2; j < tokens.length; j++)
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê if(Character.isDigit(tokens[j].charAt(0)))
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê section_coverage += Integer.parseInt(tokens[j].substring(0,tokens[j].length()-1));
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê else
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê section_coverage += tokens[j].length();
Ê Ê Ê Ê Ê Ê Ê Ê } else
Ê Ê Ê Ê Ê Ê Ê Ê Ê section_coverage = tokens[4].length()-1;
Ê Ê Ê Ê Ê Ê Ê Ê for(; current_pos < start + num_pos; current_pos++){
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê result[current_pos] += section_coverage;
Ê Ê Ê Ê Ê Ê Ê Ê Ê Ê result[0] = Math.max(result[0], result[current_pos]);
Ê Ê Ê Ê Ê Ê Ê Ê }
Ê Ê Ê Ê Ê Ê }
Ê Ê Ê Ê Ê }
Ê Ê Ê Ê Ê convertToStartRead(result, read_length);
Ê Ê Ê Ê Ê return result;
Ê Ê }

Ê Ê static void convertToStartRead(int[] pileup, int read_length)
Ê Ê {
Ê Ê Ê Ê Ê pileup[0] = 0;
Ê Ê Ê Ê Ê for(int previous_coverage = 0, current_pos = 1, current_coverage; current_pos < pileup.length; current_pos++, previous_coverage = current_coverage){
Ê Ê Ê Ê Ê Ê current_coverage = pileup[current_pos];
Ê Ê Ê Ê Ê Ê pileup[current_pos] -= previous_coverage;
Ê Ê Ê Ê Ê Ê if(current_pos > read_length)
Ê Ê Ê Ê Ê Ê Ê Ê pileup[current_pos] += pileup[current_pos-read_length];
Ê Ê Ê Ê Ê Ê pileup[0] = Math.max(pileup[0], pileup[current_pos]);
Ê Ê Ê Ê Ê }
Ê Ê }


  static String report(CNV cnv) {
  	String result = cnv.toString() + " size " + cnv.size;
	  return result;
  }
}

import java.io.*;
import java.util.*;

public class qualityStats {
   
    public static void main(String[] args) {
      try {
	final float desired_quality = (float).95;
	String[] results_files = (new File(args[0])).list();
	Arrays.sort(results_files);
	Map<String, String> covlist = new HashMap<String, String>();
	int start_thresholds = 1;
	if(!Character.isDigit(args[1].charAt(0))){
	    BufferedReader covfile = new BufferedReader(new FileReader(args[start_thresholds++]));
	    for(String next = covfile.readLine(); next != null; next = covfile.readLine())
	    {
		String[] tokens = next.split("\\s+");
		covlist.put(tokens[0], tokens[1]);
	    }
	}

	for(int i=start_thresholds; i==1 || i < args.length; i++){
	    int threshold = i < args.length ? Integer.parseInt(args[i]) : 1;
	    if(i > 1)
		System.out.println();
	    if(threshold > 1)
		System.out.println("Threshold: " + threshold);
	    List<Integer> TP = new Vector<Integer>(), FP = new Vector<Integer>(), FN = new Vector<Integer>();
	    String previous_cov = null;
	    for(int result = 0; result <= results_files.length; result++){
		String cov = result == results_files.length ? "" :
				results_files[result].substring(8,
					results_files[result].indexOf('.', 8));
		if(!cov.equals(previous_cov) && !TP.isEmpty()){
		    Collections.sort(TP, Collections.reverseOrder());
		    Collections.sort(FP, Collections.reverseOrder());
		    Collections.sort(FN, Collections.reverseOrder());
		    int sensitivity_size = FN.isEmpty() ? 0 : FN.get(0).intValue(),
			specificity_size = FP.isEmpty() ? 0 : FP.get(0).intValue();
		    int TP_size = Integer.MAX_VALUE, FP_size = Integer.MAX_VALUE, FN_size = Integer.MAX_VALUE;
		    boolean found_sensitivity = false, found_specificity = false;
		    for(ListIterator<Integer> TP_iter = TP.listIterator(),
						FP_iter = FP.listIterator(),
						FN_iter = FN.listIterator();
						(!found_sensitivity || !found_specificity) && TP_iter.hasNext(); )
		    {
			TP_size = TP_iter.next().intValue();
			while(FP_iter.hasNext() && FP_size >= TP_size)
			    FP_size = FP_iter.next().intValue();
			while(FN_iter.hasNext() && FN_size >= TP_size)
			    FN_size = FN_iter.next().intValue();
			// System.out.println("previous_cov " + previous_cov + " sizes TP " + TP_size + " FP " + FP_size + " FN " + FN_size + " counts TP " + TP_iter.nextIndex() + " FP " + FP_iter.nextIndex() + " FN " + FN_iter.nextIndex());
			if(!found_sensitivity && (float)TP_iter.nextIndex() / (TP_iter.nextIndex() + FN_iter.nextIndex() - 1) < desired_quality)
			{
			    sensitivity_size = TP_size;
			    found_sensitivity = true;
			}
			if(!found_specificity && (float)TP_iter.nextIndex() / (TP_iter.nextIndex() + FP_iter.nextIndex() - 1) < desired_quality)
			{
			    specificity_size = TP_size;
			    found_specificity = true;
			}
		    }
		    TP.clear();
		    FP.clear();
		    FN.clear();
		    String cov_display = covlist.get(previous_cov);
		    if(cov_display == null)
			cov_display = previous_cov;
		    System.out.println(cov_display + ' ' + sensitivity_size/1000. + ' ' + specificity_size/1000.);
		}
		previous_cov = cov;
		if(result < results_files.length){
		    BufferedReader infile = new BufferedReader(new FileReader(args[0] + File.separator + results_files[result]));
		    for(String line = infile.readLine(); line != null; line = infile.readLine())
		    {
			String[] tokens = line.split("[-: \\]%,]+");
			if(tokens[0].equals("false") || tokens[0].equals("undetected")){
			    (tokens[0].equals("false") ? FP : FN) . add(
				new Integer(tokens[Arrays.asList(tokens).indexOf("size")+1]));
			    continue;
			}
			Collection<Integer> orig = new Vector<Integer>(),
					    detected = new Vector<Integer>();
			int index = -1;
			boolean reading_orig = true;
			while(!tokens[++index].equals("overlap"))
			if(tokens[index].startsWith("chr"))
			    (reading_orig ? orig : detected) . add
				(Integer.parseInt(tokens[index+2]) -
				 Integer.parseInt(tokens[index+1]) + 1);
			else if(tokens[index].equals("found"))
			    reading_orig = false;
			int overlap = Integer.parseInt(tokens[index+1]);
			if(overlap >= threshold){
			    int sum = 0;
			    for(Iterator<Integer> iter = orig.iterator(); iter.hasNext(); )
				sum += iter.next().intValue();
			    TP.add(sum);
			} else {
			    FN.addAll(orig);
			    FP.addAll(detected);
			}
		    }
		    infile.close();
		}
	    }
	}
      } catch (Exception e) {
	e.printStackTrace();
      }
    }
}

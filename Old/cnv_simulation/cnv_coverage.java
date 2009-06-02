package cnv_simulation;
import java.io.*;
import java.util.*;

public class cnv_coverage {

  public static void main(String[] args) {
    try {
      BufferedReader infile = new BufferedReader(new FileReader(args[1]));
      List<CNV> CNVs = new Vector<CNV>();
      for(String next = infile.readLine(); next != null; next = infile.readLine())
	CNVs.add(new CNV(next));
      Collections.sort(CNVs, Collections.reverseOrder(CNV.ByStart));
      infile = new BufferedReader(new FileReader(args[0]));
      String next = infile.readLine();
      for(Iterator<CNV> iter = CNVs.iterator(); iter.hasNext(); ){
	CNV cnv = iter.next();
	System.out.println(cnv.toString() + ':');
	while(next != null && !cnv.overlapsLine(next))
	    next = infile.readLine();
	compute_distribution distr = new compute_distribution();
	while(next != null && cnv.overlapsLine(next)){
	    distr.add_to_counts(next);
	    next = infile.readLine();
	}
	distr.printCounts();
	System.out.println();
      }
    } catch (Throwable e) {
	e.printStackTrace();
    }
  }

}

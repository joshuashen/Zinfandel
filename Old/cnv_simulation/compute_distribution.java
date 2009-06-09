package cnv_simulation;
import java.io.*;
import java.util.*;

public class compute_distribution {

  static String[] tokens = null;
  public int[] counts = new int[1000];

  public compute_distribution(){
      Arrays.fill(counts, 0);
  }

  public void add_to_counts(String next){
	tokens = next.split("[ :]+");
	int pos_num = 1;
	if(tokens[1].indexOf('-') > 0){
	    String[] pos = tokens[1].split("-");
	    pos_num = Integer.parseInt(pos[1]) - Integer.parseInt(pos[0]) + 1;
	}
	int coverage = 0;
	for(int i=2; i < tokens.length; i++)
	if(Character.isDigit(tokens[i].charAt(0)))
	    coverage += Integer.parseInt(tokens[i].substring(0,tokens[i].length()-1));
	else
	    coverage += tokens[i].length();
	counts[coverage] += pos_num;
  }

  public void printCounts() {
      for(int i=0; i < counts.length; i++)
      if(counts[i] > 0)
	System.out.println(String.valueOf(i) + ' ' + counts[i]);
  }

  public static void main(String[] args) {
    try {
      BufferedReader infile = new BufferedReader(new FileReader(args[0]));
      compute_distribution distr = new compute_distribution();
      for(String next = infile.readLine(); next != null; next = infile.readLine())
	distr.add_to_counts(next);
      distr.printCounts();
    } catch (Throwable e) {
	System.err.println(Arrays.toString(tokens));
	e.printStackTrace();
    }
  }

}

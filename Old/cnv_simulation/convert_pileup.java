import java.io.*;
import java.util.*;

public class convert_pileup {
 private static String chr_id = null, block_start = null,
			block_end = null, pos_info = null,
			last_block_info = null;

  public static void main(String[] args) {
    try {
      BufferedReader infile = new BufferedReader(new InputStreamReader(System.in));
      for(String next = infile.readLine(); next != null; next = infile.readLine())
      {
	String[] tokens = next.split("\\t");
	pos_info = decode(tokens[4]);
	if(!pos_info.equals(last_block_info) || !tokens[0].equals(chr_id)){
	  print_info();
	  chr_id = tokens[0];
	  block_start = tokens[1];
	  last_block_info = pos_info;
	}
	block_end = tokens[1];
      }
      print_info();
    } catch (IOException e) {
	e.printStackTrace();
    }
  }

  private static void print_info() {
    if(chr_id == null)
	return;
    System.out.print(chr_id + ':' + block_start);
    if(!block_end.equals(block_start))
	System.out.print("-" + block_end);
    System.out.println(" :" + last_block_info);
  }

  static private String chars = ",.ACTGNactgn";
  static private int[] counts = new int['u'];
  private static String decode(String maq_string)
  {
    int num_hits = maq_string.length()-1;
    if(num_hits == 0)
	return "";
    Arrays.fill(counts, 0);
    String result = " ";
    for(int i=1; i < maq_string.length(); i++)
	counts[maq_string.charAt(i)]++;
    char ch;
    boolean multi = false;
    for(int i=0; i < chars.length(); i++)
    if(counts[ch = chars.charAt(i)] > 0)
    {
      if(result.length() > 1 && (multi || counts[ch] > 1))
	result += ' ';
      if(multi = counts[ch] > 1)
	result += String.valueOf(counts[ch]);
      result += ch;
    }
    return result;
  }
}

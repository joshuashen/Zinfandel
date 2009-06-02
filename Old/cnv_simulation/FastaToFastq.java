package cnv_simulation;
import java.io.*;
import java.util.*;

public class FastaToFastq {
    public static void main(String[] args) {
      try {
	if(args.length != 2){
	    System.err.println("usage: FastaToFastq param-file fasta-file");
	    System.exit(0);
	}

	char quality = 53;
	String qual_string = "";
	BufferedReader paramfile = new BufferedReader(new FileReader(args[0]));
	for(String paramline = paramfile.readLine(); paramline != null; paramline = paramfile.readLine())
	{
	    String[] tokens = paramline.split("\\s+");
	    if(tokens[0].equalsIgnoreCase("quality"))
		quality = (char)(33 + Integer.parseInt(tokens[1]));
	}

	BufferedReader results_file = new BufferedReader(new FileReader(args[1]));
	for(String read_name = results_file.readLine(); read_name != null; read_name = results_file.readLine())
	{
	    System.out.println("@" + read_name.substring(1));
	    String seq = results_file.readLine();
	    System.out.println(seq);
	    System.out.println("+");
	    if(qual_string.length() > seq.length())
		qual_string.substring(0, seq.length());
	    while(qual_string.length() < seq.length())
		qual_string += quality;
	    System.out.println(qual_string);
	}
      } catch (IOException e) {
	e.printStackTrace();
      }
    }
}

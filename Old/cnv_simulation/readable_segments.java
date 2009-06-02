import java.io.*;

public class readable_segments {

    public static void main(String[] args) {
      try {
	for(int chrn = 1; chrn <= 24; chrn++){
	    String chr = String.valueOf(chrn);
	    if(chrn == 23)
		chr = "X";
	    else if(chrn == 24)
		chr = "Y";
	    BufferedReader infile = new BufferedReader(new FileReader("chr" + chr + ".fa"));
	    infile.readLine();
	    String next = infile.readLine();
	    int start = 1;
	    for(int segnum = 1; next != null; segnum++){
		while(next != null && next.matches("N*")){
		    start += next.length();
		    next = infile.readLine();
		}
		if(next == null)
		    break;
		while(next.charAt(0) == 'N'){
		    next = next.substring(1);
		    start++;
		}
		PrintStream segmentfile = new PrintStream(new FileOutputStream("chr" + chr + ".segment" + segnum));
		String header = ">c" + chr + ".s" + segnum + ' ' + start;
		while(header.length() < 19)
		    header += ' ';
		segmentfile.println(header);
		String segment = next;
		start += next.length();
		for(next = infile.readLine(); next != null && next.indexOf('N') < 0; next = infile.readLine()){
		  segment += next;
		  start += next.length();
		  if(segment.length() > 50){
		    segmentfile.println(segment.substring(0,50));
		    segment  = segment.substring(50);
		  }
		}
		if(next != null){
		  segment += next.substring(0, next.indexOf('N'));
		  next = infile.readLine();
		}
		while(segment.length() > 50){
			segmentfile.println(segment.substring(0,50));
			segment  = segment.substring(50);
		}
		segmentfile.println(segment);
		segmentfile.close();
	    }
	    infile.close();
	}
      } catch (IOException e) {
	e.printStackTrace();
      }
    }
}

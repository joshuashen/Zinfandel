package cnv_simulation;
import java.io.*;
import java.util.*;

public class CNV implements Cloneable {
    public static enum CNV_type { DUP, DEL, INV, UNKNOWN};

    public CNV_type type;
    public GenomeSegments.chrID chromosome;
    int relative_start = -1;
    public int absolute_start = -1;
    public int size;
    public boolean heterozygous = false;

    public CNV(String CNVline)
    {
      String[] tokens = CNVline.split("\\t");
      if(tokens.length == 1){
	tokens = CNVline.split("[ :\\-]");
	boolean old_format = tokens[3].equals("start");
	for(CNV_type t : CNV_type.values())
	  if(t.toString().equals(tokens[0]))
	    type = t;
	for(GenomeSegments.chrID c : GenomeSegments.chrID.values())
	  if(c.toString().equals(tokens[old_format ? 4 : 1]))
	    chromosome = c;
	absolute_start = Integer.parseInt(tokens[old_format ? 5 : 2]);
	heterozygous = !tokens[tokens.length-1].equals("homozygous");
	if(old_format)
	    size = Integer.parseInt(tokens[2]);
	else
	    size = Integer.parseInt(tokens[3]) - absolute_start + 1;
      } else {
	chromosome = GenomeSegments.getchrID(tokens[2].substring(3));
	absolute_start = Integer.parseInt(tokens[3]);
	size = Integer.parseInt(tokens[4]) - absolute_start + 1;
	heterozygous = GenomeSegments.isNumbered(chromosome);
	type = tokens[5].startsWith("Invers") ?  CNV_type.INV :
		(Simulation.generator.nextFloat() < .5 ? CNV_type.DUP : CNV_type.DEL);
      }
    }

    public CNV(CNV_type type, GenomeSegments.chrID chromosome, int segment_length, int size, boolean het)
    { this.type = type; this.chromosome = chromosome;
	this.size = size;
	heterozygous = het && GenomeSegments.isNumbered(chromosome);
	relative_start = Simulation.generator.nextInt(segment_length - size);
    }

    String apply(String segment)
    {
      try {
	if(segment == null)
	    return null;
	String result = segment.substring(0, relative_start);
	String duplicated_area = segment.substring(relative_start, relative_start+size);
	if(type == CNV_type.DUP)
	    result += duplicated_area + duplicated_area;
	else if(type == CNV_type.INV)
	    result += sequence_read.complement(duplicated_area);
	result += segment.substring(relative_start+size, segment.length());
	return result;
      } catch (RuntimeException e) {
	System.err.println("While applying " + toString() + " relative start " + relative_start + " to segment length " + segment.length());
	throw e;
      }
    }

    public String toString()
    {
	String result = type.toString() + ' ' + chromosome + ':' + absolute_start + '-' + (absolute_start+size-1);
	if(!heterozygous && GenomeSegments.isNumbered(chromosome))
	    result += " homozygous";
	return result;
    }

    public float overlap(CNV other)
    {
	if(chromosome != other.chromosome ||
		type != other.type && type != CNV_type.UNKNOWN &&
		other.type != CNV_type.UNKNOWN)
	    return 0;
	int start_overlap = Math.max(absolute_start, other.absolute_start);
	int end_overlap = Math.min(absolute_start+size, other.absolute_start+other.size);
	if(end_overlap < start_overlap)
	    return 0;
	return (float)(end_overlap - start_overlap) / Math.max(size, other.size);
    }

    public boolean equals(Object o)
    {
	if(!(o instanceof CNV))
	    return false;

	CNV other = (CNV)o;
	boolean result = chromosome == other.chromosome && size == other.size &&
		(relative_start == -1 || other.relative_start == -1 || relative_start == other.relative_start) &&
		(absolute_start == -1 || other.absolute_start == -1 || absolute_start == other.absolute_start);
	return result;
    }

    public CNV clone() {
	try {
	    return (CNV)super.clone();
	} catch (CloneNotSupportedException e) {
	    return null;
	}
    }

    GenomeSegments.Segment find_segment()
    {
	GenomeSegments group = GenomeSegments.isNumbered(chromosome)
			? Simulation.numbered : Simulation.XY;
	int result = Arrays.binarySearch(group.segments,
		new GenomeSegments.Segment(chromosome, absolute_start));
	if(result < 0)
	   result = - result - 2;
	if(result < 0 ||
	   group.segments[result].start + group.segments[result].length < absolute_start+size)
	   return null;
	relative_start = absolute_start - group.segments[result].start;
	return group.segments[result];
    }

    public boolean overlapsLine(String line){
	String[] tokens = line.split("[ :]+");
	if(!chromosome.name().equals(tokens[0]))
	    return false;
	int lineStart, lineEnd, dash = tokens[1].indexOf('-');
	if(dash > 0){
	    lineStart = Integer.parseInt(tokens[1].substring(0,dash));
	    lineEnd = Integer.parseInt(tokens[1].substring(dash+1));
	} else
	    lineStart = lineEnd = Integer.parseInt(tokens[1]);
	return lineStart < absolute_start+size &&
	       lineEnd >= absolute_start;
    }

    public static class BySizeComparator implements Comparator<CNV> {
      ByStartComparator alternate = new ByStartComparator();
      public int compare(CNV c1, CNV c2)
      {
	int result = c2.size - c1.size;
	return result == 0 ? alternate.compare(c1,c2) : result;
      }
    }

    public static class ByStartComparator implements Comparator<CNV> {
      public int compare(CNV c1, CNV c2)
      {
	int result = c1.chromosome.compareTo(c2.chromosome);
	if(result == 0)
	    result = c1.absolute_start >= 0 && c2.absolute_start >= 0 ?
			c2.absolute_start - c1.absolute_start :
			c2.relative_start - c1.relative_start;
	if(result == 0)
	    result = c2.size - c1.size;
	return result;
      }
    }

    public static BySizeComparator BySize = new BySizeComparator();
    public static ByStartComparator ByStart = new ByStartComparator();
}

public class Chromosome {
    String name;        //Chromosome Name
    int size;           //Seqence Length of the Chromosome
    byte[] coverage;     //Depth-Coverage: Number of reads that begin at the Base-Pair
    int insertSize;     //Ignored for now

    public Chromosome(String n, int s){
        name = n;
        size = s;
        coverage = new byte[size];
    }

    //Increments the depth coverage for the base-pair position
    public void incrementCoverage(int pos){
        coverage[pos-1]++;
    }
}

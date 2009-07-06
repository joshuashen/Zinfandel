public class Chromosome {
    String name;        //Chromosome Name
    int size;           //Seqence Length of the Chromosome
    byte[] coverage;     //Depth-Coverage: Number of reads that begin at the Base-Pair
    int[] distances;    //Paired Ends Distance
    int insertSize;     //Ignored for now

    public Chromosome(String n, int s){
        name = n;
        size = s;
        coverage = new byte[size];
        distances = new int[size];

        //initialize all distance to -1, signifies no read at this position i if distance[i] remains -1
        for (int i = 0; i<size; i++){
            distances[i] = -1;
        }
    }

    //Increments the depth coverage for the base-pair position
    public void incrementCoverage(int pos){
        coverage[pos-1]++;
    }

    //Sets Paired Distance for base-pair position
    public void setDistance(int pos, int distance){
        distances[pos-1] = distance;
    }
}
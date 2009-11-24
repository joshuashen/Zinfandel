package cnv_hmm;
public class Chromosome {
    String name; //Chromosome Name
    int size; //Seqence Length of the Chromosome
    int maxDisPerPosition;
    byte[] coverage; //Depth-Coverage: Number of reads that begin at the Base-Pair
    int[][] distances; //Paired Ends Distance



    int[] flags; //flags(used for explicit breakpoint modelling)
    int insertSize; //Ignored for now


    public Chromosome(String n, int s, int d){
        name = n;
        size = s;
        maxDisPerPosition = d;
        coverage = new byte[size];
        distances = new int[size][maxDisPerPosition];


        //flags = new int[size];

        //initialize all distance to -1, signifies no read at this position i if distance[i] remains -1
        for (int i = 0; i<size; i++){
            for (int j = 0; j<maxDisPerPosition; j++){
                distances[i][j] = 0;
            }
        }
    }

    //Increments the depth coverage for the base-pair position
    public void incrementCoverage(int pos){
if (pos > size -1 || pos < 0) {
return;
}
        coverage[pos-1]++;
    }

    //Sets Paired Distance for base-pair position
    // if the position is already set, slip the distance to the next position
    public void setDistance(int pos, int distance){
        if (pos > size -1 || pos < 0) {
            return;
        }
        //If number of reads at position exceeds maximum set, then read distance is ignored
        for (int i=0; i<maxDisPerPosition; i++){
            if (distances[pos-1][i] == 0){
                distances[pos-1][i] = distance;
                break;
            }
        }
    }
/*
//Sets flag for base-pair position
public void setFlag(int pos, int flag){
flags[pos-1] = flag;
}
*/
}
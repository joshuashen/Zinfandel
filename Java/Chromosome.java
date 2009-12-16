
import java.util.LinkedList;

public class Chromosome {
    String name; //Chromosome Name
    int size; //Seqence Length of the Chromosome

    byte[] coverage; //Depth-Coverage: Number of reads that begin at the Base-Pair
    LinkedList<Integer>[] distances; //Paired Ends Distance

    int[] flags; //flags(used for explicit breakpoint modelling)
    int insertSize; //Ignored for now


    public Chromosome(String n, int s){
        name = n;
        size = s;

        coverage = new byte[size];
        distances = new LinkedList[size];
        for (int i=0; i<size; i++){
            distances[i] = new LinkedList<Integer>();
        }
        //flags = new int[size];
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
        distances[pos-1].add(distance);
    }
/*
//Sets flag for base-pair position
public void setFlag(int pos, int flag){
flags[pos-1] = flag;
}
*/
}
package cnv_hmm;
import java.util.ArrayList;

public class EmissionDistanceMatrix {
    double[][] emissionDistanceMatrix;

    public EmissionDistanceMatrix(ArrayList<State> states, String genome, double avgDistance,
                                  double standardDevDistance, int maxDistance){
        int numStates = states.size();
        //Diploid Case
        if (genome.equalsIgnoreCase("d")){
            emissionDistanceMatrix = new double[numStates][maxDistance];

            //Set emission distance maxtrix for normal state
            for (int i = 0; i<maxDistance; i++){
                 emissionDistanceMatrix[0][i] = logNormal(i, avgDistance, standardDevDistance);
                 emissionDistanceMatrix[1][i] = logNormal(i, avgDistance, standardDevDistance);
                 emissionDistanceMatrix[2][i] = logNormal(i, avgDistance, standardDevDistance);
                 emissionDistanceMatrix[3][i] = logNormal(i, avgDistance, standardDevDistance);
                 emissionDistanceMatrix[4][i] = logNormal(i, avgDistance, standardDevDistance);
                 //Breakpoint State
                 emissionDistanceMatrix[numStates-1][i] = logNormal(i, avgDistance, standardDevDistance);
            }

            for (int i = 5; i<numStates-1; i++){
                for (int j = 0; j<maxDistance; j++){
                    if (states.get(i).isInitialGridState || states.get(i).isFinalGridState){                       
                        emissionDistanceMatrix[i][j] = logNormal(j, states.get(i).delSize + avgDistance, standardDevDistance);
                    }
                    else{
                        emissionDistanceMatrix[i][j] = logNormal(j, avgDistance, standardDevDistance);
                    }
                }
            }
        }
        //Haploid Case
        else{
            emissionDistanceMatrix = new double[numStates][maxDistance];

            //Set emission distance maxtrix for normal state, set to negative infinity for other states
            for (int i = 0; i<maxDistance; i++){
                emissionDistanceMatrix[0][i] = logNormal(i, avgDistance, standardDevDistance);
                emissionDistanceMatrix[1][i] = logNormal(i, avgDistance, standardDevDistance);
                emissionDistanceMatrix[2][i] = logNormal(i, avgDistance, standardDevDistance);
                emissionDistanceMatrix[numStates-1][i] = logNormal(i, avgDistance, standardDevDistance);
            }

            for (int i = 3; i<numStates-1; i++){
                for (int j = 0; j<maxDistance; j++){
                    if (states.get(i).isInitialGridState || states.get(i).isFinalGridState){                        
                        emissionDistanceMatrix[i][j] = logNormal(j, states.get(i).delSize + avgDistance, standardDevDistance);
                    }
                    else{
                        emissionDistanceMatrix[i][j] = logNormal(j, avgDistance, standardDevDistance);
                    }
                }
            }
        }
        System.out.println("Emission Matrix Created");        
    }
    
    //log(Normal, k, average, std deviation) = exp^(-(x-avg)^2) / 2*dev^2 - log(dev * sqrt(2pi))
    public double logNormal(int dis, double avg, double dev){
        double p = (-1) * ((double)dis - avg) * ((double)dis - avg);
        p = p/(2 * dev * dev);
        p -= Math.log(dev * Math.sqrt(2*Math.PI));
        return p;
    }
}

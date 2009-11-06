
import java.util.ArrayList;

public class EmissionDistanceMatrix {
    double[][] emissionDistanceMatrixPosStrand;
    double[][] emissionDistanceMatrixNegStrand;

    double minEmissionfromDistance;

    public EmissionDistanceMatrix(ArrayList<State> states, String genome, double avgDistance,
                                  double standardDevDistance, int maxDistance, double minEmission){
        int numStates = states.size();
        minEmissionfromDistance = minEmission;
        //haploid Case
emissionDistanceMatrixPosStrand = new double[numStates][maxDistance+1];
emissionDistanceMatrixNegStrand = new double[numStates][maxDistance+1];

        if (genome.equalsIgnoreCase("h") || genome.equalsIgnoreCase("d")){
            for (int dist = 0; dist<= maxDistance; dist++){
for (int state = 0; state < numStates; state++ ) {
// initial grid state: positive strand distance is larger; negative strand avg
if (states.get(state) instanceof FiveFlankingState) {
emissionDistanceMatrixPosStrand[state][dist] = logNormal(dist, ((GridState)states.get(state)).delSize + avgDistance, standardDevDistance);
emissionDistanceMatrixNegStrand[state][dist] = logNormal(dist, avgDistance, standardDevDistance);

}

// final grid state: negative strand distance is larger; positive strand avg
else if (states.get(state) instanceof ThreeFlankingState) {
emissionDistanceMatrixPosStrand[state][dist] = logNormal(dist, avgDistance, standardDevDistance);
emissionDistanceMatrixNegStrand[state][dist] = logNormal(dist, ((GridState)states.get(state)).delSize + avgDistance, standardDevDistance);
}
else {
emissionDistanceMatrixPosStrand[state][dist] = logNormal(dist, avgDistance, standardDevDistance);
emissionDistanceMatrixNegStrand[state][dist] = logNormal(dist, avgDistance, standardDevDistance);
}
}
}
        }
        //diploid Case, heterozygous
// else{
// for (int dist = 0; dist<= maxDistance; dist++){
// for (int state = 0; state < numStates; state++ ) {
// initial grid state: positive strand distance is larger; negative strand avg
// if (states.get(state) instanceof FiveFlankingState) {

// }
// }

// }
//}
        System.out.println("Emission Matrix Created");
    }

    //log(Normal, k, average, std deviation) = exp^(-(x-avg)^2) / 2*dev^2 - log(dev * sqrt(2pi))
    public double logNormal(int dis, double avg, double dev){
        double p = (-1) * ((double)dis - avg) * ((double)dis - avg);

        p = p/(2 * dev * dev);
        p -= Math.log(dev * Math.sqrt(2*Math.PI));
        if (p < this.minEmissionfromDistance) {
         p = this.minEmissionfromDistance;
}
        return p;
    }

    //log(Normal, k, average, std deviation) = exp^(-(x-avg)^2) / 2*dev^2 - log(2 * dev * sqrt(2pi))
    public double halfLogNormal(int dis, double avg, double dev){
        double p = (-1) * ((double)dis - avg) * ((double)dis - avg);
        p = p/(2 * dev * dev);
        p -= Math.log(2* dev * Math.sqrt(2*Math.PI));
        if (p < this.minEmissionfromDistance) {
         p = this.minEmissionfromDistance;
}
        return p;
    }
}
import java.util.ArrayList;

//Create Emission Matrix
//emission = Pr(depth-coverage | copy-number) * Pr(distance | state)
//add a uniform component unif:
//Pr(distance | CNV neutral) = unif + (1-unif) * Norm(distance | CNV-neutral)

public class EmissionCoverageMatrix {
    double[][] emissionCoverageMatrix;

    //Necessary Parameters: Arraylist of states, genome(haploid/diploid), average coverage, maximum coverage
    public EmissionCoverageMatrix(ArrayList<State> states, String genome, double avgCov, int maxCoverage){
        double lambda = avgCov;
        int numStates = states.size();
        double del1; double del2; double dup1; double dup2;
        //Diploid Case
        if (genome.equalsIgnoreCase("d")){
            del1 = lambda/2; del2 = lambda/100;
            dup1 = lambda * 1.5; dup2 = lambda * 2.0;
            emissionCoverageMatrix = new double[numStates][maxCoverage];

            for (int i = 0; i<maxCoverage; i++){
                emissionCoverageMatrix[0][i] = logPoisson(i, lambda);
                emissionCoverageMatrix[1][i] = logPoisson(i, del1);
                emissionCoverageMatrix[2][i] = logPoisson(i, del2);
                emissionCoverageMatrix[3][i] = logPoisson(i, dup1);
                emissionCoverageMatrix[4][i] = logPoisson(i, dup2);
                //Breakpoint States
                //emissionCoverageMatrix[numStates-1][i] = logPoisson(i, lambda);
            }
            //Set all emission coverage values for grid states- treat the same as Deletion1
            for (int i = 5; i<numStates-1; i++){
                for (int j = 0; j<maxCoverage; j++){
                    emissionCoverageMatrix[i][j] = logPoisson(j, del1);
                }
            }
        }
        //Haploid Case
        else{
            del1 = lambda/100; dup1 = lambda * 2;
            del2 = del1; dup2 = dup1;
            emissionCoverageMatrix = new double[numStates][maxCoverage];

            for (int i = 0; i<maxCoverage; i++){
                emissionCoverageMatrix[0][i] = logPoisson(i, lambda);
                emissionCoverageMatrix[1][i] = logPoisson(i, del1);
                emissionCoverageMatrix[2][i] = logPoisson(i, dup1);
                //Breakpoint state
                //emissionCoverageMatrix[numStates-1][i] = logPoisson(i, dup1);
            }
            //Set all emission coverage values for grid states- treat the same as Deletion1
            for (int i = 3; i<numStates-1; i++){
                for (int j = 0; j<maxCoverage; j++){
                    emissionCoverageMatrix[i][j] = logPoisson(j, del1);
                }
            }
        }
        System.out.println("Emission Coverage Matrix Created");
    }

    //log(Poisson(k,lambda)) = klog(lambda) - lambda - SUM(j)_1..k
    public double logPoisson(double k, double l){
        double p = k * Math.log(l) - l;
        for (int i = 1; i<=k; i++){
            p += -1.0 * Math.log(i);
        }
        return p;
    }
}

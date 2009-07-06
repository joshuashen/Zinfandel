package cnv_hmm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.ArrayList;

public class Cnv_Hmm {
    ArrayList<String> states = new ArrayList<String>();    //States for either Diploid or Haploid

    /*  Diploid: 0-Normal 1-del1 2-del2 3-dup1 4-dup2
     *  Haploid: 0-Normal 1-del1 2-dup1
     */
    double[][] transitionMatrix;    //Transition Probabilities
    double[][] emissionCoverageMatrix;      //Emission Coverage Probabilities
    double[][] emissionDistanceMatrix;       //Emission Distance Proabbilities
    
    //Default Values: Potentially overidden by parameter file
    int maxEmission = 100;  //Maximum emission, if larger, emission is set to 100
    int maxDistance = 1000; //Maximum paired distance, if larger, set to 1000
    double avgCov = 1.0;    //Average Depth Coverage
    String genome = "d";    //Default genome is Diploid
    int numCNVs = 3000;     //3000 CNVs assumed
    double avgCNVSize = 10000.0;    //Assumed size of CNVs
    double readSize = 35.0;               //Read Size
    double depthCov = 1.0;          //Depth Coverage
    double[] initProb;              //Initial Probabilities for each State
    double chrLength = 3000000000.0;
    int avgDistance;            //Paired distance(estimation)
    int standardDevDistance;    //Standard deviation for paired distance(estimation)
    int maxDeletionSize = 1000;    //Maximum Deletion Size
    int numGridStates;
    double factor = 1.5;
    int StatesPerDelSize = 10;
    ArrayList<Integer> DelSizes = new ArrayList<Integer>();

    //Default Constructor without Parameter File
    public Cnv_Hmm(){
        avgCov = depthCov / readSize;
        int count = 0;
        double value = 100;
        DelSizes.add((int)value);
        while(value < maxDeletionSize){
            value = value * factor;
            DelSizes.add((int)value);
            count++;
        }
        numGridStates = count;
    }
    //Parameter File Provided
    public Cnv_Hmm(File params){
        int count = 0;
        double value = 100;
        DelSizes.add((int)value);
        while(value < maxDeletionSize){
            value = value * factor;
            DelSizes.add((int)value);
            count++;
        }
        numGridStates = count;
        readParameterFile(params);
        avgCov = depthCov / readSize;
    }

    public void initializeHaploidStates(){
        states.add(0, "normal"); states.add(1, "del1"); states.add(2, "dup1");
        for (int i = 0; i<numGridStates; i++){
            for (int j = 0;j<StatesPerDelSize; j++){
                states.add("delSize" + DelSizes.get(i).toString() + "-" + j);
            }
        }
    }

    public void initializeDiploidStates(){
        states.add(0, "normal"); states.add(1, "del1"); states.add(2, "del2");
        states.add(3, "dup1"); states.add(4, "dup2");
        for (int i = 0; i<numGridStates; i++){
            for (int j = 0;j<StatesPerDelSize; j++){
                states.add("delSize" + DelSizes.get(i).toString() + "-" + j);
            }
        }
    }

    public void readParameterFile(File params){
        try{
           FileReader reader = new FileReader(params);
           BufferedReader in = new BufferedReader(reader);
           System.out.println("Loading Parameter File.");
           String line;    //stores line
                while((line = in.readLine()) != null){
                    StringTokenizer st = new StringTokenizer(line, " ");
                    String parameterName = st.nextToken();
                    String parameterData = st.nextToken();
                    if (parameterName.equalsIgnoreCase("genome")){
                        genome = parameterData;
                        //Haploid
                        if (genome.equalsIgnoreCase("h")){
                            initializeHaploidStates();
                        }
                        //Diploid
                        else{
                            initializeDiploidStates();
                        }
                    }
                    //Set other parameters
                    else if (parameterName.equalsIgnoreCase("numCNVs")){
                        numCNVs = Integer.valueOf(parameterData);
                    }
                    else if (parameterName.equalsIgnoreCase("avgCNVSize")){
                        avgCNVSize = Double.valueOf(parameterData);
                    }
                    else if (parameterName.equalsIgnoreCase("depthCov")){
                        depthCov = Double.valueOf(parameterData);
                    }
                    else if (parameterName.equalsIgnoreCase("readSize")){
                        readSize = Integer.valueOf(parameterData);
                    }
                }
           System.out.println("Parameter File Processed");
           in.close();
        }
        catch(IOException e){
           e.printStackTrace();
        }
    }

    //Creates Transition Matrix and Sets Initial Probabilities for each State
    public void initializeTransitionMatrix(){
        int numStates = states.size();
        initProb = new double[numStates];

        //Create Transition Matrix and initialize all values to negative infinity
        transitionMatrix = new double[numStates][numStates];
        for (int i=0; i<numStates; i++){
            for (int j=0; j<numStates; j++){
                transitionMatrix[i][j] = Double.NEGATIVE_INFINITY;
            }
        }

        for (int i=0; i<numStates; i++){
            initProb[i] = Math.log(1.0/numStates);
            for (int j=0; j<numStates; j++){
                if (i == j){
                    if (states.get(i).equalsIgnoreCase("normal")){
                        transitionMatrix[i][j] = Math.log(1 - (numCNVs/chrLength));
                    }
                    else if (states.get(i).equalsIgnoreCase("del1") || states.get(i).equalsIgnoreCase("dup1")){
                        transitionMatrix[i][j] = Math.log(1 - (1.0/avgCNVSize));
                    }                    
                    else{
                        int index = (i-3)/10;
                        //Pr(Transition from grid state to itself)
                        transitionMatrix[i][j] = Math.log(1-(double)StatesPerDelSize/DelSizes.get(index));
                    }                     
                }
                else{
                    if (states.get(i).equalsIgnoreCase("normal") && (states.get(j).contains("-0") || 
                        states.get(j).equalsIgnoreCase("del1") || states.get(j).equalsIgnoreCase("dup1"))){
                        //Assumed equal probability of transition from normal to del1, dup1, each initial grid state
                        transitionMatrix[i][j] = Math.log((numCNVs/chrLength)/(3 + DelSizes.size() - 1.0));
                    }
                    else if ((states.get(i).equalsIgnoreCase("del1") || states.get(i).equalsIgnoreCase("dup1"))
                             && (!states.get(j).contains("-"))){
                        //Assumed equal probability of transition from del1/dup1 to normal/del1/dup1 grid states
                        transitionMatrix[i][j] = Math.log((1.0/avgCNVSize)/(3.0 - 1.0));
                    }
                    else if(states.get(i).contains("-9") && (states.get(j).equalsIgnoreCase("normal"))){
                        //Assumed equal probability of transition from final grid state back to normal/dup1/del1
                        int index = (i-3)/10;
                        transitionMatrix[i][j] = Math.log((double)StatesPerDelSize/DelSizes.get(index));
                    }
                    else if(states.get(i).contains("-") && states.get(j).contains("-")){
                        String firstState = states.get(i);
                        String secondState = states.get(j);
                        StringTokenizer token1 = new StringTokenizer(firstState, "-");
                        StringTokenizer token2 = new StringTokenizer(secondState, "-");
                        String size1 = token1.nextToken();
                        String size2 = token2.nextToken();
                        int grid1 = Integer.valueOf(token1.nextToken());
                        int grid2 = Integer.valueOf(token2.nextToken());
                        if (grid1+1==grid2 && size1.equalsIgnoreCase(size2)){
                            int index = (i-3)/10;
                            transitionMatrix[i][j] = Math.log((double)StatesPerDelSize/DelSizes.get(index));
                        }
                    }
                }
            }
        }
        System.out.println("breakpoint");
    }

    //Create Emission Matrix
    //emission = Pr(depth-coverage | copy-number) * Pr(distance | state)
    //add a uniform component unif:
    //Pr(distance | CNV neutral) = unif + (1-unif) * Norm(distance | CNV-neutral)
    public void initializeEmissionMatrix(int avg, int dev){
        double lambda = avgCov;
        avgDistance = avg;
        standardDevDistance = dev;
        int numStates = states.size();
        double del1; double del2; double dup1; double dup2;
        //Diploid Case
        if (genome.equalsIgnoreCase("d")){
            del1 = lambda/2; del2 = lambda/100;
            dup1 = lambda * 1.5; dup2 = lambda * 2.0;
            emissionCoverageMatrix = new double[numStates][maxEmission];
            emissionDistanceMatrix = new double[numStates][maxDistance];
            for (int i = 0; i<maxEmission; i++){
                emissionCoverageMatrix[0][i] = logPoisson(i, lambda);
                emissionCoverageMatrix[1][i] = logPoisson(i, del1);
                emissionCoverageMatrix[2][i] = logPoisson(i, del2);
                emissionCoverageMatrix[3][i] = logPoisson(i, dup1);
                emissionCoverageMatrix[4][i] = logPoisson(i, dup2);
            }
            //Set all emission coverage values for grid states
            for (int i = 5; i<numStates; i++){
                for (int j = 0; j<maxEmission; j++){
                    emissionCoverageMatrix[i][j] = logPoisson(j, del1);
                }
            }
            //Set emission distance maxtrix for normal state
            for (int i = 0; i<maxDistance; i++){
                 emissionDistanceMatrix[0][i] = logNormal(i, avgDistance, standardDevDistance);
                 emissionDistanceMatrix[1][i] = logNormal(i, avgDistance, standardDevDistance);
                 emissionDistanceMatrix[2][i] = logNormal(i, avgDistance, standardDevDistance);
            }

            for (int i = 5; i<numStates; i++){
                for (int j = 0; j<maxDistance; j++){
                    if (states.get(i).contains("-0") || states.get(i).contains("-9")){
                        int index = (i-5)/10;
                        emissionDistanceMatrix[i][j] = logNormal(j, DelSizes.get(index) + avgDistance, standardDevDistance);
                    }
                    else{
                        emissionDistanceMatrix[i][j] = logNormal(j, avgDistance, standardDevDistance);
                    }
                }
            }
        }
        //Haploid Case
        else{
            del1 = lambda/100; dup1 = lambda * 2;
            del2 = del1; dup2 = dup1;
            emissionCoverageMatrix = new double[numStates][maxEmission];
            emissionDistanceMatrix = new double[numStates][maxDistance];
            for (int i = 0; i<maxEmission; i++){
                emissionCoverageMatrix[0][i] = logPoisson(i, lambda);
                emissionCoverageMatrix[1][i] = logPoisson(i, del1);
                emissionCoverageMatrix[2][i] = logPoisson(i, dup1);
            }
            //Set all emission coverage values for grid states
            for (int i = 3; i<numStates; i++){
                for (int j = 0; j<maxEmission; j++){
                    emissionCoverageMatrix[i][j] = logPoisson(j, del1);
                }
            }
            //Set emission distance maxtrix for normal state, set to negative infinity for other states
            for (int i = 0; i<maxDistance; i++){
                emissionDistanceMatrix[0][i] = logNormal(i, avgDistance, standardDevDistance);
                emissionDistanceMatrix[1][i] = logNormal(i, avgDistance, standardDevDistance);
                emissionDistanceMatrix[2][i] = logNormal(i, avgDistance, standardDevDistance);
            }

            for (int i = 3; i<numStates; i++){
                for (int j = 0; j<maxDistance; j++){
                    if (states.get(i).contains("-0") || states.get(i).contains("-9")){
                        int index = (i-3)/10;
                        emissionDistanceMatrix[i][j] = logNormal(j, DelSizes.get(index) + avgDistance, standardDevDistance);
                    }
                    else{
                        emissionDistanceMatrix[i][j] = logNormal(j, avgDistance, standardDevDistance);
                    }
                }
            }
        }
        System.out.println("breakpoint");
    }
    
    //log(Poisson(k,lambda)) = klog(lambda) - lambda - SUM(j)_1..k
    public double logPoisson(double k, double l){
        double p = k * Math.log(l) - l;
        for (int i = 1; i<=k; i++){
            p += -1.0 * Math.log(i);
        }
        return p;
    }
    //log(Normal, k, average, std deviation) = exp^(-(x-avg)^2) / 2*dev^2 - log(dev * sqrt(2pi))
    public double logNormal(int dis, int avg, int dev){
        double p = (-1) * (dis - avg) * (dis - avg);
        p = p/(2 * dev * dev);
        p -= Math.log(dev * Math.sqrt(2*Math.PI));
        return p;
    }

    //Accessor for Emission Matrix
    public double Emission(byte state, byte cov, int distance){
        if (cov > maxEmission){
            cov = (byte) maxEmission;
        }
        if (distance > maxDistance){
            distance = maxDistance;
        }
        //Read at position
        if (distance != -1){
            return emissionCoverageMatrix[state][cov] + emissionDistanceMatrix[state][distance];
        }
        //No read at position
        else{
            return emissionCoverageMatrix[state][cov];
        }
    }

    //Viterbi Algorithm for Optimal Path    
    public void runViterbiAlgorithm(Chromosome chr, String name){
        int chrSize = chr.size; //Chromosome length
        int numStates = states.size();  //Number of States
        double[] delta;
        if (genome.equalsIgnoreCase("h")){
            delta = new double[3];
        }
        else{
            delta = new double[5];
        }
        byte[][] prevMap = new byte[chrSize + 1][numStates];    //Holds Optimal Previous state
        //Set Initial Probabilities and prevMap
        for (byte i = 0; i<numStates; i++){
            delta[i] = initProb[i] + Emission(i, chr.coverage[0], chr.distances[0]);
            prevMap[0][i] = i;
        }

        double[] prob = new double[numStates];

        for (int i = 1; i<chrSize; i++){
            for (int j = 0; j<numStates; j++){     //Current State
                for (int k = 0; k<numStates; k++){  //Previous State
                    prob[k] = delta[k] + transitionMatrix[k][j];
                }
                //Find Optimal Probability
                double max = Double.NEGATIVE_INFINITY;  //MUST use NEGATIVE_INFINITY
                byte maxState = 0;
                for (byte x = 0; x<numStates; x++){
                    if (max < prob[x]){
                        max = prob[x];
                        maxState = x;
                    }
                }                
                prevMap[i][j] = maxState;   //Record Optimal State
                delta[j] = prob[prevMap[i][j]] + Emission((byte)j,chr.coverage[i],chr.distances[i]);
            }
        }
        //Find Final Optimal State
        double max = Double.NEGATIVE_INFINITY;  //MUST use NEGATIVE_INFINITY
        byte lastBest = 0;
        for (byte x = 0; x<numStates; x++){
            if (max < prob[x]){
                max = prob[x];
                lastBest = x;
            }
        }
        prevMap[chr.size][0] = lastBest;

        int flag = 0;   //Record if currently in a CNV
        int[] bestPath = new int[chr.size]; //Optimal Path

        int s = 0;  //Start index
        int e = 0;  //End index

        //Maintains CNV locations and type
        ArrayList<int[]> cnvs = new ArrayList<int[]>();        

        lastBest = 0;
        for (int i = chr.size; i>0; i--){
            byte state = prevMap[i][lastBest];
            s = i;
            if (state != 0){
                if (flag == 0){
                    flag = 1;
                    e =  i;
                }
                else if(state != lastBest){
                    int[] cnv = new int[3];
                    cnv[0] = s;
                    cnv[1] = e;
                    cnv[2] = lastBest;
                    cnvs.add(cnv);
                    e = i;
                }
            }
            else if (flag != 0){
                int[] cnv = new int[3];
                cnv[0] = s;
                cnv[1] = e;
                cnv[2] = lastBest;
                cnvs.add(cnv);
                flag = 0;
                s = i;
                e = i;
            }            
            lastBest = state;
        }

        //Prints out result either to specified file or to standard output
        System.out.println("\nCNVs: chromosome name, type, start point, end point, cnv length\n");
        for (int i = cnvs.size()-1; i >= 0; i--){
            int [] cnv = cnvs.get(i);
            System.out.print(name+"\t");
            //Haploid Case
            if (genome.equalsIgnoreCase("h")){
                switch(cnv[2]){
                    case 1: System.out.print("del1"); break;
                    case 2: System.out.print("dup1"); break;
                }
            }
            //Diploid Case
            else{
                switch(cnv[2]){
                    case 1: System.out.print("del1"); break;
                    case 2: System.out.print("del2"); break;
                    case 3: System.out.print("dup1"); break;
                    case 4: System.out.print("dup2"); break;
                }
            }
            System.out.printf("\t%10d", cnv[0] + 1);  //Increment by 1 to get correct start location
            System.out.printf("\t%10d", cnv[1] + 1);  //Increment by 1 to get correct end location
            System.out.printf("\t%10d", cnv[1] - cnv[0]);   //Length of CNV
            System.out.println();            
        }

    }

}


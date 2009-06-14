import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.ArrayList;

public class Cnv_Hmm {
    String[] states;    //States for either Diploid or Haploid
    /*  Diploid: 0-Normal 1-del1 2-del2 3-dup1 4-dup2
     *  Haploid: 0-Normal 1-del1 2-dup1
     */
    double[][] transitionMatrix;    //Transition Probabilities
    double[][] emissionMatrix;      //Emission Probabilities
    
    //Default Values: Potentially overidden by parameter file
    int maxEmission = 100;  //Maximum emission, if larger, emission is set to 100
    double avgCov = 1.0;    //Average Depth Coverage
    String genome = "d";    //Default genome is Diploid
    int numCNVs = 3000;     //3000 CNVs assumed
    double avgCNVSize = 10000.0;    //Assumed size of CNVs
    double readSize = 35.0;               //Read Size
    double depthCov = 1.0;          //Depth Coverage
    double[] initProb;              //Initial Probabilities for each State
    double chrLength = 3000000000.0;

    //Default Constructor without Parameter File
    public Cnv_Hmm(){
        avgCov = depthCov / readSize;
    }

    //Parameter File Provided
    public Cnv_Hmm(File params){
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
                            states = new String[3];
                            states[0] = "normal"; states[1] = "del1"; states[2] = "dup1";                                                        
                        }
                        //Diploid
                        else{
                            states = new String[5];
                            states[0] = "normal"; states[1] = "del1"; states[2] = "del2"; 
                            states[3] = "dup1"; states[4] = "dup2";
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
        avgCov = depthCov / readSize;
    }

    //Creates Transition Matrix and Sets Initial Probabilities for each State
    public void initializeTransitionMatrix(){
        int numStates = states.length;
        initProb = new double[numStates];
        
        //Create Transition Matrix
        transitionMatrix = new double[numStates][numStates];
        for (int i=0; i<numStates; i++){
            initProb[i] = Math.log(1.0/numStates);
            for (int j=0; j<numStates; j++){
                if (i == j){
                    if (states[i].equalsIgnoreCase("normal")){
                        transitionMatrix[i][j] = Math.log(1 - (numCNVs/chrLength));
                    }
                    else{
                        transitionMatrix[i][j] = Math.log(1 - (1.0/avgCNVSize));
                    }
                }
                else{
                    if (states[i].equalsIgnoreCase("normal")){
                        transitionMatrix[i][j] = Math.log((numCNVs/chrLength)/(numStates - 1.0));
                    }
                    else{
                        transitionMatrix[i][j] = Math.log((1.0/avgCNVSize)/(numStates - 1.0));
                    }
                }
            }
        }        
    }

    //Create Emission Matrix
    //emission = Pr(depth-coverage | copy-number) * Pr(distance | state)
    //add a uniform component unif:
    //Pr(distance | CNV neutral) = unif + (1-unif) * Norm(distance | CNV-neutral)
    public void initializeEmissionMatrix(){
        double lambda = avgCov;
        int numStates = states.length;
        double del1; double del2; double dup1; double dup2;
        if (genome.equalsIgnoreCase("d")){
            del1 = lambda/2; del2 = lambda/100;
            dup1 = lambda * 1.5; dup2 = lambda * 2.0;
            emissionMatrix = new double[numStates][maxEmission+1];
            for (int i = 0; i<maxEmission; i++){
                emissionMatrix[0][i] = logPoisson(i, lambda);
                emissionMatrix[1][i] = logPoisson(i, del1);
                emissionMatrix[2][i] = logPoisson(i, del2);
                emissionMatrix[3][i] = logPoisson(i, dup1);
                emissionMatrix[4][i] = logPoisson(i, dup2);
            }
        }
        else{
            del1 = lambda/100; dup1 = lambda * 2;
            del2 = del1; dup2 = dup1;
            emissionMatrix = new double[numStates][maxEmission];
            for (int i = 0; i<maxEmission; i++){
                emissionMatrix[0][i] = logPoisson(i, lambda);
                emissionMatrix[1][i] = logPoisson(i, del1);
                emissionMatrix[2][i] = logPoisson(i, dup1);
            }
        }
    }
    
    //log(Poisson(k,lambda)) = klog(lambda) - lambda - SUM(j)_1..k
    public double logPoisson(double k, double l){
        double p = k * Math.log(l) - l;
        for (int i = 1; i<=k; i++){
            p += -1.0 * Math.log(i);
        }
        return p;
    }

    //Accessor for Emission Matrix
    public double Emission(byte state, byte cov){
        if (cov > maxEmission){
            cov = (byte) maxEmission;
        }
        return emissionMatrix[state][cov];
    }

    //Viterbi Algorithm for Optimal Path    
    public void runViterbiAlgorithm(Chromosome chr, String name){
        int chrSize = chr.size; //Chromosome length
        int numStates = states.length;  //Number of States
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
            delta[i] = initProb[i] + Emission(i, chr.coverage[0]);
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
                delta[j] = prob[prevMap[i][j]] + Emission((byte)j,chr.coverage[i]);
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


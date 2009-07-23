package cnv_hmm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.util.StringTokenizer;
import java.util.ArrayList;
import java.math.BigDecimal;

public class Cnv_Hmm {
    ArrayList<State> states = new ArrayList<State>();    //States for either Diploid or Haploid

    /*  Diploid: 0-Normal 1-del1 2-del2 3-dup1 4-dup2
     *  Haploid: 0-Normal 1-del1 2-dup1
     */
    double[][] transitionMatrix;    //Transition Probabilities
    double[][] emissionCoverageMatrix;      //Emission Coverage Probabilities
    double[][] emissionDistanceMatrix;       //Emission Distance Proabbilities
    
    //Default Values: Potentially overidden by parameter file
    int maxCoverage = 100;  //Maximum emission, if larger, emission is set to 100
    int maxDistance = 1000; //Maximum paired distance, if larger, set to 1000
    double avgCov = 1.0;    //Average Depth Coverage
    String genome = "d";    //Default genome is Diploid
    int numCNVs = 3000;     //3000 CNVs assumed
    double avgCNVSize = 10000.0;    //Assumed size of CNVs
    double readSize = 35.0;               //Read Size
    double depthCov = 1.0;          //Depth Coverage
    double[] initProb;              //Initial Probabilities for each State
    double chrLength = 3000000000.0;
    int maxDeletionSize = 1000;    //Maximum Deletion Size
    int numGridStates;
    double factor = 1.5;
    int StatesPerDelSize = 10;
    double remainbreakpointProb = 0.1;   
    ArrayList<Integer> DelSizes = new ArrayList<Integer>();

    //Default Constructor without Parameter File
    public Cnv_Hmm(){
        initializeGridParams();
        avgCov = depthCov / readSize;   //Set average coverage
    }
    //Parameter File Provided
    public Cnv_Hmm(File params){
        initializeGridParams();
        readParameterFile(params);
        avgCov = depthCov / readSize;   //Set average coverage
    }

    public void initializeGridParams(){
        //Initialize Grid State Parameters- number of Grid States, deletion sizes
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

    public void initializeHaploidStates(){
        //Original States
        states.add(new State("normal"));
        states.add(new State("del1"));
        states.add(new State("dup1"));
        //Grid States
        for (int i = 0; i<numGridStates; i++){
            for (int j = 0;j<StatesPerDelSize; j++){
                //name, isGrid, delSize, gridNumber
                states.add(new State("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j, numGridStates));
            }
        }
        //Breakpoint state
        states.add(new State("bp"));
        //Set initial probabilities: for each state- 1/number of states
        int numStates = states.size();
        initProb = new double[numStates];
        for (int i = 0; i<numStates; i++){
            initProb[i] = Math.log(1.0/numStates);
        }
    }

    public void initializeDiploidStates(){
        //Original States
        states.add(new State("normal"));
        states.add(new State("del1"));
        states.add(new State("del2"));
        states.add(new State("dup1")); 
        states.add(new State("dup2"));
        //Grid States
        for (int i = 0; i<numGridStates; i++){            
            for (int j = 0;j<StatesPerDelSize; j++){
                //name, isGrid, delSize, gridNumber
                states.add(new State("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j, numGridStates));
            }
        }
        //Breakpoint state
        states.add(new State("bp"));
        //Set initial probabilities: for each state- 1/number of states
        int numStates = states.size();
        initProb = new double[numStates];
        for (int i = 0; i<numStates; i++){
            initProb[i] = Math.log(1.0/numStates);
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
    
    public void printTransitionMatrix(){
        try{     
            FileWriter fstream = new FileWriter("TransitionMatrix.txt");
            BufferedWriter out = new BufferedWriter(fstream);
            for (int i=0; i<transitionMatrix.length; i++){
                for (int j=0; j<transitionMatrix[0].length; j++){
                    double prob = transitionMatrix[i][j];
                    if (prob != Double.NEGATIVE_INFINITY){
                        BigDecimal bd = new BigDecimal(prob);
                        bd = bd.setScale(6, BigDecimal.ROUND_DOWN);
                        out.write(bd.toString());
                    }
                    else{
                        out.write(Double.toString(prob));
                    }
                    out.write("\t");
                    if (j == transitionMatrix.length-1){
                        out.write("\n");
                    }
                }
            } 
            out.close();
        }catch (IOException e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }   
    }
   
    //Accessor for Emission Matrix
    public double Emission(byte state, byte cov, int distance){
        if (cov > maxCoverage){
            cov = (byte) maxCoverage;
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
            delta = new double[numStates];
        }
        else{
            delta = new double[numStates];
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
                //Breakpoints/Grid States
                if (cnv[2] > 2){
                    if (cnv[2] == numStates - 1){
                        System.out.println("BP");
                    }
                    else{
                        int index = (cnv[2]-3)/10;
                        System.out.println(DelSizes.get(index));
                    }
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
                //Breakpoints/Grid States
                if (cnv[2] > 4){
                    if (cnv[2] == numStates - 1){
                        System.out.println("BP");
                    }
                    else{
                        int index = (cnv[2]-3)/10;
                        System.out.println(DelSizes.get(index));
                    }
                }
            }
            System.out.printf("\t%10d", cnv[0] + 1);  //Increment by 1 to get correct start location
            System.out.printf("\t%10d", cnv[1] + 1);  //Increment by 1 to get correct end location
            System.out.printf("\t%10d", cnv[1] - cnv[0]);   //Length of CNV
            System.out.println();            
        }

    }
}


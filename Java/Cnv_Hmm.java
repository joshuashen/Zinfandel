

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
    //ArrayList<States> states = new ArrayList<States>();
    ArrayList<State> states = new ArrayList<State>();       //States for either Diploid or Haploid

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
    double[] initProb;              //Initial Probabilities for each States
    double chrLength = 3000000000.0;
//    int maxDeletionSize = 1000;    //Maximum Deletion Size
    int numGridStates;
    double factor = 2;
    int StatesPerDelSize = 10;
    double remainbreakpointProb = 0.10;
    ArrayList<Integer> DelSizes = new ArrayList<Integer>();

    //Default Constructor without Parameter File
    public Cnv_Hmm(){
        initializeGridParams();
        avgCov = depthCov / readSize;   //Set average coverage
        //Initialize States
        if (genome.equalsIgnoreCase("h")){
            initializeHaploidStates();
        }
        else{
            initializeDiploidStates();
        }
    }
    //Parameter File Provided
    public Cnv_Hmm(File params){
        initializeGridParams();
        readParameterFile(params);
        avgCov = depthCov / readSize;   //Set average coverage
        //Initialize States
        if (genome.equalsIgnoreCase("h")){
            initializeHaploidStates();
        }
        else{
            initializeDiploidStates();
        }
    }

    private void initializeGridParams(){
        //Initialize Grid States Parameters- number of Grid States, deletion sizes
        int count = 0;
        double value = 100;
        DelSizes.add((int)value);
        while(value <= maxDistance){
            value = value * factor;
            DelSizes.add((int)value);
            count++;
        }
        numGridStates = count;
    }

    private void initializeHaploidStates(){
        //Original States
        states.add(new normalState());
        states.add(new Deletion1State());
        states.add(new Duplication1State());
        //Grid States
        for (int i = 0; i<numGridStates; i++){
            for (int j = 0;j<StatesPerDelSize; j++){
                //name, isGrid, delSize

                //Initial Grid State
                if (j == 0){
                    states.add(new InitialGridState("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j));
                }
                //Final Grid State
                else if(j == StatesPerDelSize-1){
                    states.add(new FinalGridState("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j));
                }
                //Middle Grid States
                else{
                    states.add(new GridState("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j));
                }
            }
        }
        //Breakpoint state
        //states.add(new BreakpointState());
        //Set initial probabilities: for each state- 1/number of states
        int numStates = states.size();
        initProb = new double[numStates];
        for (int i = 0; i<numStates; i++){
            initProb[i] = Math.log(1.0/numStates);
        }
        System.out.println("States Initialized");
    }

    private void initializeDiploidStates(){
        //Original States
        states.add(new normalState());
        states.add(new Deletion1State());
        states.add(new Deletion2State());
        states.add(new Duplication1State());
        states.add(new Duplication2State());
        //Grid States
        for (int i = 0; i<numGridStates; i++){
            for (int j = 0;j<StatesPerDelSize; j++){
                //name, isGrid, delSize

                //Initial Grid State
                if (j == 0){
                    states.add(new InitialGridState("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j));
                }
                //Final Grid State
                else if(j == StatesPerDelSize-1){
                    states.add(new FinalGridState("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j));
                }
                //Middle Grid States
                else{
                    states.add(new GridState("delSize" + DelSizes.get(i).toString() + "-" + j, DelSizes.get(i), j));
                }
            }
        }
        //Breakpoint state
        //states.add(new BreakpointState());
        //Set initial probabilities: for each state- 1/number of states
        int numStates = states.size();
        initProb = new double[numStates];
        for (int i = 0; i<numStates; i++){
            initProb[i] = Math.log(1.0/numStates);
        }
        System.out.println("States Initialized");
    }

    public void createTransitionMatrix(){
        TransitionMatrix tm = new TransitionMatrix(states, genome, numCNVs, chrLength, avgCNVSize,
                                                   remainbreakpointProb, StatesPerDelSize, DelSizes);
        transitionMatrix = tm.transitionMatrix;
    }

    public void createCoverageMatrix(){
        EmissionCoverageMatrix ecm = new EmissionCoverageMatrix(states, genome, avgCov, maxCoverage);
        emissionCoverageMatrix = ecm.emissionCoverageMatrix;
    }

    public void createDistanceMatrix(int avg, int stddev){
        EmissionDistanceMatrix edm = new EmissionDistanceMatrix(states, genome, avg, stddev, maxDistance);
        emissionDistanceMatrix = edm.emissionDistanceMatrix;
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
                    else if (parameterName.equalsIgnoreCase("GridNum")){
                        StatesPerDelSize = Integer.valueOf(parameterData);
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
            return emissionCoverageMatrix[state][cov] + Math.log(1.0/1000);
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

	 //High Memory Usage
	 System.out.println("OK before prevMap");
        byte[][] prevMap = new byte[chrSize + 1][numStates];    //Holds Optimal Previous state
        //double[] deltaValues = new double[chrSize];
	 System.out.println("OK after prevMap");

        //Set Initial Probabilities and prevMap
        for (byte i = 0; i<numStates; i++){
            delta[i] = initProb[i] + Emission(i, chr.coverage[0], chr.distances[0]);
            prevMap[0][i] = i;
            //deltaValues[0] = delta[0];
        }

        double[] prob = new double[numStates];

        for (int i = 1; i<chrSize; i++){
            for (int j = 0; j<numStates; j++){     //Current States
                    State st = states.get(j);
                    for (int k = 0; k<prob.length; k++){
                        prob[k] = Double.NEGATIVE_INFINITY;
                    }
                    for (int l = 0; l<st.transitions.size(); l++){
                        Transition tr = st.transitions.get(l);
                        int index = states.indexOf(tr.startState);
                        prob[index] = delta[index] + tr.transitionProb;
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
                    /*
                    if (j != 0 && maxState == 0){
                        System.out.println(i + "\t" + delta[0]);
                    }
                    */
                    
                    prevMap[i][j] = maxState;   //Record Optimal State
                    delta[j] = prob[prevMap[i][j]] + Emission((byte)j,chr.coverage[i],chr.distances[i]);
                    //deltaValues[i] = delta[0];
                }
        }
        //Find Final Optimal States
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
//        int[] bestPath = new int[chr.size]; //Optimal Path

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
        System.out.println("\nCNVs: chromosome name, type, start point, end point, cnv length, deltaS\n");
        for (int i = cnvs.size()-1; i >= 0; i--){
            boolean GridCNV = false;
            int [] cnv = cnvs.get(i);
            
            //Haploid Case
            if (genome.equalsIgnoreCase("h")){
                switch(cnv[2]){
                    case 1: System.out.print(name+"\t"); System.out.print("del1"); break;
                    case 2: System.out.print(name+"\t"); System.out.print("dup1"); break;
                }
                //Breakpoints/Grid States
                if (cnv[2] > 2){

	         /*
                    if (cnv[2] == numStates - 1){
                        System.out.println("BP");
                    }
		  */
                    //else{
                        //Downcast-need to ensure no ClassCastException
                        if (states.get(cnv[2]) instanceof GridState){
                            GridCNV = true;
                            GridState st = (GridState) states.get(cnv[2]);
                            if (i - (StatesPerDelSize-1) >= 0){
                                System.out.print(name+"\t");
                                System.out.print(st.delSize);
                                int [] FinalGridCNV = cnvs.get(i - (StatesPerDelSize-1));
                                int startLoc = cnv[0];
                                int endLoc = FinalGridCNV[1];
                                System.out.printf("\t%10d", startLoc + 1);  //Increment by 1 to get correct start location
                                System.out.printf("\t%10d", endLoc + 1);  //Increment by 1 to get correct end location
                                System.out.printf("\t%10d", endLoc - startLoc);   //Length of CNV
                                System.out.print("\t");
                                //System.out.print(deltaValues[startLoc] - deltaValues[endLoc]);   //delta S
                                
                                System.out.println();
                                i = i - (StatesPerDelSize-1);
                            }
                        }
                    //}
                }

            }
            //Diploid Case
            else{
                switch(cnv[2]){
                    case 1: System.out.print(name+"\t"); System.out.print("del1"); break;
                    case 2: System.out.print(name+"\t"); System.out.print("del2"); break;
                    case 3: System.out.print(name+"\t"); System.out.print("dup1"); break;
                    case 4: System.out.print(name+"\t"); System.out.print("dup2"); break;
                }
                //Breakpoints/Grid States
                if (cnv[2] > 4){
                    //GridState st = (GridState) states.get(cnv[2]);
                    //System.out.print(name+"\t");
                    //System.out.print(st.delSize + "-" + st.gridNumber);
/*
                    if (cnv[2] == numStates - 1){
                        System.out.println("BP");
                    }
*/
//                    else{
                        //Downcast-need to ensure no ClassCastException

                        if (states.get(cnv[2]) instanceof GridState){
                            GridCNV = true;
                            GridState st = (GridState) states.get(cnv[2]);
                            if (i - (StatesPerDelSize-1) >= 0){
                                System.out.print(name+"\t");
                                System.out.print(st.delSize);
                                int [] FinalGridCNV = cnvs.get(i - (StatesPerDelSize-1));
                                int startLoc = cnv[0];
                                int endLoc = FinalGridCNV[1];
                                System.out.printf("\t%10d", startLoc + 1);  //Increment by 1 to get correct start location
                                System.out.printf("\t%10d", endLoc + 1);  //Increment by 1 to get correct end location
                                System.out.printf("\t%10d", endLoc - startLoc);   //Length of CNV
                                System.out.print("\t");
                                //System.out.print(deltaValues[startLoc] - deltaValues[endLoc]);   //delta S
                                
                                System.out.println();
                                i = i - (StatesPerDelSize-1);
                            }
                        }

//                    }
                }
            }
            if (!GridCNV){
            //if (cnv[2] != 4){
                System.out.printf("\t%10d", cnv[0] + 1);  //Increment by 1 to get correct start location
                System.out.printf("\t%10d", cnv[1] + 1);  //Increment by 1 to get correct end location
                System.out.printf("\t%10d", cnv[1] - cnv[0]);   //Length of CNV
                System.out.print("\t");
                //System.out.print(deltaValues[cnv[0]] - deltaValues[cnv[1]]);   //delta S
                
                System.out.println();
            }
        }
    }
}


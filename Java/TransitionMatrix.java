package cnv_hmm;
import java.util.ArrayList;

public class TransitionMatrix {
    double[][] transitionMatrix;
    ArrayList<State> states;
    /*
    Necessary Parameters: States, number of CNVs, genome(haploid/diploid), chromosome length, average CNV Size,
                          probability of remaining in breakpoint state, number of grid states per deletion size
                          deletion sizes
    */
    public TransitionMatrix(ArrayList<State> states, String genome, int numCNVs, double chrLength, double avgCNVSize,
                           double remainbreakpointProb, int StatesPerDelSize, ArrayList<Integer> DelSizes){
        this.states = states;
        int numStates = states.size();

        int numNonGridStates = 0;

        //3 if haploid- normal/del1/dup1
        if (genome.equalsIgnoreCase("h")){
            numNonGridStates = 3;
        }//5 if diploid- normal/del1/del2/dup1/dup2
        else{
            numNonGridStates = 5;
        }

        //Create Transition Matrix and initialize all values to negative infinity
        transitionMatrix = new double[numStates][numStates];
        for (int i=0; i<numStates; i++){
            for (int j=0; j<numStates; j++){
                transitionMatrix[i][j] = Double.NEGATIVE_INFINITY;
            }
        }

        for (int i=0; i<numStates; i++){
            for (int j=0; j<numStates; j++){
                //Probability of remaining in same state
                if (i == j){
                    //Probability of remaining in Normal state
                    if (states.get(i) instanceof normalState){
                        transitionMatrix[i][j] = Math.log(1 - (numCNVs/chrLength));
                    }
                    //Probability of remaining in Del1/Del2/Dup1/Dup2
                    else if (states.get(i) instanceof Deletion1State || states.get(i) instanceof Duplication1State ||
                             states.get(i) instanceof Deletion2State || states.get(i) instanceof Duplication2State){
                        transitionMatrix[i][j] = Math.log(1 - (1.0/avgCNVSize));
                    }
                    //Probability of remaining in breakpoint state
                    /*
                    else if (states.get(i) instanceof BreakpointState){
                        transitionMatrix[i][j] = Math.log(remainbreakpointProb);
                    }
                    */
                    //Pr(Transition from grid state to itself)
                    else if (states.get(i) instanceof GridState){
                        transitionMatrix[i][j] = Math.log(1 - ((double)StatesPerDelSize/((GridState)states.get(i)).delSize));
                    }
                }
                //Probability of transitioning from one state to another
                else{
                    //Normal -> Initial Grid States/Del1/Del2/Dup1/Dup2
                    if (states.get(i) instanceof normalState &&
                       (states.get(j) instanceof InitialGridState || states.get(j) instanceof Deletion1State || 
                        states.get(j) instanceof Duplication1State || states.get(j) instanceof Deletion2State ||
                        states.get(j) instanceof Duplication2State)){
                        //Assumed equal probability of transition from normal to del1, dup1, each initial grid state
                        transitionMatrix[i][j] = Math.log((numCNVs/chrLength)/(numNonGridStates + DelSizes.size() - 1.0));
                    }
                    //Del1/Del1/Dup1/Dup2 -> Normal/Del1/Del2/Dup1/Dup2
                    else if ((states.get(i) instanceof Deletion1State || states.get(i) instanceof Duplication1State ||
                              states.get(i) instanceof Deletion2State || states.get(i) instanceof Duplication2State)
                             && (!(states.get(j) instanceof GridState))){
                        //Assumed equal probability of transition from del1/dup1 to normal/del1/dup1 
                        transitionMatrix[i][j] = Math.log((1.0/avgCNVSize)/(numNonGridStates));
                    }
                    //Final Grid States -> Normal(Exiting Grid States back to Normal )
                    else if(states.get(i) instanceof FinalGridState && states.get(j) instanceof normalState){
                        //Prob = number of grid states/deletion size
                        transitionMatrix[i][j] = Math.log((double)StatesPerDelSize/((FinalGridState)states.get(i)).delSize);
                    }
                    /*
                    else if(states.get(i).isBreakpoint){
                        //Equal probability of transition to Normal/Del1/Dup1
                        if (genome.equalsIgnoreCase("h")){
                            transitionMatrix[i][j] = (1-remainbreakpointProb)/3.0;
                        }
                        //Equal probability of transition to Normal/Del1/Dup1/Del2/Dup2
                        else if(genome.equalsIgnoreCase("d")){
                            transitionMatrix[i][j] = (1-remainbreakpointProb)/5.0;
                        }
                    }
                    */
                    //Probability of transitioning to next grid state
                    else if(states.get(i) instanceof GridState && states.get(j) instanceof GridState){
                        if (((GridState)states.get(i)).delSize == ((GridState)states.get(j)).delSize){
                            if (((GridState)states.get(i)).gridNumber + 1 == ((GridState)states.get(j)).gridNumber){
                                transitionMatrix[i][j] = Math.log((double)StatesPerDelSize/((GridState)states.get(i)).delSize);
                            }
                        }
                    }
                }
            }
        }
        SetTransitions();
        System.out.println("Transition Matrix Created");
    }

    public void SetTransitions(){
        for (int i = 0; i<states.size(); i++){
            for (int j = 0; j<states.size(); j++){
                if(transitionMatrix[j][i] != Double.NEGATIVE_INFINITY){
                    Transition t = new Transition(states.get(j), states.get(i), transitionMatrix[j][i]);
                    states.get(i).addTransitionToState(t);
                }
            }
        }
    }

}

package cnv_hmm;
import java.util.ArrayList;

public class TransitionMatrix {
    double[][] transitionMatrix;
    /*
    Necessary Parameters: States, number of CNVs, genome(haploid/diploid), chromosome length, average CNV Size,
                          probability of remaining in breakpoint state, number of grid states per deletion size
                          deletion sizes
    */
    public TransitionMatrix(ArrayList<State> states, String genome, int numCNVs, int chrLength, int avgCNVSize,
                           double remainbreakpointProb, int StatesPerDelSize, ArrayList<Integer> DelSizes){
        int numStates = states.size();
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
                    if (states.get(i).isNormal){
                        transitionMatrix[i][j] = Math.log(1 - (numCNVs/chrLength));
                    }
                    //Probability of remaining in Del1/Del2/Dup1/Dup2
                    else if (states.get(i).isDeletion1 || states.get(i).isDuplication1 ||
                             states.get(i).isDeletion2 || states.get(i).isDuplication2){
                        transitionMatrix[i][j] = Math.log(1 - (1.0/avgCNVSize));
                    }
                    //Probability of remaining in breakpoint state
                    else if (states.get(i).isBreakpoint){
                        transitionMatrix[i][j] = Math.log(remainbreakpointProb);
                    }
                    //Pr(Transition from grid state to itself)
                    else if (states.get(i).isGridState){
                        transitionMatrix[i][j] = Math.log(1-(double)StatesPerDelSize/states.get(i).delSize);
                    }
                }
                //Probability of transitioning from one state to another
                else{
                    //Normal -> Initial Grid State/Del1/Dup1
                    if (states.get(i).isNormal &&
                       (states.get(j).isInitialGridState || states.get(j).isDeletion1 || states.get(j).isDuplication1)){
                        //Assumed equal probability of transition from normal to del1, dup1, each initial grid state
                        transitionMatrix[i][j] = Math.log((numCNVs/chrLength)/(3 + DelSizes.size() - 1.0));
                    }
                    //Del1/Dup1 -> Normal/Del1/Dup1
                    else if ((states.get(i).isDeletion1 || states.get(i).isDuplication1)
                            && (!states.get(j).isGridState)){
                        //Assumed equal probability of transition from del1/dup1 to normal/del1/dup1 grid states
                        transitionMatrix[i][j] = Math.log((1.0/avgCNVSize)/(3.0 - 1.0));
                    }
                    //Final Grid State -> Normal(Exiting Grid States back to Normal )
                    else if(states.get(i).isFinalGridState && states.get(j).isNormal){
                        //Prob = number of grid states/deletion size
                        transitionMatrix[i][j] = Math.log((double)StatesPerDelSize/states.get(i).delSize);
                    }
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
                    //Probability of transitioning to next grid state
                    else if(states.get(i).isGridState && states.get(j).isGridState){
                        if (states.get(i).gridNumber + 1 == states.get(j).gridNumber){
                            transitionMatrix[i][j] = Math.log((double)StatesPerDelSize/states.get(i).delSize);
                        }
                    }
                }
            }
        }
        System.out.println("Transition Matrix Created");
    }
}

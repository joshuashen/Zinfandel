
public class Transition {
    State startState;
    State endState;
    double transitionProb;

    public Transition(State start, State end, double prob){
        startState = start;
        endState = end;
        transitionProb = prob;
    }

}

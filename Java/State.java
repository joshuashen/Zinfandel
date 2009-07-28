import java.util.ArrayList;

public class State {
    String name;
    ArrayList<Transition> transitions ;
    //Default State is Normal
    public State(){
        name = "normal";
        transitions = new ArrayList<Transition>();
    }

    public State(String n){
        name = n;
        transitions = new ArrayList<Transition>();
    }

    public void addTransitionToState(Transition t){
        transitions.add(t);
    }


}

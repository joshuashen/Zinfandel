
package cnv_hmm;


public class State {
    String name;
    boolean isGridState = false;
    boolean isBreakpoint = false;
    boolean isNormal = false;;
    boolean isDeletion1 = false;
    boolean isDuplication1 = false;
    boolean isDeletion2 = false;
    boolean isDuplication2 = false;
    boolean isInitialGridState = false;
    boolean isFinalGridState = false;

    int delSize = -1;
    int gridNumber = -1;
    int maxGridNumber = -1;


    //Constructor for normal/dup/del/bp
    /*Naming Conventions
     * Normal           normal
     * Deletion1        del1
     * Deletion2        del2
     * Duplication1     dup1
     * Duplication2     dup2
     * Breakpoint       bp
     * GridStates       ex. delSize100-0, delSize100-1...delSize100-9
     *                      delSize150-0, delSize150-1...delSize100-9
     *                      ...
    */
    public State(String n){
        name = n;
        setBooleans();
    }

    //Grid State contructor
    public State(String n, int del, int gridNum, int maxGrids){
        name = n;
        delSize = del;
        gridNumber = gridNum;
        maxGridNumber = maxGrids;
        setBooleans();
    }

    public void setBooleans(){
        if (name.equalsIgnoreCase("normal")){
            isNormal = true;
        }
        else if(name.equalsIgnoreCase("del1")){
            isDeletion1 = true;
        }
        else if(name.equalsIgnoreCase("del2")){
            isDeletion2 = true;
        }
        else if(name.equalsIgnoreCase("dup1")){
            isDuplication1 = true;
        }
        else if(name.equalsIgnoreCase("dup2")){
            isDuplication2 = true;
        }
        else if(name.contains("-")){
            isGridState = true;
            if (gridNumber == 0){
                isInitialGridState = true;
            }
            else if(gridNumber == maxGridNumber-1){
                isFinalGridState = true;
            }
        }
        else if(name.equalsIgnoreCase("bp")){
            isBreakpoint = true;
        }
    }

}

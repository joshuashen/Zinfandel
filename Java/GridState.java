package cnv_hmm;
public class GridState extends State{
    int delSize = -1;
    int gridNumber = -1;

    public GridState(){}

    public GridState(String name, int delSize, int gridNumber){
        this.name = name;
        this.delSize = delSize;
        this.gridNumber = gridNumber;                    
    }
}

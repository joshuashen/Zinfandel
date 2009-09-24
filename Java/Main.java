package cnv_hmm;
import java.io.File;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) {
        
        String usage = "USAGE: java Main -r ref.fa -m mapview -p param";
        
        //Given defaults: will be overidden by args
        File mapview = new File("default");
        File reference = new File("default");
        File parameters = new File("default");
//        File output = new File("default");
        
        if (args.length < 6){
            System.out.println(usage);
        }
        else{
            for (int i = 0; i<args.length; i++){
                if (args[i].equalsIgnoreCase("-r")){
                    reference = new File(args[i+1]);
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-m")){
                    mapview = new File(args[i+1]);
                    i++;
                }
                else if(args[i].equalsIgnoreCase("-p")){
                    parameters = new File(args[i+1]);
                    i++;
                }
            }
            if (mapview.exists() && reference.exists() && parameters.exists()){
                Mapping map = new Mapping(mapview, reference);
                map.processFastaFile();
                map.processMapViewFile();
                Cnv_Hmm cnv = new Cnv_Hmm(parameters);
                //cnv.printTransitionMatrix();

                int avg = map.getAverageDistance();
                int dev = map.getStandardDeviationDistance();

                cnv.createTransitionMatrix();
                cnv.createCoverageMatrix();
                cnv.createDistanceMatrix(avg, dev);

                ArrayList<String> keys = map.chromosomeNames;

                for (int i = 0; i<keys.size(); i++){
                    cnv.runViterbiAlgorithm(map.chromosomes.get(keys.get(i)), keys.get(i));
                }
            }
            else{
                System.out.println("Invalid Arguments");
                System.out.println(usage);
            }
        }

    }

}

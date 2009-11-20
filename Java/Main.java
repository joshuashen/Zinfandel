
import java.io.File;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) {

        String usage = "USAGE: java Main -r ref.fa -m mapview -p param -c coverage [-s avgCNVSize] [ -d maxDistance] [ -d minDistance] [-q qualCutoff] [-h head] [-b backgroundDistancePenalty]";

int head4debug = 300000000;
int qualCutoff = 1;
double depthCov = 2;
int avgCNVSize = 1000;
double backgroundDistancePenalty = -4.5;
int maxDisPerPos = 5;
int maxGridDistance = 2000;
int minGridDistance = 150;
int interval = 100;
int lowerBound = 0;
int upperBound = 1000000;

        //Given defaults: will be overidden by args
        File mapview = new File("default");
        File reference = new File("default");
        File parameters = new File("default");
// File output = new File("default");


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
else if(args[i].equalsIgnoreCase("-d")) {
maxGridDistance = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-d2")) {
minGridDistance = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-int")) {
interval = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-q")){
qualCutoff = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-h")){
head4debug = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-b")){
backgroundDistancePenalty = Double.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-c")) {
depthCov = Double.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-s")) {
avgCNVSize = Integer.valueOf(args[i+1]);
i++;
}
//maximum distances per position
else if(args[i].equalsIgnoreCase("-mdpp")) {
maxDisPerPos = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-l")) {
lowerBound = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-u")) {
upperBound = Integer.valueOf(args[i+1]);
i++;
}
            }

            if (mapview.exists() && reference.exists() && parameters.exists()){
                Mapping map = new Mapping(mapview, reference, maxGridDistance, qualCutoff, head4debug, maxDisPerPos, lowerBound, upperBound);
                map.processFastaFile();
                map.processMapViewFile();
                Cnv_Hmm cnv = new Cnv_Hmm(parameters, depthCov, avgCNVSize, maxGridDistance, minGridDistance, interval);
                //cnv.printTransitionMatrix();

                int avg = map.getAverageDistance();
                int dev = map.getStandardDeviationDistance();

                cnv.setbackgroundDistancePenalty(backgroundDistancePenalty);
                cnv.createTransitionMatrix(avg);
                cnv.createCoverageMatrix();
                cnv.createDistanceMatrix(avg, dev);
                cnv.setMaxDisPerDos(maxDisPerPos);

                ArrayList<String> keys = map.chromosomeNames;

                for (int i = 0; i<keys.size(); i++){
                    cnv.runViterbiAlgorithm(map.chromosomes.get(keys.get(i)), keys.get(i), lowerBound);
                }
            }
            else{
                System.out.println("Invalid Arguments");
                System.out.println(usage);
            }
        }

    }

}




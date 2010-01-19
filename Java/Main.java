
import java.io.File;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) {

        String usage = "USAGE: java Main -r ref.fa -m mapview -p param -c coverage [-s avgCNVSize] [ -d maxDistance] [ -d2 minDistance] [ -dMax maxSolid] [-q qualCutoff] [-h head] [-b backgroundDistancePenalty]";

int head4debug = 300000000;
int qualCutoff = 1;
double depthCov = 2;
int avgCNVSize = 1000;
double backgroundDistancePenalty = -4.5;
int maxDisPerPos = 5;
int maxDistance = 10000;
int minGridDistance = 150;
int factor = 2;	//TEMPORARY
int interval = 100;
int gapSize = 10000000;
int overlap = 1000000;
int maxGridDistance = 2000;

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
maxDistance = Integer.valueOf(args[i+1]);
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
else if(args[i].equalsIgnoreCase("-g")) {
gapSize = Integer.valueOf(args[i+1]);
i++;
}
else if(args[i].equalsIgnoreCase("-o")) {
overlap = Integer.valueOf(args[i+1]);
i++;
}
            }

            int lowerBound = 0;
            int upperBound = lowerBound + gapSize;

            if (mapview.exists() && reference.exists() && parameters.exists()){
                System.out.println(lowerBound + "," + upperBound);
                Mapping map = new Mapping(mapview, reference, maxDistance, qualCutoff, head4debug, maxDisPerPos, lowerBound, upperBound);
                map.processFastaFile();
                map.processMapViewFile();
                Cnv_Hmm cnv = new Cnv_Hmm(parameters, depthCov, avgCNVSize, maxDistance, minGridDistance, interval, maxGridDistance);
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
                lowerBound = lowerBound - overlap;
                int chrSize = map.getChrSize();
		  boolean finished = false;
           while (lowerBound<chrSize && !finished){
		      lowerBound = lowerBound + gapSize;
		      if (upperBound + gapSize > chrSize){
		          upperBound = chrSize;
			   finished = true;
		      }
		      else{
		      	   upperBound = upperBound + gapSize;
		      }
                    System.out.println(lowerBound + "," + upperBound);
                    map = new Mapping(mapview, reference, maxDistance, qualCutoff, head4debug, maxDisPerPos, lowerBound, upperBound);
                    map.processFastaFile();
                    map.processMapViewFile();
                    for (int i = 0; i<keys.size(); i++){
                        cnv.runViterbiAlgorithm(map.chromosomes.get(keys.get(i)), keys.get(i), lowerBound);
                    }
            }

            }
            else{
                System.out.println("Invalid Arguments");
                System.out.println(usage);
            }
        }

    }

}




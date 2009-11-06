
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.ArrayList;
//import java.util.regex.*;

public class Mapping {
    Hashtable<String, Chromosome> chromosomes; //HashTable of Chromosome Objects:
    ArrayList<String> chromosomeNames = new ArrayList<String>(); //Names of each chromosome
                                        //indexed in order of sequence read from reference file
    File mapview; //ecoli.map.mapview.random_27k
    File reference; //E.coli.K12-MG1655.fasta
    int averageDistance; //Average distance
    int standardDevDistance; //Standard deviation of distance
    int limit; //Max bp number
    int maxDistance;
    int pairMappingQualCutOff;
    int maxDisPerPos;
    // int minDelta; // minimum delta between avg and

    //Initialize Files and new Chromosomes ArrayList
    public Mapping(File mv, File ref, int maxD, int cutoff, int limit, int maxDisPerPos){
        mapview = mv;
        reference = ref;
        this.limit = limit;
        chromosomes = new Hashtable<String, Chromosome>();
        maxDistance = maxD;
        pairMappingQualCutOff = cutoff;
        this.maxDisPerPos = maxDisPerPos;
    }

    public int getAverageDistance(){
        return averageDistance;
    }

    public int getStandardDeviationDistance(){
        return standardDevDistance;
    }

    //Process Chromosomes from the reference file
    public void processFastaFile(){
        String name = ""; //Chromosome Name
        int seqLength = 0; //Sequence Length
        try{
           FileReader reader = new FileReader(reference);
           BufferedReader in = new BufferedReader(reader);
           System.out.println("Loading Reference Chromosomes.");
           String line; //Current line
           line = in.readLine();
           if (line.startsWith(">")){
               StringTokenizer st = new StringTokenizer(line.substring(1), " ");
               name = st.nextToken(); //Name is the string immediately following the ">" symbol
           }
           while ((line = in.readLine()) != null && seqLength < limit){
               if (line.startsWith(">")){ //Start of new Chromosome
                   chromosomes.put(name, new Chromosome(name, seqLength, maxDisPerPos)); //Add Chromosome to HashTable
                   chromosomeNames.add(name);
                   StringTokenizer st = new StringTokenizer(line.substring(1), " ");
                   name = st.nextToken();
                   seqLength = 0; //Reset seqLength for next Chromosome
               }
               else{
                   seqLength += line.length();
               }
           }
// chromosomes.put(name, new Chromosome(name, limit)); //Add Final Chromosome
           chromosomes.put(name, new Chromosome(name, seqLength, maxDisPerPos)); //Add Final Chromosome
           System.out.println(seqLength);
           chromosomeNames.add(name);
           System.out.println("Reference Chromosomes Loaded.");
           in.close();
        }
        catch(IOException e){
           e.printStackTrace();
        }
    }

    public void processMapViewFile(){
        try{
           FileReader reader = new FileReader(mapview);
           BufferedReader in = new BufferedReader(reader);
Hashtable<String, Integer> positiveBuffer = new Hashtable<String, Integer>(); //HashTable of pairs
Hashtable<String, Integer> startPos = new Hashtable<String, Integer>();
// Hashtable<String, Integer> negativeBuffer; //HashTable of pairs
// ArrayList<String> positiveBuffer = new ArrayList<String>();
// ArrayList<String> negativeBuffer = new ArrayList<String>();

           System.out.println("Loading MapView File.");
           String line; //stores line
                int count = 0; //number of reads counted
                long total = 0; //total distance
                String ref = ""; //reference name;
int mappingQual;
int minPairMappingQual;
int flag;
int readSize;
                while((line = in.readLine()) != null){

// StringTokenizer is deapprieciated.
// StringTokenizer st = new StringTokenizer(line, "\t");
String[] cols = line.split("\\s");
String read = cols[0];
ref = cols[1];
int start = Integer.valueOf(cols[2]);
                    String direction = cols[3];
                    int distance = Integer.valueOf(cols[4]);
                    mappingQual = Integer.valueOf(cols[7]); // single read mapping quality, which is a good indicator of repetitiveness

readSize = Integer.valueOf(cols[13]);
                    if (start < limit){
                        chromosomes.get(ref).incrementCoverage(start);
                    }
else {
break;
}

// !!!! incomplete: should deal with negative distance caused by insertions or tandem duplications
// if (distance >= 0 && distance <= maxDistance){
if (direction.equals("+") && distance > 0 && distance <= maxDistance) {
                        String readString = read.substring(0, read.length() - 2);
// if (!negativeBuffer.containsKey(readString)){
positiveBuffer.put(readString, mappingQual);
startPos.put(readString, start);

/*
if (start < limit){
chromosomes.get(ref).setDistance(start, distance);
total+=Math.abs(distance);
count++;
}
*/
// }
                        //else{
// negativeBuffer.remove(readString);
                        //}
                    }
                    else if (direction.equals("-") && distance < 0 && Math.abs(distance) <= maxDistance){
                        String readString = read.substring(0, read.length() - 2);
                        if (positiveBuffer.containsKey(readString)){

minPairMappingQual = Math.min(mappingQual, positiveBuffer.get(readString));

if (minPairMappingQual > pairMappingQualCutOff) {
// chromosomes.get(ref).setDistance(start + readSize + distance, Math.abs(distance));
chromosomes.get(ref).setDistance(startPos.get(readString), Math.abs(distance));
chromosomes.get(ref).setDistance(start, distance);
total += Math.abs(distance);
count ++;
}
                            positiveBuffer.remove(readString);
startPos.remove(readString);

                        }
// else{
// negativeBuffer.add(readString);
                        // }
                    }
                    //int flag = Integer.valueOf(st.nextToken());
                    //chromosomes.get(ref).setFlag(start, flag);
                }

                double avgDistance = ((double) total)/count;
                averageDistance = (int) avgDistance;
                System.out.println(avgDistance + "\t" + averageDistance);
                Chromosome c = chromosomes.get(ref);
                int[][] distances = c.distances;
                long sum = 0;
                int num = 0;
                for (int i = 0; i<distances.length; i++){
                    for (int j =0; j<maxDisPerPos; j++)
                    if (distances[i][j] != 0){
                        sum += ((Math.abs(Math.abs(distances[i][j]) - averageDistance)) *
                               (Math.abs(Math.abs(distances[i][j]) - averageDistance)));
                        num++;
                    }
                }
                double Dev = Math.sqrt((double) Math.abs(sum)/num);
                standardDevDistance = (int) Dev;
                System.out.println(Dev + "\t" + standardDevDistance);

           System.out.println("MapView File Processed");
           in.close();
positiveBuffer.clear();
startPos.clear();
      }
      catch(IOException e){
           e.printStackTrace();
      }
    }
}


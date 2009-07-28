import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.ArrayList;

public class Mapping {
    Hashtable<String, Chromosome> chromosomes;  //HashTable of Chromosome Objects:
    ArrayList<String> chromosomeNames = new ArrayList<String>();    //Names of each chromosome
                                        //indexed in order of sequence read from reference file
    File mapview;                       //ecoli.map.mapview.random_27k
    File reference;                     //E.coli.K12-MG1655.fasta
    int averageDistance;                //Average distance
    int standardDevDistance;            //Standard deviation of distance

    //Initialize Files and new Chromosomes ArrayList
    public Mapping(File mv, File ref){
        mapview = mv;
        reference = ref;
        chromosomes = new Hashtable<String, Chromosome>();
    }

    public int getAverageDistance(){
        return averageDistance;
    }

    public int getStandardDeviationDistance(){
        return standardDevDistance;
    }

    //Process Chromosomes from the reference file
    public void processFastaFile(){
        String name = "";        //Chromosome Name
        int seqLength = 0;  //Sequence Length
        try{
           FileReader reader = new FileReader(reference);
           BufferedReader in = new BufferedReader(reader);
           System.out.println("Loading Reference Chromosomes.");
           String line;    //Current line
           line = in.readLine();
           if (line.startsWith(">")){
               StringTokenizer st = new StringTokenizer(line.substring(1), " ");
               name = st.nextToken();   //Name is the string immediately following the ">" symbol
           }
           while ((line = in.readLine()) != null){
               if (line.startsWith(">")){   //Start of new Chromosome
                   chromosomes.put(name, new Chromosome(name, seqLength));  //Add Chromosome to HashTable
                   chromosomeNames.add(name);
                   StringTokenizer st = new StringTokenizer(line.substring(1), " ");
                   name = st.nextToken();
                   seqLength = 0;   //Reset seqLength for next Chromosome
               }
               else{
                   seqLength += line.length();
               }
           }
           chromosomes.put(name, new Chromosome(name, seqLength));  //Add Final Chromosome
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
           System.out.println("Loading MapView File.");
           String line;    //stores line
                int count = 0;      //number of reads counted
                int total = 0;      //total distance
                String ref = "";    //reference name;
                while((line = in.readLine()) != null){
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    String read = st.nextToken();
                    ref = st.nextToken();
                    int start = Integer.valueOf(st.nextToken());
                    String direction = st.nextToken();
                    int distance = Integer.valueOf(st.nextToken());
                    chromosomes.get(ref).incrementCoverage(start);
                    if (distance > 0 && distance < 1000){
                        chromosomes.get(ref).setDistance(start, distance);
                        total+=Math.abs(distance);
                        count++;
                    }
                    int flag = Integer.valueOf(st.nextToken());
                    chromosomes.get(ref).setFlag(start, flag);
                }
                
                double avgDistance = ((double) total)/count;
                averageDistance = (int) avgDistance;
                Chromosome c = chromosomes.get(ref);
                int[] distances = c.distances;
                long sum = 0;
                int num = 0;
                for (int i = 0; i<distances.length; i++){
                    if (distances[i] != -1){
                        sum += ((Math.abs(distances[i] - averageDistance)) *
                               (Math.abs(distances[i] - averageDistance)));
                        num++;
                    }
                }
                double Dev = Math.sqrt((double) Math.abs(sum)/num);
                standardDevDistance = (int) Dev;

           System.out.println("MapView File Processed");
           in.close();
      }
      catch(IOException e){
           e.printStackTrace();
      }
    }
}

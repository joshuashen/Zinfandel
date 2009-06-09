package cnv_hmm;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Hashtable;
import java.util.StringTokenizer;

public class Mapping {
    Hashtable<String, Chromosome> chromosomes;  //HashTable of Chromosome Objects:
                                        //indexed in order of sequence read from reference file
    File mapview;                       //ecoli.map.mapview.random_27k
    File reference;                     //E.coli.K12-MG1655.fasta

    //Initialize Files and new Chromosomes ArrayList
    public Mapping(File mv, File ref){
        mapview = mv;
        reference = ref;
        chromosomes = new Hashtable<String, Chromosome>();
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
                   StringTokenizer st = new StringTokenizer(line.substring(1), " ");
                   name = st.nextToken();
                   seqLength = 0;   //Reset seqLength for next Chromosome
               }
               else{
                   seqLength += line.length();
               }
           }
           chromosomes.put(name, new Chromosome(name, seqLength));  //Add Final Chromosome
           
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
                while((line = in.readLine()) != null){
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    String read = st.nextToken();
                    String ref = st.nextToken();
                    int start = Integer.valueOf(st.nextToken());
                    if (start < 300000){    //Temporary limit to reduce memory usage
                        chromosomes.get(ref).incrementCoverage(start);
                    }
                }
           System.out.println("MapView File Processed");
           in.close();
      }
      catch(IOException e){
           e.printStackTrace();
      }
    }
}

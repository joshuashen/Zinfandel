package resultstesting;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import java.util.ArrayList;

public class Main {

    public static void main(String[] args) {
      ArrayList<cnv> trueCNVs = new ArrayList<cnv>();
      try{
           File truecnvfile = new File(args[0]);
           FileReader reader = new FileReader(truecnvfile);
           BufferedReader in = new BufferedReader(reader);
           System.out.println("Loading True CNVs.");
           String line;    //stores line
           line = in.readLine();    //ignore first line
                while((line = in.readLine()) != null){
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    String genomeName = st.nextToken();
                    if(genomeName.equalsIgnoreCase("EcoliGenome")){
                        String type = st.nextToken();
                        int startPos = Integer.valueOf(st.nextToken());
                        int endPos = Integer.valueOf(st.nextToken());
                        int size = Integer.valueOf(st.nextToken());
                        trueCNVs.add(new cnv(type, size, startPos, endPos));
                    }
                }
           System.out.println(trueCNVs.size());
           System.out.println("True CNVs Loaded");
           in.close();
      }
      catch(IOException e){
           e.printStackTrace();
      }

      ArrayList<cnv> resultCNVs = new ArrayList<cnv>();
      try{
           File truecnvfile = new File(args[1]);
           FileReader reader = new FileReader(truecnvfile);
           BufferedReader in = new BufferedReader(reader);
           System.out.println("Loading Result CNVs.");
           String line;    //stores line
           line = in.readLine();    //ignore first line
                while((line = in.readLine()) != null){
                    StringTokenizer st = new StringTokenizer(line, "\t");
                    String genomeName = "";
                    if (st.hasMoreTokens()){
                         genomeName = st.nextToken();
                    }
                    if(genomeName.equalsIgnoreCase("EcoliGenome")){
                        String type = st.nextToken();
                        int startPos = Integer.valueOf(st.nextToken().replaceAll(" ", ""));
                        int endPos = Integer.valueOf(st.nextToken().replaceAll(" ", ""));
                        int size = Integer.valueOf(st.nextToken().replaceAll(" ", ""));
                        resultCNVs.add(new cnv(type, size, startPos, endPos));
                    }
                }
           System.out.println("Result CNVs Loaded");
           in.close();
      }
      catch(IOException e){
           e.printStackTrace();
      }

      int numDeletions = Integer.valueOf(args[2]);
      int errorRange = Integer.valueOf(args[3]);

      int Matches = 0;
      int falsePositives = 0;
      int delSizeMatches = 0;


      for (int i = 0; i<trueCNVs.size(); i++){
          boolean valid = false;
          boolean validSize = false;
          for (int j = 0; j<resultCNVs.size(); j++){
              if ((resultCNVs.get(j).startPos >= trueCNVs.get(i).startPos - errorRange) &&
                   resultCNVs.get(j).endPos <= trueCNVs.get(i).endPos + errorRange){
                  if (trueCNVs.get(i).type.equalsIgnoreCase("DEL")){
                    valid = true;
                    String type = resultCNVs.get(j).type;
                    if (!type.equalsIgnoreCase("del1") && !type.equalsIgnoreCase("dup1")){
                        int x = type.indexOf("-");
                        int delSize = Integer.valueOf(type.substring(0, x));
                        if (delSize == trueCNVs.get(i).size){
                            validSize = true;
                        }
                    }
                  }
              }
          }
          if (valid){
              Matches++;
              if (validSize){
                  delSizeMatches++;
              }
          }
      }

      for (int i = 0; i<resultCNVs.size(); i++){
          boolean found = false;
          for (int j = 0; j<trueCNVs.size(); j++){
            if ((resultCNVs.get(i).startPos >= trueCNVs.get(j).startPos - errorRange) &&
                 resultCNVs.get(i).endPos <= trueCNVs.get(j).endPos + errorRange &&
                 trueCNVs.get(j).type.equalsIgnoreCase("DEL")){
                found = true;
            }
          }
          if (!found){
              falsePositives++;
          }
      }

      System.out.println((double)Matches/numDeletions);
      System.out.println((double)falsePositives/resultCNVs.size());
      System.out.println((double)delSizeMatches/Matches);

    }

}

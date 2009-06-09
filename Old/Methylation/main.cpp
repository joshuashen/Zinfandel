//====================================
// ProportionOfCompleteness_AAC_gap_regroup.cpp
/*************************************
*************** NOTE *****************
**************************************
This application do not compute the number of F3 position duplicates. 
The number of R3 position duplicates is the number within each chromosome type. 
The gap-size ranges are determined by the rule: 

Within each chromosome group:
Sort all pairs according to gap size.
Go through all pairs P_i:
if |Gap(P_i) - Gap(P_i+1)| <= 15k,
	then P_i+1 should be in the same cluster as P_i
else
	create a new cluster for P_i+1

Find consisted clusters according to F_prev - R_curr >= F_curr - R_curr - L (L~4kb) within each gap-size group.

After finding potential clusters, check if there exists a gap-size difference larger than RANGE(15000) within each cluster.
If there exists one, then split the cluster into two parts. The splitting method is rearranging the data according to gap-size (from small to large),
and the largest distance(NOT gap-size) between two data in the F-R plot is the splitting point(breaking point).
Do the iterating checking and splitting until there  does not exist any gap-size difference larger than RANGE(15000) within each sub cluster.

NOTE: For the purpose of not overpassing the integer's upbound, I specified that:
	  If ( |F1 - F2| > square root of the integer's upbound ) or ( |R1 - R2| > square root of the integer's upbound ) or ( the integer's upbound - |F1 - F2|^2 < |R1 - R2|^2)
		then there should be a splitting point(breaking point) between the two data.
*/

//====================================
#include"GeneralHeader.h"
#include"SetOfCandidates.h"
#include"time.h"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>

using namespace std;

//std::ofstream fout("/Users/liyefu/Xcode_Code/C++_Code/ProportionOfCompleteness_AAC_gap_regroup/PrintResult_ProportionOfCompleteness_AAC_gap_regroup.txt");

bool openInputJudge = false;

int main (int argc, char * const argv[]) {
	
	// Must provide the whole path of the specified file even if it has been imported into the project
	
	//char* DataFile[2];
	//char* DataFile[1];
	
	//DataFile[0] = "/Users/liyefu/Temp/SOLiD/data/GRACE20070703_1_F3_R3_0_15000.mates"; //original data file 1
	//DataFile[1] = "/Users/liyefu/Temp/SOLiD/data/GRACE20070703_2_F3_R3_0_15000.mates"; //original data file 2
	
	//char* AAAFile = "/Users/liyefu/Xcode_Code/Xcode_atlas-cnv/AAA_lines/AAA_lines.txt"; // the text file of combined AAA lines 
	string ENDING_INPUT = "y";
	
	string inputFileName;
	vector<string> dataFile;
	
	std::cout << "Please type the input files' paths with names included. e.g. /Users/.../GRACE20070703_1_F3_R3_0_15000.mates\n";
	std::cout << "When finishing typing input files, plese type \"" << ENDING_INPUT << "\" and ENTER key.\n\n";
	
	while(1)
	{
		std::cin >> inputFileName;
	
		if(inputFileName == ENDING_INPUT)
		{
			break;
		}
		
		dataFile.push_back(inputFileName);
	}
	
	std::cout <<"\n";
	std::cout << "Please type the output file's path with name included. e.g. /Users/.../Result.txt\n\n";
	
	std::cin >> inputFileName;
	dataFile.push_back(inputFileName);
	std::ofstream fout(dataFile[dataFile.size() - 1].c_str());
	
	
	if(!fout.is_open())
	{
		std::cout <<"\n";
		std::cout << "Wrong path for the output file.\n";
		return 0;
	}
	
	std::cout <<"\n";
	//-------------------------------------------------------------------

	SetOfCandidates can(dataFile, "AAC");
	
	if(!openInputJudge)
	{
		std::cout << "Input files error.\n";
		return 0;
	}
	
	//can.Subset(1, 100000000, 200000000);
	
	/*
	can.Grouping(ByChr(can._groups));
	can.Grouping(ByGapSize(can._groups));
	can.Grouping(ByConsistence(can._groups));
	*/
	
	ByChr c1(can._groups);
	can.grouping(c1);
	
	ByGapSize c2(can._groups);
	can.grouping(c2);
	
	ByConsistence c3(can._groups);
	can.grouping(c3);
	
	//can.ReadAAALines(AAAFile);
	
	can.display(fout);
	
	
	//-------------------------------------------------------------------
	
	//TIME
	time_t timeval;
	timeval=time(NULL);
	fout<<"Time as local time is "<<ctime(&timeval)<<endl;
	std::cout<<"Time as local time is "<<ctime(&timeval)<<endl; 
	//TIME
	
    return 0;
}

//====================================
// SetOfCandidates.cpp
//====================================

#include"GeneralHeader.h"
#include"SetOfCandidates.h"
#include"criterion.h"
#include"LeftAndRight.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

extern bool openInputJudge;

bool judge_for_ignore(char * seq, char fig);

void eat_comments(ifstream& fin) {

	string str; // The line read from the data file
	
	while(1)
	{
		getline(fin, str);
		if(str[0] == '#' && str[1] == '#'){
			break;
		}
	}
}//-----------------------------------

SetOfCandidates::SetOfCandidates(vector<string>& fileName, char* category){ // Select out candidates with specified category, and sort them according to R3 position
	
	string str; // The line read from the data file
	char* str1;
	
	//char seps[] = " "; // space
	char seps[] = "	"; // This is Tab, different from the one above
	
	char* token;
	char* stringArray[COLUMN_NUM]; // There are 11 data items in the original file 
	int arrayCount = 0;
	int vecCount = 0;
	
	bool judge = false; // judge if the number of consecutive 0 or 1 is over the specified number limit 
	                    // in order to ignore reads with 00000 or 11111
	
	vector<Candidate> _data;
						
	//------------------------------------
	
	//for(int file_index = 0; file_index < 2; file_index++)
	for(int file_index = 0; file_index < fileName.size() - 1; file_index++)
	{
	
		ifstream fin(fileName[file_index].c_str());
		
		if(!fin.is_open())
		{
			std::cout << "Unable to open file " << fileName[file_index] << "\n";
			continue;
		}
		else
		{
			openInputJudge = true;
		}
		
		eat_comments(fin);
		
		while(getline(fin, str)){ // Begin with the first data line
	
			judge = false;
	
			str1 = &str[0];
			token = strtok(str1, seps);
			arrayCount = 0;
		
			while(token != NULL)
			{
				stringArray[arrayCount] = token;
				arrayCount++;
				token = strtok(NULL, seps);
			}
		
			if(strcmp(stringArray[CATEGORY_COL], category) == 0){ // specify the category of the candidate, such as AAA, AAB, ...
			
				// if the F3 sequence contains 000000, then ignore the data line ------------
				judge = judge_for_ignore(stringArray[SEQUENCE_F_COL], '0');
			
				if(judge == true)
				{
					continue;
				}
			
				// if the R3 sequence contains 00000, then ignore the data line ------------
				judge = judge_for_ignore(stringArray[SEQUENCE_R_COL], '0');
			
				if(judge == true)
				{
					continue;
				}
			
				// if the F3 sequence contains 111111, then ignore the data line ------------
				judge = judge_for_ignore(stringArray[SEQUENCE_F_COL], '1');
			
				if(judge == true)
				{
					continue;
				}
			
				// if the R3 sequence contains 111111, then ignore the data line ------------
				judge = judge_for_ignore(stringArray[SEQUENCE_R_COL], '1');
			
				if(judge == true)
				{
					continue;
				}
			
			
				_data.push_back(Candidate());
				_data[vecCount].set(stringArray);
			
				vecCount++;
			
			} // if
		} // while
		
	} // for (file_index)
	
	_groups.push_back(_data);
	
}//-----------------------------------

bool judge_for_ignore(char * seq, char fig) // ignore reads with 00000 or 11111
{
	int cnt = 0;
	bool jud = false;
		
	for(int i = 1; i < strlen(seq); i++)
	{
		if(seq[i] == fig)
		{
			cnt++;
		}
		else
		{
			cnt = 0;
		}
			
			
		if(cnt == LIMIT_NUM)
		{
			jud = true;
			break;
		}
			
	}// for
	
	return jud;

}//-----------------------------------


//NEW************************************************

void SetOfCandidates::grouping(Criterion& ByCertainCriterion)
{
	if(_groups.size() == 0)
	{
		std::cout << "There is no data available.\n";
		return;
	}
	
	ByCertainCriterion.criterionGrouping();
	
	ByCertainCriterion.updateGroupList();

}


/*

//read in AAA lines from AAA_lines.txt and classify those lines according to chromosome types
void SetOfCandidates::ReadAAALines(char* AAAFile)
{
	vector<Candidate> unit1;
	
	string str; // The line read from the data file
	char* str1;
	
	//char seps[] = " "; // space
	char seps[] = "	"; // This is Tab, different from the one above
	
	char* token;
	char* stringArray[COLUMN_NUM]; // There are 11 data items in the original file 
	int arrayCount = 0;

	int chr_type = 0;
	
	
	//Initialize the vector for 24 chromosome types ------ start
	for(int i = 0; i <= CHROMOSOME_NUM; i++)
	{
		AAAVector.push_back(unit1);
	}
	//Initialize the vector for 24 chromosome types ------ end
	
	ifstream fin(AAAFile);

	while(getline(fin, str)){ // Begin with the first data line
	
		str1 = &str[0];
		token = strtok(str1, seps);
		arrayCount = 0;
		
		while(token != NULL)
		{
			stringArray[arrayCount] = token;
			arrayCount++;
			token = strtok(NULL, seps);
		}
		
		chr_type = atoi(stringArray[REFERENCE_F_COL]);
		
		AAAVector[chr_type].push_back(Candidate());
		AAAVector[chr_type][AAAVector[chr_type].size() - 1].set(stringArray);
		
	}

}

//return true if the cluster has coincidence area; if true, invoke proportion() to calculate the proportion of completeness

bool SetOfCandidates::FindCoincidenceArea(int cnt_m, int cnt_n, int cnt_i)
{
	int cnt_j; // the counting variant for the number of each cluster
	int first_F = 0; // the coordinate(position) of the F of the first element in the cluster
	int first_R = 0; // the coordinate(position) of the R of the first element in the cluster
	bool differ_from_first = false; // judge if the other elements' coordinates(positions) are different from the first one. If there is any,
									// this variant should be marked as true
	
	int most_left = 0; // the smallest position of R within the cluster
	int most_right = 0; // the largest position of F within the cluster
	int coincidence_left = 0; // the largest position of R in the coincidence area within the cluster
	int coincidence_right = 0; // the smallest position of F in the coincidence area within the cluster
	
	if(_consistedGroups[cnt_m][cnt_n][cnt_i].size() > 1) // judge if this is a cluster with more than one element
	{
		first_F = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate;
		first_R = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g2->coordinate;
		most_left = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g2->coordinate;
		most_right = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate;
		coincidence_left = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g2->coordinate;
		coincidence_right = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate;
		
		for(cnt_j = 1; cnt_j < _consistedGroups[cnt_m][cnt_n][cnt_i].size(); cnt_j++)
		{
			if( (first_F !=  _consistedGroups[cnt_m][cnt_n][cnt_i][cnt_j].g1->coordinate) || (first_R != _consistedGroups[cnt_m][cnt_n][cnt_i][cnt_j].g2->coordinate) )
			{
				differ_from_first = true;
			}
			
			if(most_left > _consistedGroups[cnt_m][cnt_n][cnt_i][cnt_j].g2->coordinate)
			{
				most_left = _consistedGroups[cnt_m][cnt_n][cnt_i][cnt_j].g2->coordinate;
			}
			
			if(most_right < _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate)
			{
				most_right = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate;
			}
			
			if(coincidence_left < _consistedGroups[cnt_m][cnt_n][cnt_i][cnt_j].g2->coordinate)
			{
				coincidence_left = _consistedGroups[cnt_m][cnt_n][cnt_i][cnt_j].g2->coordinate;
			}
			
			if(coincidence_right > _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate)
			{
				coincidence_right = _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->coordinate;
			
			}
		
		}
		
		if(differ_from_first == true) // judge if the other elements' coordinates(positions) are different from the first one
		{
			if(coincidence_left <= coincidence_right) // judge if the area has a real coincidence
			{
				proportion(most_left, most_right, coincidence_left, coincidence_right, _consistedGroups[cnt_m][cnt_n][cnt_i][0].g1->chr);
				
				return true;
			}
			else
			{
				return false;
			}
		
		
		}
		else
		{
			return false;
		}
		
	
	}
	else
	{
		return false;
	}

}

// Reduce the vacancy area by checking AAA lines which could cover the area
void SetOfCandidates::proportion(int most_left, int most_right, int coincidence_left, int coincidence_right, int chromosome)
{
	int cnt_chr; // the counting variant within each chromosome type
	int cnt_v; // the counting variant for VacancyArea
	
	int vacancy = 0; // the vacant length
	
	int coincidence_middle = 0;
	coincidence_middle = (coincidence_left + coincidence_right) / 2;
	
	vector<LeftAndRight> VacancyArea;
	LeftAndRight LeftRightTemp;
	
	// initialize the vacancy area
	VacancyArea.push_back(LeftAndRight());
	VacancyArea[0].LeftRightSet(coincidence_left, coincidence_right);
	//
	
	cover_num_coincidence_AAA = 0;
	cover_num_outer_left = 0;
	cover_num_outer_right = 0;
	cover_num_inner_left = 0;
	cover_num_inner_right = 0;
	cover_num_middle = 0;
	
	outer_left = most_left;
	outer_right = most_right;
	inner_left = coincidence_left;
	inner_right = coincidence_right;
	middle = coincidence_middle;
	
	for(cnt_chr = 0; cnt_chr < AAAVector[chromosome].size(); cnt_chr++)
	{
		
		if( (AAAVector[chromosome][cnt_chr].g2->coordinate <=  most_left) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= most_left) )
		{
			cover_num_outer_left++;
		}
		
		if( (AAAVector[chromosome][cnt_chr].g2->coordinate <=  most_right) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= most_right) )
		{
			cover_num_outer_right++;
		}
		
		if( (AAAVector[chromosome][cnt_chr].g2->coordinate <=  coincidence_left) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= coincidence_left) )
		{
			cover_num_inner_left++;
		}
		
		if( (AAAVector[chromosome][cnt_chr].g2->coordinate <=  coincidence_right) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= coincidence_right) )
		{
			cover_num_inner_right++;
		}
		
		if( (AAAVector[chromosome][cnt_chr].g2->coordinate <=  coincidence_middle) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= coincidence_middle) )
		{
			cover_num_middle++;
		}
		
		// if the AAA line has coincidence area with the specified area (from coincidence_left to coincidence_right), then [record this AAA line into vector AAARelated: NO] then count the number of AAA data
		
		if( (AAAVector[chromosome][cnt_chr].g1->coordinate < coincidence_left) || (AAAVector[chromosome][cnt_chr].g2->coordinate > coincidence_right) )
		{
			continue;
		}
		else	
		{
			cover_num_coincidence_AAA++;
		}
		
		//
		
		for(cnt_v = 0; cnt_v < VacancyArea.size(); cnt_v++)
		{
			if(VacancyArea[cnt_v].right_point < VacancyArea[cnt_v].left_point)
			{
				continue;
			}
		
			if( (AAAVector[chromosome][cnt_chr].g1->coordinate < VacancyArea[cnt_v].left_point) || (AAAVector[chromosome][cnt_chr].g2->coordinate > VacancyArea[cnt_v].right_point) )
			{
				continue;
			}
			
			if( (AAAVector[chromosome][cnt_chr].g2->coordinate < VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate < VacancyArea[cnt_v].right_point) )
			{
				LeftRightTemp.LeftRightSet(AAAVector[chromosome][cnt_chr].g1->coordinate + 1, VacancyArea[cnt_v].right_point);
				VacancyArea.insert(VacancyArea.begin() + cnt_v, LeftRightTemp);
				VacancyArea.erase(VacancyArea.begin() + cnt_v + 1);
				continue;
			}
			
			if( (AAAVector[chromosome][cnt_chr].g2->coordinate <= VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate >= VacancyArea[cnt_v].right_point) )
			{
				VacancyArea.erase(VacancyArea.begin() + cnt_v);
				cnt_v--;
				continue;
				
			}
			
			if( ( (AAAVector[chromosome][cnt_chr].g2->coordinate > VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate < VacancyArea[cnt_v].right_point) )
				|| ( (AAAVector[chromosome][cnt_chr].g2->coordinate == VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate < VacancyArea[cnt_v].right_point) )
				|| ( (AAAVector[chromosome][cnt_chr].g2->coordinate > VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate == VacancyArea[cnt_v].right_point) )
			)
			{
				LeftRightTemp.LeftRightSet(AAAVector[chromosome][cnt_chr].g1->coordinate + 1, VacancyArea[cnt_v].right_point);
				VacancyArea.insert(VacancyArea.begin() + cnt_v, LeftRightTemp);
				
				LeftRightTemp.LeftRightSet(VacancyArea[cnt_v + 1].left_point, AAAVector[chromosome][cnt_chr].g2->coordinate - 1);
				VacancyArea.insert(VacancyArea.begin() + cnt_v, LeftRightTemp);
				
				VacancyArea.erase(VacancyArea.begin() + cnt_v + 2);
				
				break;
			}
			
			if( (AAAVector[chromosome][cnt_chr].g2->coordinate > VacancyArea[cnt_v].left_point) && (AAAVector[chromosome][cnt_chr].g2->coordinate <= VacancyArea[cnt_v].right_point) && (AAAVector[chromosome][cnt_chr].g1->coordinate > VacancyArea[cnt_v].right_point) )
			{
				LeftRightTemp.LeftRightSet(VacancyArea[cnt_v].left_point, AAAVector[chromosome][cnt_chr].g2->coordinate - 1);
				VacancyArea.insert(VacancyArea.begin() + cnt_v, LeftRightTemp);
				VacancyArea.erase(VacancyArea.begin() + cnt_v + 1);
				continue;
			}
			
		} // for (cnt_v)
	
	} // for (cnt_chr)
	
	
	for(cnt_v = 0; cnt_v < VacancyArea.size(); cnt_v++)
	{
		vacancy = vacancy + VacancyArea[cnt_v].right_point - VacancyArea[cnt_v].left_point + 1;
	}
	
	proportion_of_completeness = (double)(coincidence_right - coincidence_left + 1 - vacancy) / (coincidence_right - coincidence_left + 1);
	
}

*/


void SetOfCandidates::display(std::ofstream& fout){

	fout << "group#   key(id)	estimated_gap_size	chr_F	coordinate(position)_F	orientation_F	duplicates_num_F	chr_R	coordinate(position)_R	orientation_R	duplicates_num_R	category(label)\n\n";

/*
	//for(int i = 0; i < _consistedGroups.size(); i++)
	for(int i = 0; i < 1; i++) // temperately print out the result of chromosome 1
	{
		fout << "\n";
		fout << "\n";
		fout << "------------------------------------\n";
		fout << "Chromosome Type : " << (_consistedGroups[i][0][0][0].g1)->chr << "\n";
		
		for(int j = 0; j < _consistedGroups[i].size(); j++)
		{
			for(int k = 0; k < _consistedGroups[i][j].size(); k++)
			{
				
				if(FindCoincidenceArea(i, j, k))
				{
				
					fout << "\n";
				
					fout << "Chromosome Type : " << (_consistedGroups[i][0][0][0].g1)->chr <<  "   ======> Gap-size Group : " << j <<  "   ======> Candidate Group " << k << " :" << "\n";
				
					for(int l = 0; l < _consistedGroups[i][j][k].size(); l++)
					{
						_consistedGroups[i][j][k][l].print(fout);
						
					} // for (int l)
					
					fout << "\n";
					
					fout << "Proportion of completeness : " << proportion_of_completeness << "\n";
					
					fout << "The number of AAA data related to the coincidence : " << cover_num_coincidence_AAA << "\n";
					
					fout << "outer_left_position : " << outer_left << "		coverage number : " << cover_num_outer_left << "\n";
					
					fout << "inner_left_position : " << inner_left << "		coverage number : " << cover_num_inner_left << "\n";
					
					fout << "coincidence_middle_position : " << middle << "		coverage number : " << cover_num_middle << "\n";
					
					fout << "inner_right_position : " << inner_right << "		coverage number : " << cover_num_inner_right << "\n";
					
					fout << "outer_right_position : " << outer_right << "		coverage number : " << cover_num_outer_right << "\n\n";
					
				}
					
			} // for (int k)
			
		} // for (int j)
		
	} // for (int i)
	
	*/
	
	int cnt = 0;
	
	for(lvCanIter it = _groups.begin(); it != _groups.end(); it++)
	{
	
		//fout << "Group " << cnt << " :\n";
		for(int i = 0; i < (*it).size(); i++)
		{
			fout << cnt << "   ";
			(*it)[i].print(fout);
		}
	
		cnt++;
		
		fout << "\n";
	}
}

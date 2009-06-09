#include"GeneralHeader.h"
#include"criterion.h"
#include"candidate.h"
#include <iostream>
#include <vector>
#include <list>

using namespace std;

void Criterion::addNewCandidates(lvCanIter& lit, int ele_index)
{
	(*(--_new_group_list.end())).push_back(Candidate());
	*((*(--_new_group_list.end())).end() - 1) = (*lit)[ele_index];
	
} //-----------------------------------

void Criterion::updateGroupList()
{
	(*_group_list).clear();
	(*_group_list) = _new_group_list;
	_new_group_list.clear();

} //-----------------------------------

//ByChr: construct groups according to chromosome type, and
		 //compute R3 position duplicates within each chromosome type

//ByGapSize: construct groups according to gap sizes within each chromosome group

//ByConsistence: construct groups according to data consistence within each chromosome group

void Criterion::criterionGrouping()
{
	//ByChr: construct groups according to chromosome type
	//ByGapSize: construct data groups according to gap-size difference within each chromosome group(when |Gap(P_i) - Gap(P_i+1)| > 15k, create a new gap-size group)
	for(lvCanIter lit = (*_group_list).begin(); lit != (*_group_list).end(); lit++)
	{
		criterionInitialization(lit);
	
		for(int cnt_i = 1; cnt_i < (*lit).size(); cnt_i++)
		{
			groupingDetails(lit, cnt_i);
		
		} // for (cnt_i)
	
		extraGroupingWork(lit); // this function is special for ByConsistence
	
	} // for (lit)

	//ByChr: compute the number of R3 position duplicates within each chromosome type
	//ByGapSize: // sort candidates according to R3 position(coordinate)  within each gap_size range of each chromosome type, and
				 //sort candidates according to F3 positions(coordinate) which have the same R3 position(coordinate) within each gap_size range of each chromosome type
	afterGroupingProcess(); // this function is for calculating work for ByChr and sorting work for ByGapSize

} //-----------------------------------

//------------------------------------

// criterion initialization ------ start

void BySharing::criterionInitialization(lvCanIter& lit)
{
	vector<Candidate> unit;
	
	_new_group_list.push_back(unit);
	sort((*lit).begin(), (*lit).end(), getSortingCriterion());
	curr_control_variant = getCurrentControlVariant(lit, 0);
	addNewCandidates(lit, 0);

} //-----------------------------------

void ByConsistence::criterionInitialization(lvCanIter& lit)
{
	F_prev = (*lit)[0].g1->coordinate;
	start_point = 0;
	end_point  = 0;
	
} //-----------------------------------

// criterion initialization ------ end

// grouping details ------ start

void BySharing::groupingDetails(lvCanIter& lit, int ele_index)
{
	judgeIfAddNewVector(lit, ele_index);
	addNewCandidates(lit, ele_index);
	
} //-----------------------------------

void ByConsistence::groupingDetails(lvCanIter& lit, int ele_index)
{
	DelSize = (*lit)[ele_index].g1->coordinate - (*lit)[ele_index].g2->coordinate - LENGTH;
		
	dis = F_prev - (*lit)[ele_index].g2->coordinate;
		
	if(dis >= 0 && dis >= DelSize)
	{
		end_point++;
	}
	else
	{	
		regroupCluster(lit, start_point, end_point);
					
		start_point = end_point + 1;
		end_point  = start_point;
	}
		
	F_prev = (*lit)[ele_index].g1->coordinate;

} //-----------------------------------

// grouping details ------ end

// extra grouping work ------ start

void BySharing::extraGroupingWork(lvCanIter& lit) { }

void ByConsistence::extraGroupingWork(lvCanIter& lit)
{
	regroupCluster(lit, start_point, end_point);
	
} //-----------------------------------

// extra grouping work ------ end

// after grouping process ------ start

void BySharing::afterGroupingProcess()
{
	for(lvCanIter nlIter = _new_group_list.begin(); nlIter != _new_group_list.end(); nlIter++)
	{
		//ByChr: sort each vector in _new_group_list according to R3 position(coordinate)
		//ByGapSize: //sort each vector in _new_group_list according to R3 position(coordinate)
		sort((*nlIter).begin(), (*nlIter).end(), sort_by_R3_position);
		
		element_begin = 0;
		element_end = 0;
		curr_coordinate = ((*nlIter)[0].g2)->coordinate;
		
		chrDuplicateInitialization();
		
		for(int cnt_i = 1; cnt_i < (*nlIter).size(); cnt_i++)
		{
			if(((*nlIter)[cnt_i].g2)->coordinate == curr_coordinate)
			{
				element_end++;
				
				chrDuplicateAccumulation();
	
			}
			else
			{
				dupliSortProcess(nlIter);
				
				element_begin = cnt_i;
				element_end = cnt_i;
				
				chrDuplicateInitialization();
				
				curr_coordinate = ((*nlIter)[cnt_i].g2)->coordinate;

			}
		
		}// for (cnt_i)
		
		dupliSortProcess(nlIter);
	
	}//for (nlIter)

} //-----------------------------------

void ByConsistence::afterGroupingProcess() { }

// after grouping process ------ end

//------------------------------------

SortFuncPtr ByChr::getSortingCriterion()
{
	return sort_by_F3_reference;
}

SortFuncPtr ByGapSize::getSortingCriterion()
{
	return sort_by_gap_size;
}

//------------------------------------

int ByChr::getCurrentControlVariant(lvCanIter& lit, int ele_index)
{
	return ((*lit)[ele_index].g1)->chr;	
}


int ByGapSize::getCurrentControlVariant(lvCanIter& lit, int ele_index)
{
	return ((*lit)[ele_index].g1)->coordinate - ((*lit)[ele_index].g2)->coordinate - LENGTH;
}

//------------------------------------

void ByChr::judgeIfAddNewVector(lvCanIter& lit, int ele_index)
{
	vector<Candidate> unit;

	if(((*lit)[ele_index].g1)->chr != curr_control_variant)
	{
		_new_group_list.push_back(unit);
	}
					
	curr_control_variant = ((*lit)[ele_index].g1)->chr;

}

void ByGapSize::judgeIfAddNewVector(lvCanIter& lit, int ele_index)
{
	vector<Candidate> unit;
	
	if((((*lit)[ele_index].g1)->coordinate - ((*lit)[ele_index].g2)->coordinate - LENGTH - curr_control_variant) > RANGE )
	{
		_new_group_list.push_back(unit);
	}
					
	curr_control_variant = ((*lit)[ele_index].g1)->coordinate - ((*lit)[ele_index].g2)->coordinate - LENGTH;
	
}

//------------------------------------

void ByChr::chrDuplicateInitialization()
{
	duplicates_num = 1;
}

void ByGapSize::chrDuplicateInitialization() { }

//------------------------------------

void ByChr::chrDuplicateAccumulation()
{
	duplicates_num++;
}

void ByGapSize::chrDuplicateAccumulation() { }

//------------------------------------

void ByChr::dupliSortProcess(lvCanIter& nlIter)
{
	for(int cnt_i =element_begin; cnt_i <= element_end; cnt_i++)
	{
		((*nlIter)[cnt_i].g2)->coordinate_duplicates = duplicates_num;
	}

}

void ByGapSize::dupliSortProcess(lvCanIter& nlIter)
{
	sort((*nlIter).begin() + element_begin, (*nlIter).begin() + element_end + 1, sort_by_F3_position);
}

//------------------------------------


// Cluster Regrouping if certain gap-size difference is larger than RANGE(15000)

void ByConsistence::regroupCluster(lvCanIter& lit, int start, int end)
{
	vector<Candidate> unit;
	
	// sort the potential cluster according to gap-size
	sort((*lit).begin() + start, (*lit).begin() + end + 1, sort_by_gap_size);
	
	if(((*lit)[end].g1->coordinate - (*lit)[end].g2->coordinate) 
		- ((*lit)[start].g1->coordinate - (*lit)[start].g2->coordinate) > RANGE ) // the difference of gap-size is larger than RANGE
	{
		recursiveCount = 0;
		
		recursiveRegrouping(lit, start, end);
	
	
		// import the clusters from tempForRecursion
		
		for(int i = 0; i < tempForRecursion.size(); i++)
		{
			_new_group_list.push_back(unit);
			
			for(int j = 0; j < tempForRecursion[i].size(); j++)
			{
				(*(--_new_group_list.end())).push_back(Candidate());
				*((*(--_new_group_list.end())).end() - 1) = tempForRecursion[i][j];
				
			}
			
		}
		
	}
	else
	{
		_new_group_list.push_back(unit);
		
		for(int i = start; i <= end; i++)
		{
			(*(--_new_group_list.end())).push_back(Candidate());
			*((*(--_new_group_list.end())).end() - 1) = (*lit)[i];
		}
	
	}
	
	tempForRecursion.clear();

} //-----------------------------------


// Carry out the task of recursive regrouping clusters

void ByConsistence::recursiveRegrouping(lvCanIter& lit, int start, int end)
{
	int break_point_l = 0;
	int break_point_r = 0;
	int longest_distance_sqr = 0;

	int f1 = 0;
	int f2 = 0;
	int r1 = 0;
	int r2 = 0;
	
	vector<Candidate> unit;
	
	if((start == end) || (((*lit)[end].g1->coordinate - (*lit)[end].g2->coordinate) 
		- ((*lit)[start].g1->coordinate - (*lit)[start].g2->coordinate) <= RANGE)) // The condition of terminating recursion
	{
		tempForRecursion.push_back(unit);
		
		for(int i = start; i <= end; i++)
		{
			tempForRecursion[recursiveCount].push_back(Candidate());
			*(tempForRecursion[recursiveCount].end() - 1) = (*lit)[i];
		}
		
		recursiveCount++;
	
	}
	else
	{
		// Search for the breakpoint where distance of two points in F-R plot is the longest 
		for(int i = start + 1; i <= end; i++)
		{
			f1 = (*lit)[i].g1->coordinate;
			f2 = (*lit)[i - 1].g1->coordinate;
			r1 = (*lit)[i].g2->coordinate;
			r2 = (*lit)[i - 1].g2->coordinate;
			
			if( abs(f1 - f2) > SQUARE_ROOT_OF_INT_UPBOUND || abs(r1 - r2) > SQUARE_ROOT_OF_INT_UPBOUND
				|| (INT_UPBOUND - (f1 - f2) * (f1 - f2) < (r1 - r2) * (r1 - r2))) // guarantee that the value of distanceSquare is legal
				
			{
				break_point_l = i - 1;
				break_point_r = i;
				break;
			}
			else
			{
				if( distanceSquare(lit, i, i - 1) >= longest_distance_sqr)
				{
					longest_distance_sqr = distanceSquare(lit, i, i - 1);
					break_point_l = i - 1;
					break_point_r = i;
				}
			}
		
		} // for (i)
		
		recursiveRegrouping(lit, start, break_point_l);
				
		recursiveRegrouping(lit, break_point_r, end);
		
	}

} //-----------------------------------

// Calculate the square distance of two points in the F-R plot

int ByConsistence::distanceSquare(lvCanIter& lit, int a, int b)
{
	return (((*lit)[a].g1->coordinate - (*lit)[b].g1->coordinate) * ((*lit)[a].g1->coordinate - (*lit)[b].g1->coordinate) + ((*lit)[a].g2->coordinate - (*lit)[b].g2->coordinate) * ((*lit)[a].g2->coordinate - (*lit)[b].g2->coordinate));
	
} //-----------------------------------


bool sort_by_F3_reference(const Candidate& a, const Candidate& b){ 
	
	return (a.g1)->chr < (b.g1)->chr;
	
} //-----------------------------------

bool sort_by_R3_position(const Candidate& a, const Candidate& b){ 
	
	return (a.g2)->coordinate < (b.g2)->coordinate;
	
} //-----------------------------------


bool sort_by_gap_size(const Candidate& a, const Candidate& b){ 
	
	return (a.g1)->coordinate - (a.g2)->coordinate < (b.g1)->coordinate - (b.g2)->coordinate;
	
} //-----------------------------------

bool sort_by_F3_position(const Candidate& a, const Candidate& b){

	return (a.g1)->coordinate < (b.g1)->coordinate;
	
} //-----------------------------------
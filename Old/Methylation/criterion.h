//====================================
// criterion.h
//====================================
#ifndef CRITERION_H
#define CRITERION_H

#include"GeneralHeader.h"
#include"candidate.h"
#include <vector>
#include <list>
#include <algorithm>

using namespace std;

typedef vector<Candidate> vCandidate;
typedef vCandidate::iterator vCanIter;
typedef list<vCandidate> lvCandidate;
typedef lvCandidate::iterator lvCanIter;

typedef bool (*SortFuncPtr)(const Candidate& a, const Candidate& b); // SortFuncPtr is a pointer to sorting functions

bool sort_by_F3_reference(const Candidate& a, const Candidate& b);
bool sort_by_R3_position(const Candidate& a, const Candidate& b);
bool sort_by_F3_position(const Candidate& a, const Candidate& b);
bool sort_by_gap_size(const Candidate& a, const Candidate& b);

class Criterion{

  protected:
  
	lvCandidate* _group_list;
	lvCandidate _new_group_list;
	
	virtual void criterionInitialization(lvCanIter& lit) = 0;
	virtual void groupingDetails(lvCanIter& lit, int ele_index) = 0;
	virtual void extraGroupingWork(lvCanIter& lit) = 0; // only for ByConsistence
	virtual void afterGroupingProcess() = 0;
	
	void addNewCandidates(lvCanIter& lit, int ele_index);

  public:
	
	Criterion(lvCandidate& _groups) 
	{
		_group_list = &_groups;
	}
	
	void criterionGrouping();
	void updateGroupList();
};

class BySharing : public Criterion{

  protected:

	int curr_control_variant; // ByChr: indicate the current chromosome type 
							  // ByGapSize: the current gap-size
							  
	int element_begin; // ByChr: specify the starting element within each chromosome group
					   // ByGapSize: specify the starting element within each gap_size group of each chromosome group
						   
	int element_end; // ByChr: specify the end element within each chromosome group
					 // ByGapSize: specify the end element within each gap_size group of each chromosome group
						 
	int curr_coordinate; // ByChr: indicate the current R3 position(coordinate) within each chromosome type
						 // ByGapSize: indicate the current R3 position(coordinate) within each gap_size group of each chromosome type

  public:
	
	BySharing(lvCandidate& _groups) : Criterion(_groups) { }
	
	void criterionInitialization(lvCanIter& lit);
	void groupingDetails(lvCanIter& lit, int ele_index);
	void extraGroupingWork(lvCanIter& lit); // only for ByConsistence
	void afterGroupingProcess();
	
	virtual SortFuncPtr getSortingCriterion() = 0;
	virtual int getCurrentControlVariant(lvCanIter& lit, int ele_index) = 0;
	virtual void judgeIfAddNewVector(lvCanIter& lit, int ele_index) = 0;
	virtual void chrDuplicateInitialization() = 0;
	virtual void chrDuplicateAccumulation() = 0;
	virtual void dupliSortProcess(lvCanIter& nlIter) = 0;
};

class ByChr : public BySharing{

  private:
  
	int duplicates_num; // ByChr: indicate the number of duplicates of certain R3 position within each chromosome group

  public:
  
	ByChr(lvCandidate& _groups) : BySharing(_groups), duplicates_num(1) { }
	
	SortFuncPtr getSortingCriterion();
	int getCurrentControlVariant(lvCanIter& lit, int ele_index);
	void judgeIfAddNewVector(lvCanIter& lit, int ele_index);
	void chrDuplicateInitialization();
	void chrDuplicateAccumulation();
	void dupliSortProcess(lvCanIter& nlIter);
};

class ByGapSize : public BySharing{

  public:
  
	ByGapSize(lvCandidate& _groups) : BySharing(_groups) { }
	
	SortFuncPtr getSortingCriterion();
	int getCurrentControlVariant(lvCanIter& lit, int ele_index);
	void judgeIfAddNewVector(lvCanIter& lit, int ele_index);
	void chrDuplicateInitialization();
	void chrDuplicateAccumulation();
	void dupliSortProcess(lvCanIter& nlIter);

};

class ByConsistence : public Criterion{

  private:
  
	int DelSize; // estimated deletion size: DelSize = F_curr - R_curr - L (L~4kb)
	int F_prev;
	int dis; // dis = F_prev - R_curr
	int start_point; // remember the starting point of a potential cluster for regrouping
	int end_point; // remember the ending point of a potential cluster for regrouping
  
  public:
  
	int recursiveCount; // When regrouping clusters, this counting variant is used in recursion
	vector<vector<Candidate> > tempForRecursion;
	
	ByConsistence(lvCandidate& _groups) : Criterion(_groups), DelSize(0),  F_prev(0), dis(0), start_point(0), end_point(0) { }
	
	void criterionInitialization(lvCanIter& lit);
	void groupingDetails(lvCanIter& lit, int ele_index);
	void extraGroupingWork(lvCanIter& lit);
	void afterGroupingProcess();
	
	void regroupCluster(lvCanIter& lit, int start, int end);
	void recursiveRegrouping(lvCanIter& lit, int start, int end);
	int distanceSquare(lvCanIter& lit, int a, int b);
	
};

#endif
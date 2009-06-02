//====================================
// SetOfCandidates.h
//====================================
#ifndef SETOFCANDIDATES_H
#define SETOFCANDIDATES_H

#include"GeneralHeader.h"
#include"candidate.h"
#include"criterion.h"
#include <vector>
#include <list>

using namespace std;

class SetOfCandidates {
  
  public:
  
	lvCandidate _groups;
  
	SetOfCandidates(vector<string>& fileName, char* category);
  
	void grouping(Criterion& ByCertainCriterion);
  
	//read in AAA lines from AAA_lines.txt and classify those lines according to chromosome types
	//void ReadAAALines(char* AAAFile);
  
	void display(std::ofstream& fout);
  
 private:
  
	int outer_left; // the left bound position within each cluster
	int outer_right; // the right bound position within each cluster
	int inner_left; // the left bound position of coincidence area within each cluster
	int inner_right; // the right bound position of coincidence area within each cluster
	int middle; // the middle position between inner_left and inner_right
  
	int cover_num_coincidence_AAA; // the number of AAA data which have coincidence area with the [inner_left, inner_right] of each cluster
	int cover_num_outer_left; // the number of outer_left position coverage
	int cover_num_outer_right; // the number of outer_right position coverage
	int cover_num_inner_left; // the number of inner_left position coverage
	int cover_num_inner_right; // the number of inner_right position coverage
	int cover_num_middle; // the number of middle position coverage
  
	double proportion_of_completeness;
  
	vector<vector<Candidate> > AAAVector;
  
	//vector<Candidate> AAARelated;
  
	//void RegroupCluster(int cnt_m, int cnt_n, int start, int end);
	//void IterRegrouping(int cnt_m, int cnt_n, int start, int end);
	//int DistanceSquare(int cnt_m, int cnt_n, int a, int b);
  
	//return true if the cluster has coincidence area; if true, invoke proportion() to calculate the proportion of completeness
	//bool FindCoincidenceArea(int cnt_m, int cnt_n, int cnt_i);
  
	//void proportion(int most_left, int most_right, int coincidence_left, int coincidence_right, int chromosome); // calculate the proportion of completeness
  
};

#endif

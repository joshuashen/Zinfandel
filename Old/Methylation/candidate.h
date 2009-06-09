//====================================
// candidate.h
//====================================
#ifndef CANDIDATE_H
#define CANDIDATE_H

#include"GeneralHeader.h"
#include"genomic.h"

enum Category{AAA, AAB, AAC, BAA, BAB, BAC, ABA, ABB, ABC, CSS, NNN};
enum ORIGINAL_INPUT_COLUMN_NAMES {ID_COL, SEQUENCE_F_COL, SEQUENCE_R_COL, REFERENCE_F_COL = 6, REFERENCE_R_COL, POSITION_F_COL, POSITION_R_COL, CATEGORY_COL};

class Candidate{
	//char* key;
	//Genomic *g1, *g2;
	Category cat;
	
public:
	char* key;
	Genomic *g1, *g2;
	//Category cat;
	Candidate();
	//~Candidate();
	void set(char* []);
	void print(std::ofstream&);
};	

#endif


//====================================
// candidate.cpp
//====================================

#include"GeneralHeader.h"
#include "candidate.h"
#include <iostream>
#include <string>
#include <iomanip>

using namespace std;

Candidate::Candidate(){
	key = new char[20];
	g1 = new Genomic;
	g2 = new Genomic;
	cat = AAA;
}//------------------------------------

void Candidate::set(char* array[]){
	
	bool judge_F = true;
	bool judge_R = true;
	
	strcpy(key, array[ID_COL]); // beadId; Assign a certain value to a dynamic variable that has been allocated
	judge_F = g1->genomicSet(array[REFERENCE_F_COL], array[POSITION_F_COL]); // the parameters are: F3 reference and F3 position
	judge_R = g2->genomicSet(array[REFERENCE_R_COL], array[POSITION_R_COL]); // the parameters are: R3 reference and R3 position
	
	if(judge_F == false && judge_R == false)
	{
		g1->genomic_neg_FR_exchange(array[POSITION_R_COL]);
		g2->genomic_neg_FR_exchange(array[POSITION_F_COL]);
	
	}
		
	if(strcmp(array[CATEGORY_COL], "AAA") == 0){
		cat = AAA;
	}
	else if(strcmp(array[CATEGORY_COL], "AAB") == 0){
		cat = AAB;
	}
	else if(strcmp(array[CATEGORY_COL], "AAC") == 0){
		cat = AAC;		
	}
	else if(strcmp(array[CATEGORY_COL], "BAA") == 0){
		cat = BAA;		
	}
	else if(strcmp(array[CATEGORY_COL], "BAB") == 0){
		cat = BAB;		
	}
	else if(strcmp(array[CATEGORY_COL], "BAC") == 0){
		cat = BAC;		
	}
	else if(strcmp(array[CATEGORY_COL], "ABA") == 0){
		cat = ABA;		
	}
	else if(strcmp(array[CATEGORY_COL], "ABB") == 0){
		cat = ABB;		
	}
	else if(strcmp(array[CATEGORY_COL], "ABC") == 0){
		cat = ABC;		
	}
	else if(strcmp(array[CATEGORY_COL], "C**") == 0){
		cat = CSS;		
	}
	else{
		cat = NNN;
	}
	
}//------------------------------------

void Candidate::print(std::ofstream& fout){
	
	fout << key;
	
	//---
	fout << setfill(' ');
	fout << setw(15) << g1->coordinate - g2->coordinate - LENGTH;
	//---
	
	g1->genomicPrint(fout);
	g2->genomicPrint(fout);
	
	fout << setfill(' ');
	fout << setw(9) << cat << "\n";
	
}//------------------------------------

/*Candidate::~Candidate(){
	delete g1;
	delete g2;
}//------------------------------------*/
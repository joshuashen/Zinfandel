//====================================
// genomic.h
//====================================
#ifndef GENOMIC_H
#define GENOMIC_H

#include"GeneralHeader.h"
#include <iostream>
#include <fstream>

class Genomic{
	//int chr; // 1, 2, 3, ..., 22, 23(X), 24(Y)
	//int coordinate; // max: 200,000,000
	int orientation; // 1 or -1
	//int coordinate_duplicates; //duplicate within each chromosome type
public:
	int chr; // 1, 2, 3, ..., 22, 23(X), 24(Y)
	int coordinate; // max: 200,000,000
	int coordinate_duplicates; //duplicate within each chromosome type
	Genomic();
	bool genomicSet(char*, char*);
	void genomic_neg_FR_exchange(char*);
	void genomicPrint(std::ofstream&);
};	

#endif


//====================================
// genomic.cpp
//====================================

#include"GeneralHeader.h"
#include "genomic.h"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

Genomic::Genomic(){
	chr = 1;
	coordinate = 100000000;
	orientation = 1;
	coordinate_duplicates = 0;
}//------------------------------------

bool Genomic::genomicSet(char* chrType, char* position){

	chr = atoi(chrType);
	coordinate = atoi(position);
	if(*position == '-')
	{
		orientation = -1; // when position is negative
		return false;
	}
	else
	{
		orientation = 1; // when position is positive
		return true;
	}
	
}//------------------------------------

void Genomic::genomic_neg_FR_exchange(char* position)
{
	coordinate = - atoi(position);

}//------------------------------------

void Genomic::genomicPrint(std::ofstream& fout){

	//std::cout << setfill(' ');
	fout << setfill(' ');
	//std::cout << setw(8) << chr << setw(16) << coordinate << setw(8) << orientation;
	fout << setw(8) << chr << setw(16) << coordinate << setw(8) << orientation << setw(8) << coordinate_duplicates;
	
}//------------------------------------

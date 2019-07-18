#include "6mer_index.h"
#include <iostream>

#include "6mer_index_official.rc"
#include <cstring>

/* A = 0; C = 1; G = 2; T = 3 */
int g::Mer2Signal::SixMer2Index(char const* sixmer)
{
// 	const short atag = 0, ctag = 1, gtag = 2, ttag = 3;
	short idx = 0;
	
	for(int i = 0; i < 6; i++){
		short tag;
		switch(sixmer[i]){
		case 'A':
			tag = 0;
			break;
		case 'C':
			tag = 1;
			break;
		case 'G':
			tag = 2;
			break;
		case 'T':
			tag = 3;
			break;
		}
		tag <<= (5-i)*2;
		idx |= tag;
	}
	
	return idx;
}

int g::Mer2Signal::SixMer2Index(char g0, char g1, char g2, char g3, char g4, char g5)
{
	char tmp[6] = {g0, g1, g2, g3, g4, g5};
	return SixMer2Index(tmp);
}

int g::Mer2Signal::SixMer2Index(const std::vector<char>& sixmer)
{
	return SixMer2Index(&(sixmer[0]));
}

int g::Mer2Signal::SixMer2Index(const std::string& sixmer)
{
	return SixMer2Index(sixmer.c_str());
}

double g::Mer2Signal::LevelMeanAt(int index)
{
	return index_table_6mer[index][0];
// 	return index_table_6mer[index];
}

double g::Mer2Signal::LevelStdvAt(int index)
{
	return index_table_6mer[index][1];
}

double g::Mer2Signal::SdMeanAt(int index)
{
	return index_table_6mer[index][2];
}

double g::Mer2Signal::SdStdvAt(int index)
{
	return index_table_6mer[index][3];
}

double g::Mer2Signal::WeightAt(int index)
{
	return index_table_6mer[index][4];
}


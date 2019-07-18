#ifndef FIVEMER_INDEX_H__
#define FIVEMER_INDEX_H__

#include <vector>
#include <string>

#define KMER_NUM		4096
#define TABLE_ITEMS		5

namespace g{
	
class Mer2Signal
{
private:
	/* avg, var, len, #, #, #, # (defined in 5mer_index_table.rc) */
	const static double index_table_6mer[KMER_NUM][TABLE_ITEMS];
// 	const static double index_table_6mer[4096];
	
public:
	static int SixMer2Index(const std::vector<char>& sixmer);
	static int SixMer2Index(char g0, char g1, char g2, char g3, char g4, char g5);
	static int SixMer2Index(const std::string& sixmer);
	static int SixMer2Index(char const* sixmer);
	static double LevelMeanAt(int index);
	static double LevelStdvAt(int index);
	static double SdMeanAt(int index);
	static double SdStdvAt(int index);
	static double WeightAt(int index);
};
	

}


#endif
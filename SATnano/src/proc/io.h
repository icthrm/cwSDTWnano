#ifndef IO_H__
#define IO_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <6mer/6mer_index.h>
#include "hdf5.h"

namespace g{
namespace io{
	
bool ReadATCG(const char* name, std::vector<char>& genomes);

bool WriteATCG(const char* name, const std::vector<char>& genomes);

bool ReadSignalSequence(const char* name, std::vector<double>& signals);

bool WriteSignalSequence(const char* name, const std::vector<double>& signals);

void Genomes2SignalSequence(const std::vector<char>& genomes, std::vector<double>& signals, int scale = 1);

void Genomes2SignalSequence(const std::vector<char>& genomes, std::vector<double>& level_means, std::vector<double>& level_stdv, int scale = 1);

void GetGenomeStdvInPositions(const std::vector<char>& genomes, std::vector<double>& level_stdv);

std::string GetFileExtension(const std::string& filename);

std::string GetFileName(const std::string& filename);

bool GetFilesName(const char* dirname, std::vector<std::string>& filenames);

/** interface for fast5 data **/
bool ReadFast5ATCG(const char* name, std::vector<char>& genomes);

bool ReadFast5SignalSequence(const char* name, std::vector<double>& signals);

bool PrintFast5Info(const char* name);

std::string GetFast5ID(const char* name);

}

}

#endif

#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <proc/io.h>
#include "util/exception.h"
#include "6mer/6mer_index.h"

int main(int argc, char **argv)
{
	struct options opts;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}
	
	EX_TRACE("Transform genomes to signal sequence...\n");
	
	std::vector<char> genomes;
	std::vector<double> signals;
	
	g::io::ReadATCG(opts.input, genomes);
	EX_TRACE("%ld genomes are readed.\n", genomes.size());
	
	g::io::Genomes2SignalSequence(genomes, signals, 1);
	
	g::io::WriteSignalSequence("genome2sig.result", signals);
	
// 	genome::Mer2Signal::FiveMer2Index("ATAAA");
	
	return 0;
}


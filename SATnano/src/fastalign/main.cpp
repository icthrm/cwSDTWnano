#include "opts.h"
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "proc/proc.h"
#include "util/exception.h"
#include "6mer/6mer_index.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;

void WriteSequenceAlignment(const char* output, const vector<char>& genomes, const std::vector<double>& reference_orig, const std::vector<double>& peer_orig, 
							const std::vector<double>& reference, const std::vector<double>& peer, vector<pair<int,int> >& alignment)
{
	std::ofstream o(output);
	
	double diff, avgdiff = 0;
	
	for(int i = 0; i < alignment.size() && alignment[i].first+6 < genomes.size(); i++){
		o<<setw(10)<<alignment[i].first<<" "<<setw(9)<<alignment[i].second<<" | "
		<<setw(15)<<reference_orig[alignment[i].first]<<" "<<setw(15)<<peer_orig[alignment[i].second]<<" | "
		<<setw(15)<<reference[alignment[i].first]<<" "<<setw(15)<<peer[alignment[i].second]<<"          diff:"
		<<setw(15)<<(diff = std::fabs(reference[alignment[i].first]-peer[alignment[i].second]))<<"          "; avgdiff += diff;
		for(int j = 0; j < 6 && j+alignment[i].first < genomes.size(); j++){
			o<<genomes.at(alignment[i].first+j);
		}
		o<<std::endl;
	}
}

void WriteFundamentalAlignment(const char* output, double ndist, const vector<char>& genomes, const std::vector<double>& reference_orig, 
							   const std::vector<double>& peer_orig, vector<pair<int,int> >& alignment)
{
	std::ofstream o(output);
	
	o<<ndist<<std::endl;
	o<<setw(10)<<"#query_idx"<<setw(15)<<"signal_idx"<<" | "
		<<setw(15)<<"query_value"<<" "<<setw(15)<<"signal_value"<<setw(15)<<"6-mer"<<std::endl;
	
	for(int i = 0; i < alignment.size() && alignment[i].first+6 < genomes.size(); i++){
		o<<setw(10)<<alignment[i].first<<setw(15)<<alignment[i].second<<" | "
		<<setw(15)<<reference_orig[alignment[i].first]<<" "<<setw(15)<<peer_orig[alignment[i].second]<<setw(10);
		for(int j = 0; j < 6 && j+alignment[i].first < genomes.size(); j++){
			o<<genomes.at(alignment[i].first+j);
		}
		o<<std::endl;
	}
}

int main(int argc, char **argv)
{
	struct options opts;
	opts.radius = 15;
	opts.scale0 = sqrt(2);
	opts.level = 3;
	
	if(GetOpts(argc, argv, &opts) < 0){
		EX_TRACE("**WRONG INPUT!**\n")
		return 0;
	}
		
	std::vector<char> genomes;
	std::vector<double> reference;	//reference: genome
	
	EX_TIME_BEGIN("\nTransform genomes to signal sequence...");
	
	if(strcmp(g::io::GetFileExtension(opts.input).c_str(), "fast5") == 0 
		&& g::io::ReadFast5ATCG(opts.input, genomes)){
		EX_TRACE("%ld genomes are readed.\n", genomes.size());
	}
	else if(g::io::ReadATCG(opts.input, genomes)){
		EX_TRACE("%ld genomes are readed.\n", genomes.size());
	}
	else{
		EX_TRACE("Cannot open %s.\n", opts.input);
		return -1;
	}
	
	g::io::Genomes2SignalSequence(genomes, reference, 1);
	g::io::WriteSignalSequence("genome2sig.result", reference);
	
	EX_TIME_END("Transform genomes to signal sequence...");
	EX_TIME_BEGIN("Continous Wavelet Dynamic Time Warping...")
	
	std::vector<double> peer;		//peer: nanopore signal
	
	if(strcmp(g::io::GetFileExtension(opts.peer).c_str(), "fast5") == 0 
		&& !g::io::ReadFast5SignalSequence(opts.peer, peer)){
		EX_TRACE("Cannot open %s.\n", opts.peer);
		return -1;
	}
	else if(!g::io::ReadSignalSequence(opts.peer, peer)){
		EX_TRACE("Cannot open %s.\n", opts.peer);
		return -1;
	}
	else{
		EX_TRACE("%ld signals are readed.\n", peer.size());
	}
	
	std::vector<double> reference_orig(reference), peer_orig(peer);
	
	g::proc::ZScoreNormalize(reference);
	g::proc::ZScoreNormalize(peer);
	
	std::vector<std::pair<int,int> > alignment;
	double ndist = g::proc::CWDynamicTimeWarping(reference, peer, alignment, opts.scale0, opts.level, opts.radius, 1, true); //true);
	
// 	WriteSequenceAlignment(opts.output, genomes, reference_orig, peer_orig, reference, peer, alignment);
	
	if(!opts.fundamental){
		WriteSequenceAlignment(opts.output, genomes, reference_orig, peer_orig, reference, peer, alignment);
	}
	else{
		WriteFundamentalAlignment(opts.output, ndist, genomes, reference_orig, peer_orig, alignment);
	}
	
	EX_TIME_END("Continous Wavelet Dynamic Time Warping...")
	
	return 0;
}


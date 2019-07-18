#include "io.h"
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <cstring>

bool g::io::ReadATCG(const char* name, std::vector<char>& genomes)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }
    
    //-> skip first header
    std::string buf;
    if(!getline(in, buf)){
        return false;
    }
    
    while(in.good()){
        char item;
        in>>item;
        if(in.fail()){
            break;
        }
        genomes.push_back(item);
    }
    in.close();
	
    return true;
}

bool g::io::WriteATCG(const char* name, const std::vector<char>& genomes)
{
	std::ofstream out(name);
    if(!out.good()) {
        return false;
    }
    
    out<<">"<<name<<std::endl;
    
    for(size_t i = 0; i < genomes.size(); i++){
		out<<genomes[i]<<std::endl;
    }
    
    out.close();
    
    return true;
}

bool g::io::ReadSignalSequence(const char* name, std::vector<double>& signals)
{
    std::ifstream in(name);
    if(!in.good()) {
        return false;
    }

    while(in.good()){
        double item;
        in>>item;
        if(in.fail()){
            break;
        }
        signals.push_back(item);
    }
    in.close();
	
    return true;
}

bool g::io::WriteSignalSequence(const char* name, const std::vector<double>& signals)
{
	std::ofstream out(name);
    if(!out.good()) {
        return false;
    }
    
    for(size_t i = 0; i < signals.size(); i++){
		out<<signals[i]<<std::endl;
    }
    
    out.close();
    
    return true;
}

void g::io::Genomes2SignalSequence(const std::vector< char >& genomes, std::vector< double >& signals, int scale)
{
	size_t bound = genomes.size()-6;//genomes.size()%6;
	for(size_t i = 0; i < bound; i++){
		int idx = g::Mer2Signal::SixMer2Index(&(genomes[0])+i);
		double sigval = g::Mer2Signal::LevelMeanAt(idx);
		
		for(int c = scale; c--;){
			signals.push_back(sigval);
		}
	}
	
	for(size_t i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			signals.push_back(100);
		}
	}
}

void g::io::Genomes2SignalSequence(const std::vector< char >& genomes, std::vector< double >& level_means, std::vector< double >& level_stdv, int scale)
{
	size_t bound = genomes.size()-6;//genomes.size()%6;
	for(size_t i = 0; i < bound; i++){
		int idx = g::Mer2Signal::SixMer2Index(&(genomes[0])+i);
		double sigval = g::Mer2Signal::LevelMeanAt(idx);
		double sigstd = g::Mer2Signal::LevelStdvAt(idx);
		
		for(int c = scale; c--;){
			level_means.push_back(sigval);
			level_stdv.push_back(sigstd);
		}
	}
	
	for(size_t i = bound; i < genomes.size(); i++){
		for(int c = scale; c--;){
			level_means.push_back(100);
			level_stdv.push_back(0);
		}
	}
}

void g::io::GetGenomeStdvInPositions(const std::vector<char>& genomes, std::vector<double>& level_stdv)
{
	size_t bound = genomes.size()-6;//genomes.size()%6;
	for(size_t i = 0; i < bound; i++){
		int idx = g::Mer2Signal::SixMer2Index(&(genomes[0])+i);
		double sigstd = g::Mer2Signal::LevelStdvAt(idx);
		
		level_stdv.push_back(sigstd);
	}
	
	for(size_t i = bound; i < genomes.size(); i++){
		level_stdv.push_back(0);
	}
}

std::string g::io::GetFileExtension(const std::string& filename)
{
    if(filename.find_last_of(".") != std::string::npos)
        return filename.substr(filename.find_last_of(".")+1);
    return "";
}

std::string g::io::GetFileName(const std::string& filename)
{
    if(filename.find_last_of(".") != std::string::npos)
        return filename.substr(0, filename.find_last_of("."));
    return "";
}

bool g::io::GetFilesName(const char* dirname, std::vector<std::string>& filenames)
{
	DIR* dir;
    struct dirent* ptr;
    
    dir = opendir(dirname);
	if(!dir){
		return false;
	}
	
	while((ptr = readdir(dir)) != NULL){
		if(strcmp(ptr->d_name, ".") == 0){
			continue;
		}
		if(strcmp(ptr->d_name, "..") == 0){
			continue;
		}
		
        std::string name = ptr->d_name;
		filenames.push_back(name);
	}

    closedir(dir);
	
	return true;
}

//fast5 io
#include "fast5.hpp"
#include <sstream>

template < typename T >
void print_vector(std::ostream& os, const std::vector< T >& v, const std::string& delim)
{
    for (auto it = v.begin(); it != v.end(); ++it)
    {
        if (it != v.begin()) os << delim;
        os << *it;
    }
}

template < typename U, typename V >
void print_map(std::ostream& os, const std::map< U, V >& m, const std::string& prefix)
{
    for (const auto& p : m)
    {
        os << prefix << p.first << "=" << p.second << std::endl;
    }
}

bool g::io::ReadFast5ATCG(const char* name, std::vector< char >& genomes)
{
	fast5::File f;
    try{
        f.open(name);
		
		auto rs_int = f.get_raw_int_samples();
		auto fq = f.get_basecall_fastq(0);
		std::vector<std::string> tokens;
		std::string token;
		std::istringstream tokenStream(fq);
		std::getline(tokenStream, token, '\n');
		std::getline(tokenStream, token, '\n');
		genomes.resize(token.size());
		memcpy(&(genomes[0]), token.c_str(), sizeof(char)*token.size());
    }
    catch (hdf5_tools::Exception& e){
        return false;
    }
    f.close();
	
	return true;	
}

bool g::io::ReadFast5SignalSequence(const char* name, std::vector< double >& signals)
{
	fast5::File f;
    try{
        f.open(name);
		
		auto rs_int = f.get_raw_int_samples();
		for(unsigned i = 0; i < rs_int.size(); ++i){
			signals.push_back(rs_int[i]);
// 			std::cout << rs_int[i] << std::endl;
		}
    }
    catch (hdf5_tools::Exception& e){
        return false;
    }
    f.close();
	
	return true;
}


std::string g::io::GetFast5ID(const char* name)
{
	fast5::File f;
	std::string id;
	try{
		f.open(name);
		assert(f.is_open());
		auto rs_params = f.get_raw_samples_params();
		id = rs_params.read_id;
	}
	catch (hdf5_tools::Exception& e)
	{
		return "";
	}
	f.close();
	
	return id;
}

bool g::io::PrintFast5Info(const char* name)
{
	if(!fast5::File::is_valid_file(name)){
		std::cout << "not a fast5 file [" << name << "]" <<std::endl;
		return false;
	}
	fast5::File f;
	//
	// All fast5 operations are performed inside a try-catch block. This should
	// resist various hdf5 errors without leaking memory.
	//
	try
	{
		f.open(name);
		assert(f.is_open());
		//
		// extract version information for the ONT software used to generate this dataset
		//
		std::cout << "file_version=" << f.file_version() << std::endl;
		//
		// inspect channel_id params
		//
		if (f.have_channel_id_params())
		{
			auto channel_id_params = f.get_channel_id_params();
			std::cout << "channel_id/channel_number=" << channel_id_params.channel_number << std::endl
					<< "channel_id/digitisation=" << channel_id_params.digitisation << std::endl
					<< "channel_id/offset=" << channel_id_params.offset << std::endl
					<< "channel_id/range=" << channel_id_params.range << std::endl
					<< "channel_id/sampling_rate=" << channel_id_params.sampling_rate << std::endl;
		}
		//
		// inspect tracking_id params
		//
		if (f.have_tracking_id_params())
		{
			auto tracking_id_params = f.get_tracking_id_params();
			print_map(std::cout, tracking_id_params, "tracking_id/");
		}
		//
		// inspect sequences params
		//
		if (f.have_sequences_params())
		{
			auto sequences_params = f.get_sequences_params();
			print_map(std::cout, sequences_params, "sequences/");
		}
		//
		// inspect raw samples
		//
		if (f.have_raw_samples())
		{
			auto rs_params = f.get_raw_samples_params();
			auto rs = f.get_raw_samples();
			std::cout << "raw_samples/read_id=" << rs_params.read_id << std::endl
					<< "raw_samples/read_number=" << rs_params.read_number << std::endl
					<< "raw_samples/start_mux=" << rs_params.start_mux << std::endl
					<< "raw_samples/start_time=" << rs_params.start_time << std::endl
					<< "raw_samples/duration=" << rs_params.duration << std::endl
					<< "raw_samples/size=" << rs.size() << std::endl;
			const auto& e = rs.front();
			std::cout << "  (" << e << ")" << std::endl;
		}
		//
		// inspect eventdetection events
		//
		std::cout << "eventdetection_group_list=";
		print_vector(std::cout, f.get_eventdetection_group_list(), ",");
		std::cout << std::endl;
		if (f.have_eventdetection_events())
		{
			auto ed_params = f.get_eventdetection_params();
			print_map(std::cout, ed_params, "eventdetection/");
			auto ed_ev_params = f.get_eventdetection_events_params();
			auto ed_ev = f.get_eventdetection_events();
			std::cout << "eventdetection/events/abasic_found=" << ed_ev_params.abasic_found << std::endl
					<< "eventdetection/events/duration=" << ed_ev_params.duration << std::endl
					<< "eventdetection/events/median_before=" << ed_ev_params.median_before << std::endl
					<< "eventdetection/events/read_id=" << ed_ev_params.read_id << std::endl
					<< "eventdetection/events/read_number=" << ed_ev_params.read_number << std::endl
					<< "eventdetection/events/scaling_used=" << ed_ev_params.scaling_used << std::endl
					<< "eventdetection/events/start_mux=" << ed_ev_params.start_mux << std::endl
					<< "eventdetection/events/start_time=" << ed_ev_params.start_time << std::endl
					<< "eventdetection/events/size=" << ed_ev.size() << std::endl;
			const auto& e = ed_ev.front();
			std::cout << "  (mean=" << e.mean
					<< ", stdv=" << e.stdv
					<< ", start=" << e.start
					<< ", length=" << e.length << ")" << std::endl;
		} // if have_eventdetection_events
		//
		// inspect basecall groups
		//
		for (unsigned st = 0; st < 3; ++st)
		{
			std::cout << "basecall(" << st << ")/group_list=";
			print_vector(std::cout, f.get_basecall_strand_group_list(st), ",");
			std::cout << std::endl;
			// basecall sequence
			if (f.have_basecall_seq(st))
			{
				std::cout << "basecall(" << st << ")/seq_size=" << f.get_basecall_seq(st).size() << std::endl;
			}
			// basecall model
			if (f.have_basecall_model(st))
			{
				std::cout << "basecall(" << st << ")/model_file=" << f.get_basecall_model_file(st) << std::endl;
				auto m_params = f.get_basecall_model_params(st);
				auto m = f.get_basecall_model(st);
				std::cout << "basecall(" << st << ")/model/scale=" << m_params.scale << std::endl
						<< "basecall(" << st << ")/model/shift=" << m_params.shift << std::endl
						<< "basecall(" << st << ")/model/drift=" << m_params.drift << std::endl
						<< "basecall(" << st << ")/model/var=" << m_params.var << std::endl
						<< "basecall(" << st << ")/model/scale_sd=" << m_params.scale_sd << std::endl
						<< "basecall(" << st << ")/model/var_sd=" << m_params.var_sd << std::endl
						<< "basecall(" << st << ")/model/size=" << m.size() << std::endl;
				const auto& e = m.front();
				std::cout << "  (kmer=" << e.get_kmer()
						<< ", level_mean=" << e.level_mean
						<< ", level_stdv=" << e.level_stdv << ")" << std::endl;
			}
			// basecall events
			if (f.have_basecall_events(st))
			{
				auto ev = f.get_basecall_events(st);
				std::cout << "basecall(" << st << ")/events/size=" << ev.size() << std::endl;
				const auto& e = ev.front();
				std::cout << "  (mean=" << e.mean
						<< ", stdv=" << e.stdv
						<< ", start=" << e.start
						<< ", length=" << e.length
						<< ", model_state=" << e.get_model_state()
						<< ", p_model_state=" << e.p_model_state
						<< ", move=" << e.move << ")" << std::endl;
			}
			// basecall alignment
			if (st == 2 and f.have_basecall_alignment())
			{
				auto al = f.get_basecall_alignment();
				std::cout << "basecall(2)/alignment/size=" << al.size() << std::endl;
				const auto& e = al.front();
				std::cout << "  (template_index=" << e.template_index
						<< ", complement_index=" << e.complement_index
						<< ", kmer=" << e.get_kmer() << ")" << std::endl;
			}
		} // for st
	}
	catch (hdf5_tools::Exception& e)
	{
		std::cout << "hdf5 error: " << e.what() << std::endl;
	}
	
	return true;
}


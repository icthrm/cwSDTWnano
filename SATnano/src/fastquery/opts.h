#ifndef OPTS_H__
#define OPTS_H__

#include <iostream>
#include <sstream>
#include <cassert>
extern "C" {
#include <getopt.h>
}
#include "util/exception.h"

struct options {
    char input[255];
    char peer[255];
	char output[255];
	int radius;
	int mode;
	float scale_diff;
	int K;
	int L;
	int fundamental;
};

inline int GetOpts(int argc, char **argv, options* opts_) {

    static struct option longopts[] = {
        { "help",            no_argument,            NULL,              'h' },
        { "input",     required_argument,      NULL,    			  'i' },
        { "peer",    	     required_argument,      NULL,      	'p' },
        { "output",      	required_argument,      NULL,              'o' },
        { "radius",      required_argument,      NULL,              'r' },
		{ "mode",      required_argument,      NULL,              'm' },
		{ "scale",      required_argument,      NULL,              's' },
		{ "K",      required_argument,      NULL,              'K' },
		{ "L",      required_argument,      NULL,              'L' },
		{ "fundamental",      no_argument,      NULL,              'f' },
        { NULL,              0,                      NULL,               0 }
    };

    if((argc != 7 && argc != 8 && argc != 9 && argc != 10 && argc != 11 && argc != 12 && argc != 13 && argc != 14 &&  argc != 15 && argc != 16 && argc != 17 && argc != 18) 
		&& argc >= 3 || (argc == 2 && argv[1][0] != '-' && argv[1][1] != 'h') || argc == 1) {
        EX_TRACE("[-i INPUT (GENOME) SEQUENCE][-p PEER SIGNAL][-o OUTPUT]([-r NEIGHBOUR RADIUS])([-s SCALE DIFFERENCE])([-K SEED NUMBER])([-L SEED LENGTH])([-m MODE])([-f FUNDAMENTAL OUTPUT])\n");
        return -1;
    }

    int ch;
    while((ch = getopt_long(argc, argv, "fhi:p:o:r:s:m:K:L:", longopts, NULL))!= -1) {
        switch(ch) {

        case '?':
            EX_TRACE("Invalid option '%s'.", argv[optind-1]);
            return -1;

        case ':':
            EX_TRACE("Missing option argument for '%s'.", argv[optind-1]);
            return -1;

        case 'h':
            EX_TRACE("[-i INPUT (GENOME) SEQUENCE][-p PEER SIGNAL][-o OUTPUT]([-r NEIGHBOUR RADIUS])([-s SCALE DIFFERENCE])([-K SEED NUMBER])([-L SEED LENGTH])([-m MODE])([-f FUNDAMENTAL OUTPUT])\n"
                     "INPUT (GENOME) SEQUENCE: reference (inquiry) sequence i.e. ATCG...;\n"
					 "PEER SIGNAL: (nanopore) signal that will be transformed to align with reference;\n"
					 "OUTPUT: signal alignment result;\n"
					 "NEIGHBOUR RADIUS: warp search radius (default 50);\n"
					 "MODE: 0 for direct mode; 1 for fast mode (default 0);\n"
					 "SCALE DIFFERENCE: estimated scale diffference between genome sequence and raw signal (default 8);\n"
					 "SEED NUMBER: the number of seeds; works in fast mode (default 3);\n"
					 "SEED LENGTH: the length of the short seed; works in fast mode (default 128);\n"
					 "FUNDAMENTAL OUTPUT: set the output in the format for non-standard signal detection (no param).\n"
					 "\n#WARNING!!! SAT now accepts fast5 files as input: \n"
					 "    for \"-i xxx.fast5\", the ATCG sequence will only be read;\n"
					 "    for \"-p xxx.fast5\", the raw signal will only be read.\n"
					 );
						
			exit(0);
            return 0;

        case 'i':
        {
            std::istringstream iss(optarg);
            iss >> opts_->input;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
		case 'K':
        {
            std::istringstream iss(optarg);
            iss >> opts_->K;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;	

		case 'L':
        {
            std::istringstream iss(optarg);
            iss >> opts_->L;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
        case 'p':
        {
            std::istringstream iss(optarg);
            iss >> opts_->peer;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;

        case 'o':
        {
            std::istringstream iss(optarg);
            iss >> opts_->output;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
		case 'r':
        {
            std::istringstream iss(optarg);
            iss >> opts_->radius;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
		case 'm':
        {
            std::istringstream iss(optarg);
            iss >> opts_->mode;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
		case 's':
        {
            std::istringstream iss(optarg);
            iss >> opts_->scale_diff;
            if(iss.fail())
                EX_TRACE("Invalid argument '%s'.", optarg);
        }
        break;
		
		case 'f':
		{
            opts_->fundamental = 1;		
		}
		break;

        case 0:
            break;

        default:
            assert(false);
        }
    }
    return 1;
}

#endif
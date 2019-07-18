#ifndef CWDTW_SOURCE
#define CWDTW_SOURCE

#include "proc.h"
#include "wavelib.h"
#include <malloc.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <float.h>
#include <limits.h>
#include <util/exception.h>
#include <util/qsort.h>
#include <numeric>
#include <algorithm>
#include <map>

void g::proc::BoundGeneration(std::vector<std::pair<int,int> >& cosali, int radius, std::vector<std::pair<int,int> >& bound)
{
	bool firstorder = true;
	
	if(cosali[cosali.size()-1].first > cosali[cosali.size()-1].second){
		firstorder = false;
	}
	
	if(!firstorder){		//genome first
		for(int i = cosali.size(); i--;){
			std::swap(cosali[i].first, cosali[i].second);
		}
	}
	
	std::vector<std::pair<int,int> > cosali2;
	
	cosali2.push_back(cosali[0]);
	for(int i = 1; i < cosali.size(); i++){
		if(cosali[i].first != cosali2[cosali2.size()-1].first){
			cosali2.push_back(cosali[i]);
		}
		else{
			cosali2[cosali2.size()-1].second = cosali[i].second;
		}
	}
	
	cosali = cosali2;
	
	bound.resize(cosali[cosali.size()-1].first+1);
	bound.assign(cosali[cosali.size()-1].first+1, std::make_pair(-1,-1));
	
	for(int i = 1; i < cosali.size(); i++){
		int pre_idx = cosali[i-1].first, cur_idx = cosali[i].first;
		int pre_anchor = cosali[i-1].second, cur_anchor = cosali[i].second;
		
		double anchor_diff = (cur_anchor-pre_anchor)/double(cur_idx-pre_idx);
		
		int neighbor = int(anchor_diff+0.5)+radius;
		
		for(int k = pre_idx, count = 0; k < cur_idx; k++, count++){
			int mid = pre_anchor+std::round(count*anchor_diff);		//assign point relationship
			bound[k].first = mid-neighbor < 0 ? 0 : mid-neighbor;
			bound[k].second = mid+neighbor > cosali[cosali.size()-1].second ? cosali[cosali.size()-1].second : mid+neighbor;
		}
	}
	
	bound[0].first = 0;
	bound[bound.size()-1].first = bound[bound.size()-2].first;
	bound[bound.size()-1].second = cosali[cosali.size()-1].second;
}

void g::proc::SquareBoundGeneration(std::vector<std::pair<int,int> >& cosali, int radius, std::vector<std::pair<int,int> >& bound)
{
	if(radius > std::min(cosali[cosali.size()-1].first, cosali[cosali.size()-1].second)){
		radius = int(0.25*std::min(cosali[cosali.size()-1].first, cosali[cosali.size()-1].second));
		if(radius == 0){
			radius = 1;
		}
	}
	
	bool firstorder = true;
	
	if(cosali[cosali.size()-1].first > cosali[cosali.size()-1].second){
		firstorder = false;
	}
	
	if(!firstorder){		//genome first
		for(int i = cosali.size(); i--;){
			std::swap(cosali[i].first, cosali[i].second);
		}
	}
	
	std::vector<std::pair<int,int> > cosali2;
	
	cosali2.push_back(cosali[0]);
	for(size_t i = 1; i < cosali.size(); i++){
		if(cosali[i].first != cosali2[cosali2.size()-1].first){
			cosali2.push_back(cosali[i]);
		}
		else{
			cosali2[cosali2.size()-1].second = cosali[i].second;
		}
	}
	
	cosali = cosali2;
	cosali2.resize(cosali[cosali.size()-1].first+1);
	cosali2[0].first = 0; cosali2[0].second = 0;
	for(size_t i = 1; i < cosali.size(); i++){
		int pre_idx = cosali[i-1].first, cur_idx = cosali[i].first;
		int pre_anchor = cosali[i-1].second, cur_anchor = cosali[i].second;
		double anchor_diff = (cur_anchor-pre_anchor)/double(cur_idx-pre_idx);
		
		for(int k = pre_idx, count = 0; k < cur_idx; k++, count++){
			int mid = pre_anchor+int(count*anchor_diff);		//assign point relationship
			cosali2[k].first = k;
			cosali2[k].second = mid;
		}
	}
	
	int rboundary = cosali[cosali.size()-1].second;
	
	cosali2[cosali2.size()-1].first = cosali[cosali.size()-1].first; cosali2[cosali2.size()-1].second = rboundary;
	
	bound.resize(cosali[cosali.size()-1].first+1, std::make_pair(INT_MAX, INT_MIN));
// 	bound.assign(cosali[cosali.size()-1].first+1, std::make_pair(-1,-1));
	bound[0].first = 0; bound[0].second = 0;
	
// 	std::cout<<cosali.size()<<std::endl;
	
	for(size_t k = 1; k < cosali2.size(); k++){
		double lf = std::min(cosali2[k].second-radius, bound[k].first);
		bound[k].first = lf < 0 ? 0 : lf;
		double rf = std::max(cosali2[k].second+radius, bound[k].second);
		bound[k].second = rf > rboundary ? rboundary : rf;
		
		size_t tidx = std::min(size_t(cosali2[k].first+radius), bound.size()-1);
		for(size_t j = std::max(cosali2[k].first-radius, 0); j < tidx; j++){
			bound[j].first = std::min(cosali2[k].second-radius, bound[j].first);
			if(bound[j].first < 0){
				bound[j].first = 0;
			}
			bound[j].second = std::max(cosali2[k].second+radius, bound[j].second);
			if(bound[j].second > rboundary){
				bound[j].second = rboundary;
			}
		}
	}
	
	bound[bound.size()-1].first = bound[bound.size()-2].first;
	bound[bound.size()-1].second = rboundary;
}

/** if @boundmode == 0, use rectangle boundary; else == 1, use square boundary*/
static void CoarseAlignment(std::vector<std::vector<double> >& sig1, std::vector<std::vector<double> >& sig2, 
							std::vector<std::pair<int,int> >& alignment, double* totaldiff = 0, int radius = 180, int boundmode = 0)
{
	double tdiff;
	std::vector<std::pair<int, double> > sig1peaks, sig2peaks;
	
	for(int k = 0; k < sig1.size(); k++){
		g::proc::PeakPick(sig1[k], sig1peaks);
		g::proc::PeakPick(sig2[k], sig2peaks);
// 		std::cout<<sig1peaks.size()<<"\t"<<sig2peaks.size()<<std::endl;
		std::vector<double> peak1(sig1peaks.size());
		std::vector<double> peak2(sig2peaks.size());
		for(int i = sig1peaks.size(); i--;){
			peak1[i] = sig1peaks[i].second;
		}
		for(int i = sig2peaks.size(); i--;){
			peak2[i] = sig2peaks[i].second;
		}
		if(k == 0){
			tdiff = g::proc::DynamicTimeWarping(peak1, peak2, alignment);
		}
		else{
			int c = 0;
			for(int i = 0; i < alignment.size(); i++){
				while(sig1peaks[c].first < alignment[i].first){
					c++;
				}
				alignment[i].first = c;	
			}
			
			c = 0;
			for(int i = 0; i < alignment.size(); i++){
				while(sig2peaks[c].first < alignment[i].second){
					c++;
				}
				alignment[i].second = c;	
			}
			
			std::vector<std::pair<int,int> > bound; 
			if(boundmode == 0){
				g::proc::BoundGeneration(alignment, radius, bound);
			}
			else{
				g::proc::SquareBoundGeneration(alignment, radius, bound);
			}
			tdiff = g::proc::BoundDynamicTimeWarping(peak1, peak2, bound, alignment);
		}
		for(int i = alignment.size(); i--;){
			alignment[i].first = sig1peaks[alignment[i].first].first;
			alignment[i].second = sig2peaks[alignment[i].second].first;
		}
	}
	
	if(totaldiff){
		*totaldiff = tdiff;
	}
}

#undef EX_TRACE
#define EX_TRACE(MSG,...)		//unmark trace

double g::proc::CWDynamicTimeWarping(const std::vector< double >& reference, const std::vector< double >& peer, std::vector< std::pair< int, int > >& alignment, 
									 double scale0, double l, double radius, int boundmode, bool left_constrained)
{
	double scale = (double)peer.size()/reference.size();
	
	std::vector<std::vector<double> > rcwt, pcwt;
	
	EX_TRACE("CWT Analysis...\n")
	
	int npyr = l;  //1
	double dscale = 1;
	
	g::proc::CWTAnalysis(peer, pcwt, scale0*scale, dscale, npyr);
	g::proc::CWTAnalysis(reference, rcwt, scale0, dscale, npyr);
	
	// 	g::io::WriteSignalSequence("rcwt", rcwt[0]);
	// 	g::io::WriteSignalSequence("pcwt", pcwt[0]);
	
	//if multiscale is used, pyr logical should be added.
	for(int i = 0; i < rcwt.size(); i++){
		g::proc::ZScoreNormalize(rcwt[i]);
		g::proc::ZScoreNormalize(pcwt[i]);
	}
	std::vector<std::pair<int,int> > cosali;
	double tdiff;
	
	EX_TRACE("Coarse Alignment...\n")
	CoarseAlignment(rcwt, pcwt, cosali, &tdiff, radius*2, boundmode);
	EX_TRACE("Average Deviation (%.1lf/%ld=%.3lf)\n", tdiff, cosali.size(), tdiff/cosali.size())
	
	std::vector<std::pair<int,int> > bound;
	
	BoundGeneration(cosali, radius*scale, bound);
	
	EX_TRACE("Fine Alignment...\n")
	
	if(left_constrained){
		tdiff = g::proc::BoundDynamicTimeWarpingR(reference, peer, bound, alignment);
	}
	else{
		tdiff = g::proc::BoundDynamicTimeWarping(reference, peer, bound, alignment);
	}

	EX_TRACE("Average Deviation (%.1lf/%ld=%.6lf)\n", tdiff, alignment.size(), tdiff/alignment.size())
	
	return tdiff/alignment.size();
}


//I = internal function
static void SubsequenceDTW(const std::vector< double >& seq1, const std::vector< double >& seq2, double** score)
{
	for(int i = 0; i < seq1.size(); i++){
		for(int j = 0; j < seq2.size(); j++){
			score[i][j] = std::fabs(seq1[i]-seq2[j]);
		}
	}

	for(int i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}

	for(int i = 1; i < seq1.size(); i++){
		for(int j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(std::min(score[i-1][j], score[i][j-1]), score[i-1][j-1]);
		}
	}
}

static void PathBackTrack(double** score, int xi, int yj, std::vector<std::pair<int,int> >& alignment)
{
		alignment.clear();
		int i = xi, j = yj;
		while(true){
			alignment.push_back(std::make_pair(i,j));
			int ipre = i-1 < 0 ? 0 : i-1;
			int jpre = j-1 < 0 ? 0 : j-1;
			
			double premin = std::min(std::min(score[ipre][j], score[i][jpre]), score[ipre][jpre]);
			
			if(premin == score[ipre][jpre]){
				i = ipre; j = jpre;
			}
			if(premin == score[ipre][j]){
				i = ipre;
			}
			if(premin == score[i][jpre]){
				j = jpre;
			}
	//		if(i == 0 /&& j == 0){
			if( i == 0 ){
				if(j-1 > 0){
					alignment.push_back(std::make_pair(i,j-1));
				}
				else{
					alignment.push_back(std::make_pair(i,j));
				}
				break;
			}
		}
		std::reverse(alignment.begin(), alignment.end());
}

/** if @boundmode == 0, use rectangle boundary; else == 1, use square boundary*/
static void CoarsePathRefine(std::vector<std::vector<double> >& sig1, std::vector<std::vector<double> >& sig2, 
							std::vector<std::pair<int,int> >& alignment, double* totaldiff = 0, int radius = 128, int boundmode = 0)
{
	double tdiff;
	std::vector<std::pair<int, double> > sig1peaks, sig2peaks;
	
	for(int k = 0; k < sig1.size(); k++){
		g::proc::PeakPick(sig1[k], sig1peaks);
		g::proc::PeakPick(sig2[k], sig2peaks);
// 		std::cout<<sig1peaks.size()<<"\t"<<sig2peaks.size()<<std::endl;
		std::vector<double> peak1(sig1peaks.size());
		std::vector<double> peak2(sig2peaks.size());
		for(int i = sig1peaks.size(); i--;){
			peak1[i] = sig1peaks[i].second;
		}
		for(int i = sig2peaks.size(); i--;){
			peak2[i] = sig2peaks[i].second;
		}

		int c = 0;
		for(int i = 0; i < alignment.size(); i++){
			while(sig1peaks[c].first < alignment[i].first){
				c++;
			}
			alignment[i].first = c;	
		}
		
		c = 0;
		for(int i = 0; i < alignment.size(); i++){
			while(sig2peaks[c].first < alignment[i].second){
				c++;
			}
			alignment[i].second = c;	
		}
		
		std::vector<std::pair<int,int> > bound; 
		if(boundmode == 0){
			g::proc::BoundGeneration(alignment, radius, bound);
		}
		else{
			g::proc::SquareBoundGeneration(alignment, radius, bound);
		}
		tdiff = g::proc::BoundDynamicTimeWarping(peak1, peak2, bound, alignment);
		
		for(int i = alignment.size(); i--;){
			alignment[i].first = sig1peaks[alignment[i].first].first;
			alignment[i].second = sig2peaks[alignment[i].second].first;
		}
	}
	
	if(totaldiff){
		*totaldiff = tdiff;
	}
}

static void SubsequenceCoarsePathes(const std::vector<double>& reference, const std::vector<double>& operated_peer, int level, std::vector<std::vector<std::pair<int,int> > >& path_v)
{
	std::vector<std::vector<double> > rcwt, pcwt;
	
	g::proc::CWTAnalysis(reference, rcwt, pow(2, level));
	g::proc::CWTAnalysis(operated_peer, pcwt, pow(2, level));
	
	std::vector<std::pair<int, double> > sig1peaks, sig2peaks;
	g::proc::PeakPick(rcwt[0], sig1peaks);
	g::proc::PeakPick(pcwt[0], sig2peaks);

	std::vector<double> peak1(sig1peaks.size()), peak2(sig2peaks.size());
	
	for(int i = sig1peaks.size(); i--;){
		peak1[i] = sig1peaks[i].second;
	}
	for(int i = sig2peaks.size(); i--;){
		peak2[i] = sig2peaks[i].second;
	}
	
	int XL = peak1.size(), YL = peak2.size();
	
	double* score[XL];
	for(int i = 0; i < XL; i++){
		score[i] = new double[YL];
	}
	
	SubsequenceDTW(peak1, peak2, score);
	double thre = std::accumulate(score[XL-1]+peak1.size()-2, score[XL-1]+YL-1, 0.0f)*.5;  //median_copy(YL, score[XL-1]); //
	
	for(int j = 0; j < YL; j++){
		if(score[XL-1][j] < thre){
			std::vector<std::pair<int,int> > alignment;
			PathBackTrack(score, XL-1, j, alignment);
			path_v.push_back(alignment);
		}
	}
	
	for(int i = 0; i < XL; i++){
		delete score[i];
	}
	
	for(int k = 0; k < path_v.size(); k++){
		std::vector<std::pair<int,int> >& alignment = path_v[k];
		
		for(int i = alignment.size(); i--;){
			alignment[i].first = sig1peaks[alignment[i].first].first;
			alignment[i].second = sig2peaks[alignment[i].second].first;
		}
	}	
}

static bool BasePathVerify(const std::vector<std::vector<std::pair<int,int> > >& path_v, 
						   const std::vector<std::vector<std::pair<int,int> > >& seed_align, std::vector<std::pair<int,int> >& cosali)
{
	int min_dist = INT_MAX;
	int min_index = -1;
	
	std::vector<int> min_path_idx;
	
	for(int i = 0; i < path_v.size(); i++){
		int dist = 0;
		for(int j = 0; j < seed_align.size(); j++){
			if(seed_align[j][0].second < path_v[i][0].second){
				dist += path_v[i][0].second-seed_align[j][0].second;
			}
			if(seed_align[j][seed_align[j].size()-1].second < path_v[i][0].second){
				dist += path_v[i][0].second-seed_align[j][seed_align[j].size()-1].second;
			}
			if(seed_align[j][seed_align[j].size()-1].second > path_v[i][path_v[i].size()-1].second){
				dist += seed_align[j][seed_align[j].size()-1].second - path_v[i][path_v[i].size()-1].second;
			}
			if(seed_align[j][0].second > path_v[i][path_v[i].size()-1].second){
				dist += seed_align[j][0].second - path_v[i][path_v[i].size()-1].second;
			}
		}
		if(dist < min_dist){
			min_dist = dist;
			min_index = i;
		}
		
		if(min_dist == 0){
			min_path_idx.push_back(min_index);
		}
	}
	
	if(min_path_idx.size() == 0){
		return false;			//in this case, there may exits false mapping.
// 		min_path_idx.push_back(min_index);
	}
	
	min_dist = INT_MAX;
	min_index = -1;
	
	for(int i = 0; i < min_path_idx.size(); i++){
		int dist = path_v[min_path_idx[i]][path_v[min_path_idx[i]].size()-1].second-path_v[min_path_idx[i]][path_v[min_path_idx[i]].size()-1].second;
		if(min_dist > dist){
			min_index = i;
		}
	}
	
	if(min_dist > (seed_align[seed_align.size()-1][seed_align[seed_align.size()-1].size()-1].second-seed_align[0][0].second)+45){		//the path need reestimation
		return false;
	}
	
	cosali = path_v[min_path_idx[min_index]];	//this is the path for original data index

	{		//this section try to remap the path with the seed
	for(std::vector<std::pair<int,int> >::iterator iter = cosali.begin(); iter != cosali.end();){
		if(iter->second < seed_align[0][0].second){
			iter = cosali.erase(iter);
		}
		else{
			iter++;
		}
	}
	
	if(cosali[0].second > seed_align[0][0].second){
		cosali.insert(cosali.begin(), std::make_pair(0, seed_align[0][0].second));
	}
	
	int ranchor_idx = seed_align[seed_align.size()-1][seed_align[seed_align.size()-1].size()-1].second;
	for(int i = 0; i < cosali.size(); i++){
		if(cosali[i].second > ranchor_idx){
			cosali[i].second = ranchor_idx;
		}
	}
	
	if(cosali[cosali.size()-1].second > seed_align[seed_align.size()-1][seed_align[seed_align.size()-1].size()-1].second){
		cosali.insert(cosali.end(), std::make_pair(0, seed_align[seed_align.size()-1][seed_align[seed_align.size()-1].size()-1].second));
	}	}
	return true;
}

static void BaseCoarsePath(const std::vector<double>& reference, const std::vector<double>& operated_peer, int start, int end, int level, std::vector<std::pair<int,int> >& alignment)
{
	std::vector<double> segment(end-start+1);
	memcpy(&(segment[0]), &(operated_peer[start]), sizeof(double)*segment.size());
	
	std::vector<std::vector<double> > rcwt, pcwt;
	
	g::proc::CWTAnalysis(reference, rcwt, pow(2, level));
	g::proc::CWTAnalysis(segment, pcwt, pow(2, level));
	
	std::vector<std::pair<int, double> > sig1peaks, sig2peaks;
	g::proc::PeakPick(rcwt[0], sig1peaks);
	g::proc::PeakPick(pcwt[0], sig2peaks);

	std::vector<double> peak1(sig1peaks.size()), peak2(sig2peaks.size());
	
	for(int i = sig1peaks.size(); i--;){
		peak1[i] = sig1peaks[i].second;
	}
	for(int i = sig2peaks.size(); i--;){
		peak2[i] = sig2peaks[i].second;
	}
	
	g::proc::DynamicTimeWarping(peak1, peak2, alignment);
	
	for(int i = alignment.size(); i--;){
		alignment[i].first = sig1peaks[alignment[i].first].first;
		alignment[i].second = sig2peaks[alignment[i].second].first+start;
	}	
}

double g::proc::DirectInquiry(const std::vector<double>& reference, const std::vector<double>& peer, std::vector<std::pair<int,int> >& alignment, double scale_diff, 
						double radius, int boundmode, bool left_constrained)
{
	//compress the nanopore raw signal
	std::vector<double> operated_peer;
	g::proc::Resample(1, scale_diff, peer, operated_peer);
	
	std::vector<double> local_score;
	std::vector<std::pair<int,int> > cosali;
	double tdiff = g::proc::SubsequenceDynamicTimeWarping(reference, operated_peer, local_score, &cosali);
	
	int anchor_idx = cosali[0].second;
	for(int i = 0; i < cosali.size(); i++){
		cosali[i].second -= anchor_idx;
		cosali[i].second *= scale_diff;
	}
	
	anchor_idx *= scale_diff;
	
	operated_peer.resize(cosali[cosali.size()-1].second+1);
	memcpy(&(operated_peer[0]), &(peer[anchor_idx]), sizeof(double)*operated_peer.size());
	
	std::vector<std::pair<int,int> > bound;
	
	BoundGeneration(cosali, radius*scale_diff, bound);
	
	EX_TRACE("Fine Alignment...\n")
	
	if(left_constrained){
		tdiff = g::proc::BoundDynamicTimeWarpingR(reference, operated_peer, bound, alignment);
	}
	else{
		tdiff = g::proc::BoundDynamicTimeWarping(reference, operated_peer, bound, alignment);
	}
	
	for(int i = 0; i < alignment.size(); i++){			//map to a segmentation in original data
		alignment[i].second += anchor_idx;
	}
	
	EX_TRACE("Average Deviation (%.1lf/%ld=%.6lf)\n", tdiff, alignment.size(), tdiff/alignment.size())
	
	return tdiff/alignment.size();
}
//base_function.cpp


// double data[2][2] = { //    X      Y    rows = 2; Y = a + bx
//     {187.1, 25.4},
//     {179.5, 22.8},
// }
static int LinearRegression(double *data, int rows, double* a, double* b, double* SquarePoor)
{
    int m;
    double *p, Lxx = 0.0, Lxy = 0.0, xa = 0.0, ya = 0.0;
    if(data == 0 || a == 0 || b == 0 || rows < 1){
        return -1;
	}
	
    for(p = data, m = 0; m < rows; m ++){
        xa += *p ++;
        ya += *p ++;
    }
    xa /= rows;                                     // X mean
    ya /= rows;                                     // y mean
    for(p = data, m = 0; m < rows; m ++, p += 2){
        Lxx += ((*p - xa) * (*p - xa));             // Lxx = Sum((X - Xa))
        Lxy += ((*p - xa) * (*(p + 1) - ya));       // Lxy = Sum((X - Xa)(Y - Ya))
    }
    *b = Lxy / Lxx;                                 // b = Lxy / Lxx
    *a = ya - *b * xa;                              // a = Ya - b*Xa
    if(SquarePoor == 0){
        return 0;
	}
    // deviation
    SquarePoor[0] = SquarePoor[1] = 0.0;
    for(p = data, m = 0; m < rows; m++, p++){
        Lxy = *a + *b * *p ++;
        SquarePoor[0] += ((Lxy - ya) * (Lxy - ya)); 
        SquarePoor[1] += ((*p - Lxy) * (*p - Lxy)); 
    }
    SquarePoor[2] = SquarePoor[0];                  
    SquarePoor[3] = SquarePoor[1] / (rows - 2);     
    return 0;
}

static bool LegalityCheck(std::vector<std::vector<std::pair<int,int> > >& seed_align)
{
	if(seed_align.size() == 1){
		return true;
	}
	for(int i = 1; i < seed_align.size(); i++){
		if(seed_align[i].back().second < seed_align[i-1].back().second || seed_align[i][0].second < seed_align[i-1][0].second){
			return false;
		}
	}
	
	int size = 0;
	for(int i = 0; i < seed_align.size(); i++){
		size += seed_align[i].size();
	}
	
	double data[size][2];
	int c = 0; 
	for(int i = 0; i < seed_align.size(); i++){
		for(int j = 0; j < seed_align[i].size(); j++){
			data[c][0] = seed_align[i][j].first;
			data[c][1] = seed_align[i][j].second;
			c++;
		}
	}
	double a, b; double squarepoor[4];
	
	LinearRegression( (double*)data, size, &a, &b, squarepoor);
	if(b < 0.8 || b > 1.2){			//linear feature
		return false;
	}
	if(sqrt(squarepoor[0] / (squarepoor[0] + squarepoor[1])) < 0.85){		//statistic feature 
		return false;
	}
// 	std::cout<<b<<" "<<sqrt(squarepoor[0] / (squarepoor[0] + squarepoor[1]))<<std::endl;
	
	return true;
}

double g::proc::FastInquiry(const std::vector<double>& reference, const std::vector<double>& peer, std::vector<std::pair<int,int> >& alignment, double scale_diff, 
							int seed_num, int seed_length, double radius, int boundmode, bool left_constrained)
{
	//compress the nanopore raw signal
	std::vector<double> operated_peer;
	g::proc::Resample(1, scale_diff, peer, operated_peer);
	
	
	if(seed_num <= 0){
		seed_num = 1;
	}
	
	int ankle[seed_num];
	
	ankle[0] = 0; ankle[seed_num-1] = reference.size()-1-seed_length;
	
	int seed_step = (reference.size()-seed_length)/seed_num;
	for(int i = 1; i < seed_num-1; i++){
		ankle[i] = seed_step*i;
	}
	
	std::vector<std::vector<std::pair<int,int> > > seed_align;
	
	//get seed alignment
	for(int i = 0; i < seed_num; i++){
		std::vector<double> test(seed_length);
		memcpy(&(test[0]), &(reference[ankle[i]]), sizeof(double)*seed_length);
		
		std::vector<double> local_score;
		std::vector<std::pair<int,int> > alignment;
		double v = g::proc::SubsequenceDynamicTimeWarping(test, operated_peer, local_score, &alignment);
		seed_align.push_back(alignment);
	}
	
	for(int k = 0; k < seed_align.size(); k++){
		std::vector<std::pair<int,int> >& alignment = seed_align[k];
		
		for(int i = alignment.size(); i--;){
			alignment[i].first += ankle[k];
		}
	}
	
// 	for(int i = 0; i < seed_align.size(); i++){
// 		std::cout<<seed_align[i][0].second<<" "<<seed_align[i][seed_align[i].size()-1].second<<std::endl;
// 	}
	
	if(!LegalityCheck(seed_align)){
		alignment.clear(); alignment.push_back(std::make_pair(-1, -1));
		return DBL_MAX;
	}
	
	//get coarse alignment
	int MAX_LEVEL = std::log(reference.size())-2;
	
	std::vector<std::vector<std::pair<int,int> > > path_v;
	SubsequenceCoarsePathes(reference, operated_peer, MAX_LEVEL, path_v);
	
	std::vector<std::pair<int,int> > cosali;
	if(!BasePathVerify(path_v, seed_align, cosali)){
		BaseCoarsePath(reference, operated_peer, seed_align.front().front().second, seed_align.back().back().second, MAX_LEVEL, cosali);
	}
	
// 	for(int i = 0; i < cosali.size(); i++){			//map to a segmentation in original data
// 		std::cout<<cosali[i].first<<"\t"<<cosali[i].second<<std::endl;
// 	}
	
	int anchor_idx = cosali[0].second;
	for(int i = 0; i < cosali.size(); i++){			//map to a segmentation in original data
		cosali[i].second -= anchor_idx;
		cosali[i].second *= scale_diff;
	}
	anchor_idx *= scale_diff;
	
	double dscale = 1;
	operated_peer.resize(cosali[cosali.size()-1].second+1);
	memcpy(&(operated_peer[0]), &(peer[anchor_idx]), sizeof(double)*operated_peer.size());
	
	std::vector<std::vector<double> > rcwt, pcwt;
	g::proc::CWTAnalysis(reference, rcwt, 2, dscale, MAX_LEVEL-1);
	g::proc::CWTAnalysis(operated_peer, pcwt, 2*scale_diff, dscale, MAX_LEVEL-1);
	
	//if multiscale is used, pyr logical should be added.
	for(int i = 0; i < rcwt.size(); i++){
		g::proc::ZScoreNormalize(rcwt[i]);
		g::proc::ZScoreNormalize(pcwt[i]);
	}

	double tdiff;
	
	CoarsePathRefine(rcwt, pcwt, cosali, &tdiff, radius*2, boundmode);
	
	std::vector<std::pair<int,int> > bound;
	
	BoundGeneration(cosali, radius*scale_diff, bound);
	
	EX_TRACE("Fine Alignment...\n")
	
	if(left_constrained){
		tdiff = g::proc::BoundDynamicTimeWarpingR(reference, operated_peer, bound, alignment);
	}
	else{
		tdiff = g::proc::BoundDynamicTimeWarping(reference, operated_peer, bound, alignment);
	}
	
	for(int i = 0; i < alignment.size(); i++){			//map to a segmentation in original data
		alignment[i].second += anchor_idx;
	}

	EX_TRACE("Average Deviation (%.1lf/%ld=%.6lf)\n", tdiff, alignment.size(), tdiff/alignment.size())
	
	return tdiff/alignment.size();
}

//K=3 以及 L=192 r=70 默认参数

#endif
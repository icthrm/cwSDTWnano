#ifndef PROC_H__
#define PROC_H__

#include "io.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unistd.h>
// #include <fftw3.h>
#include "wavelib.h"

namespace g{
namespace proc{
// void SpectrumAmplitude(std::vector<double>& signals, std::vector<double>& amplitude);
// 
// void SpectrumAmplitude(double* signals, int size, std::vector<double>& amplitude);

void ZScoreNormalize(std::vector<double>& signals, double* avg = NULL, double* stdev = NULL);

void Resample(int upfactor, int downfactor, const std::vector<double>& input, std::vector<double>& output);

void ScaleSignalAmplitude(std::vector<double>& signals, double scale = 1);

void MedianFilter(std::vector<double>& signals, int width = 5);

void WaveDenoise(std::vector<double>& signals, bool soft = true);

void RemoveHotSpot(std::vector<double>& signals, int thre = 5);

void LaplaceDiff(const std::vector<double>& raw, std::vector<double>& ldiff);

void Diff(const std::vector<double>& raw, std::vector<double>& diff);

double L2Distance(const std::vector<double>& sig1, const std::vector<double>& sig2);
double L1Distance(const std::vector<double>& sig1, const std::vector<double>& sig2);

void PeakPick(const std::vector<double>& raw, std::vector<std::pair<int, double> >& peaks);

void LocalMaximum(const std::vector<double>& raw, std::vector<std::pair<int, double> >& peaks);

/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void CWTAnalysis(const std::vector<double>& raw, std::vector<std::vector<double> >& output, double scale0, double dscale = 1, int npyr = 1);

/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void CWTAnalysisEsp(const std::vector<double>& raw, std::vector<std::vector<double> >& output, double scale0, double dscale = 1, int npyr = 1);

void DumpCWTSpectrogram(const std::vector<std::vector<double> >& waves, char* filename);

double DynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment);

double DynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment);

double SubsequenceDynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, 
									 std::vector<double>& local_score, std::vector<std::pair<int,int> >* alignment = NULL);

void SquareBoundGeneration(std::vector<std::pair<int,int> >& cosali, int radius, std::vector<std::pair<int,int> >& bound);

void BoundGeneration(std::vector<std::pair<int,int> >& cosali, int radius, std::vector<std::pair<int,int> >& bound);

double BoundDynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment);

/** Left path is blocked*/
double BoundDynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment);

double CWDynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment, 
							double scale0 = 1.414, double l = 3, double radius = 50, int boundmode = 1, bool left_constrained = false);

double FastInquiry(const std::vector<double>& reference, const std::vector<double>& peer, std::vector<std::pair<int,int> >& alignment, double scale_diff = 8, 
						int seed_num = 3, int seed_length = 128, double radius = 50, int boundmode = 1, bool left_constrained = false);

double DirectInquiry(const std::vector<double>& reference, const std::vector<double>& peer, std::vector<std::pair<int,int> >& alignment, double scale_diff = 8, 
						double radius = 50, int boundmode = 1, bool left_constrained = false);

// double NaiveInquiry(const std::vector<double>& reference, const std::vector<double>& peer, std::vector<std::pair<int, double> >& indexes, int k_th = 3, double ratio = 8,   
// 				   double scale0 = 1.414, double l = 3, double radius = 50, int boundmode = 1, bool left_constrained = false);

/*alignment pair order: <genome_idx, raw_signal_idx */
double EditDistanceError(const std::vector<std::pair<int,int> >& alignment, const std::vector<std::pair<int,int> >& ground_truth, double range = 0.85);

void DumpTimeSeries(const std::vector<double>& signal, const char* path);
}

}


#endif

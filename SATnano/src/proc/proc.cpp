#include "proc.h"
extern "C"{
#include <sys/stat.h>
#undef min
#undef max
}
#include <algorithm>
#include <numeric>
#include <memory.h>
#include <climits>
#include <cfloat>
#include <cmath>
#include <util/exception.h>

#ifndef CWDTW_SOURCE
#define CWDTW_SOURCE
#include "base_fun.cpp"
#endif

void g::proc::ZScoreNormalize(std::vector< double >& signals, double* avg, double* stdev)
{
	double sum = std::accumulate(signals.begin(), signals.end(), 0.0);
	double mean =  sum / signals.size();

	double acc = 0.0;
	for(size_t i = signals.size(); i--;){
		signals[i] = signals[i]-mean;
		acc += signals[i]*signals[i];
	}

	double deviation = std::sqrt(acc/signals.size());
	
	for(size_t i = signals.size(); i--;){
		signals[i] /= deviation;
	}
	
	if(avg){*avg = mean;}
	if(stdev){*stdev = deviation;}
}

#include "resample.cpp"

void g::proc::Resample(int upfactor, int downfactor, const std::vector< double >& input, std::vector< double >& output)
{
	const int n = 10;
	const double bta = 5.0;
	if (upfactor <= 0 || downfactor <= 0){
		throw std::runtime_error("factors must be positive integer");
	}
	int gcd = GetGCD(upfactor, downfactor);
	upfactor /= gcd;
	downfactor /= gcd;

	if (upfactor == downfactor) {
		output = input;
		return;
	}

	int inputSize = input.size();
	int outputSize = QuotientCeil(inputSize * upfactor, downfactor);

	int maxFactor = std::max(upfactor, downfactor);
	double firlsFreq = 1.0 / 2.0 / static_cast<double>(maxFactor);
	int length = 2 * n * maxFactor + 1;
	double firlsFreqs[] = { 0.0, 2.0 * firlsFreq, 2.0 * firlsFreq, 1.0 };
	std::vector<double> firlsFreqsV;
	firlsFreqsV.assign(firlsFreqs, firlsFreqs + 4);
	double firlsAmplitude[] = { 1.0, 1.0, 0.0, 0.0 };
	std::vector<double> firlsAmplitudeV;
	firlsAmplitudeV.assign(firlsAmplitude, firlsAmplitude + 4);
	std::vector<double> coefficients;
	Firls(length - 1, firlsFreqsV, firlsAmplitudeV, coefficients);
	std::vector<double> window;
	Kaiser(length, bta, window);
	int coefficientsSize = coefficients.size();
	for (int i = 0; i < coefficientsSize; i++)
		coefficients[i] *= upfactor * window[i];

	int lengthHalf = (length - 1) / 2;
	int nz = downfactor - lengthHalf % downfactor;
	std::vector<double> h;
	h.reserve(coefficientsSize + nz);
	for (int i = 0; i < nz; i++)
		h.push_back(0.0);
	for (int i = 0; i < coefficientsSize; i++)
		h.push_back(coefficients[i]);
	int hSize = h.size();
	lengthHalf += nz;
	int delay = lengthHalf / downfactor;
	nz = 0;
	while (QuotientCeil((inputSize - 1) * upfactor + hSize + nz, downfactor)
			- delay < outputSize)
		nz++;
	for (int i = 0; i < nz; i++)
		h.push_back(0.0);
	std::vector<double> y;
	upfirdn(upfactor, downfactor, input, h, y);
	output.clear();
	output.reserve(outputSize);
	for (int i = delay; i < outputSize + delay; i++) {
		output.push_back(y[i]);
	}
}

void g::proc::ScaleSignalAmplitude(std::vector< double >& signals, double scale)
{
	for(size_t i = signals.size(); i--;){
		signals[i] *= scale;
	}
}

double maxval(double *x, double *y, int n) {
  int i;
  double tmp, res = std::fabs(x[0] - y[0]);

  for (i=1; i<n; i++) {
    tmp = std::fabs(x[i] - y[i]);
    if (res < tmp) res = tmp;
  }

  return tmp;
}

#define PI 3.1415926535897932384626433

// void g::proc::SpectrumAmplitude(std::vector<double>& signals, std::vector<double>& amplitude)
// {
// 	SpectrumAmplitude(&(signals[0]), signals.size(), amplitude);
// }
// 
// void g::proc::SpectrumAmplitude(double* signals, int size, std::vector< double >& amplitude)
// {
// 	int n = size;
// 	amplitude.resize(n*.5+1);
// 	
// 	double* data = (double*)malloc(sizeof(double)*n);
// 	double* spectrum = (double*)malloc(sizeof(double)*n);
// 
// 	fftw_plan plan = fftw_plan_r2r_1d(n, data, spectrum, FFTW_R2HC, FFTW_PATIENT);
// 	
// 	memcpy(data, &(signals[0]), sizeof(double)*n);
// 	
// 	fftw_execute_r2r(plan, data, spectrum); // signal to spectrum
// 	fftw_destroy_plan(plan); 
// 	
// 	amplitude[0] = spectrum[0];
// 	for(int i = 1; i <= n*.5; i++){
// 		amplitude[i] = sqrt(spectrum[i]*spectrum[i]+spectrum[n-i]*spectrum[n-i]);
// 	}
// 	
// 	fftw_free(data);
// 	fftw_free(spectrum);
// }

void g::proc::RemoveHotSpot(std::vector< double >& signals, int thre)
{
	double mean =  std::accumulate(signals.begin(), signals.end(), 0.0) / signals.size();

	double acc = 0.0;
	for(size_t i = 0; i < signals.size(); i++){
		acc += signals[i]*signals[i];
	}

	double deviation = std::sqrt(acc/signals.size());

	double threshold = mean+thre*deviation;
	
	for(size_t i = 0; i < signals.size(); i++){
		if(signals[i] > threshold){
			signals[i] = threshold;
		}
		else if(signals[i] < -threshold){
			signals[i] = -threshold;
		}
	}
}

//used in WaveDenoise
static double GetThre(double* coeff, int size)  
{  
    double thr = 0.0;  
    double sigma = 0.0;  
	double copy[size];
	
    for(int i = 0; i < size; i++){  
        copy[i] = std::fabs(coeff[i]);
	}
  
    std::sort(copy, copy+size);  
  
    if (size % 2 == 0 && size >= 2){  
        sigma = (copy[size/2-1]+copy[size/2])/2/0.6745; 
	}
    else{
        sigma = copy[size/2]/0.6745;  
	}
  
    double N = size;  
    thr = sigma *sqrt(2.0*log(N));    

    return thr;  
}

// cutoff
static void WThresh(double* coeff, double thre, int size, bool soft)  
{   
	
    if (!soft){  //hard threshold
        for(int i = 0; i < size; i++){  
            if(std::fabs(coeff[i]) < thre){
                coeff[i] = 0.0;  
			}
        }  
    }  
    else{   //soft threshold
        for(int i = 0; i < size; i++){  
            if(std::fabs(coeff[i]) < thre){  
                coeff[i] = 0.0;  
            }  
            else{  
                if(coeff[i] < 0.0){  
                    coeff[i] = thre - std::fabs(coeff[i]);  
				}
                else{  
                    coeff[i] = std::fabs(coeff[i]) - thre;    
				}
            }  
        }  
	}  
} 

static double absmax(double *array, int N){
	double max;
	int i;

	max = 0.0;
	for (i = 0; i < N; ++i) {
		if (std::fabs(array[i]) >= max) {
			max = std::fabs(array[i]);
		}
	}

	return max;
}

static void PyrWThre(wt_object wt, int lvl_thre/*1,2..N*/){
	int J = wt->J;
	int t = wt->length[0];
	
	for(int i = 0; i < J; ++i){
		if(i+1 >= lvl_thre){
			int idx = t;
			int lvlen = wt->length[i+1];
			double thre = GetThre(wt->output+idx, lvlen);
			WThresh(wt->output+idx, thre, lvlen, true);  
		}
// 		printf("Level %d Access : output[%d] Length : %d \n", i + 1,t,wt->length[i+1]);
		t += wt->length[i+1];
	}
}

void g::proc::WaveDenoise(std::vector< double >& signals, bool soft)
{
	double* sigs = &signals[0];		//sst_nino3.dat

	size_t N = signals.size();
	int npyr = 10; 			// Total Number of scales

	wave_object obj = wave_init("db4");
	
	wt_object wt = wt_init(obj, "dwt", N, npyr);// Initialize the wavelet transform object
	setDWTExtension(wt, "sym");// Options are "per" and "sym". Symmetric is the default option
	setWTConv(wt, "direct");
	
	dwt(wt, sigs);// Perform DWT

	PyrWThre(wt, npyr-1);
	
	std::vector<double> output(wt->outlength);
	double* outs = &output[0];
	
	idwt(wt, outs); // Inverse Discrete Wavelet Packet Transform

// 	std::vector<double> diff;
// 	for(int i = 0; i < N; ++i){
// 		diff[i] = (sigs[i] - outs[i]);// /sigs[i];
// 	}
// 	wt_summary(wt); // Tree Summary
// 	printf("\n MAX %g \n", absmax(&diff[0], wt->siglength)); // If Reconstruction succeeded then the output should be a small value.

	wave_free(obj);
	wt_free(wt);
	
	signals = output;
}

void g::proc::MedianFilter(std::vector<double>& signals, int width)
{
	double tmp[width];
	int width_2 = width*.5;
	std::vector<double> nsignl = signals;
	
	for(size_t i = width_2; i < signals.size()-width_2-1; i++){
		memcpy(tmp, &(signals[i-width_2]), sizeof(double)*width);
		std::sort(tmp, tmp+width);
		nsignl[i] = tmp[width_2];
	}
	
	signals = nsignl;
}

void g::proc::Diff(const std::vector<double>& raw, std::vector<double>& diff)
{
	diff.resize(raw.size());
	
	diff[0] = 0;
	
	for(int i = 1; i < raw.size(); i++){
		diff[i] = raw[i]-raw[i-1];
	}
}

double g::proc::L2Distance(const std::vector< double >& sig1, const std::vector< double >& sig2)
{
	double acc = 0;
	
	if(sig1.size() != sig2.size()){
		return -1;
	}
	
	for(int i = 0; i < sig1.size(); i++){
		double dif = sig1[i]-sig2[i];
		acc += dif*dif;
	}
	
	return std::sqrt(acc);
}

double g::proc::L1Distance(const std::vector< double >& sig1, const std::vector< double >& sig2)
{
	double acc = 0;
	
	if(sig1.size() != sig2.size()){
		return -1;
	}
	
	for(int i = 0; i < sig1.size(); i++){
		acc += std::fabs(sig1[i]-sig2[i]);
	}
	
	return acc;
}

void g::proc::LaplaceDiff(const std::vector<double>& raw, std::vector<double>& ldiff)
{
	ldiff.resize(raw.size());
	double tmplt[] = {-1,2,-1};
	int tsize = 3;
	int h_tsize = tsize/2;
	int tbegin = h_tsize;
	int tend = raw.size()-h_tsize;
	
	memset(&ldiff[0], 0, sizeof(double)*ldiff.size());
	
	for(int i = tbegin; i < tend; i++){
		for(int j = 0; j < tsize; j++){
			ldiff[i] += tmplt[j]*raw[j+i-h_tsize];
		}
	}
}

void g::proc::PeakPick(const std::vector<double>& raw, std::vector<std::pair<int, double> >& peaks)
{
	peaks.clear();
	std::vector<double> diff;
	Diff(raw, diff);

	peaks.push_back(std::make_pair(0, raw[0]));
	
	for(int i = 0; i < diff.size()-1; i++){
		if((diff[i] > 0 && diff[i+1] < 0) || (diff[i] < 0 && diff[i+1] > 0)){		//local max min
			peaks.push_back(std::make_pair(i, raw[i]));
		}
	}
	
	peaks.push_back(std::make_pair(raw.size()-1, raw[raw.size()-1]));
}

void g::proc::LocalMaximum(const std::vector<double>& raw, std::vector<std::pair<int, double> >& peaks)
{
	peaks.clear();
	std::vector<double> diff;
	Diff(raw, diff);

	peaks.push_back(std::make_pair(0, raw[0]));
	
	for(int i = 0; i < diff.size()-1; i++){
		if(diff[i] > 0 && diff[i+1] < 0){		//local max min
			peaks.push_back(std::make_pair(i, raw[i]));
		}
	}
	
	peaks.push_back(std::make_pair(raw.size()-1, raw[raw.size()-1]));
}

/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void g::proc::CWTAnalysis(const std::vector<double>& raw, std::vector<std::vector<double> >& output, double scale0, double dscale, int npyr)
{
	const double* sigs = &raw[0];		//sst_nino3.dat
	cwt_object wt;

	size_t N = raw.size();
	double dt = 1;//2;		//sample rate	>  maybe we should use 2?
// 	npyr =  1; 			// Total Number of scales

	wt = cwt_init("dog", 2.0, N, dt, npyr);	//"morlet", "dog", "paul"
	setCWTScales(wt, scale0, dscale, "pow", 2.0);
	
// 	cwt_summary(wt);
	cwt(wt, sigs);

	output.resize(npyr);
	for(size_t k = npyr; k--;){
		int idx = npyr-k-1;
		
		output[idx].resize(raw.size());
		size_t offset = k*raw.size();
		for(size_t i = 0; i < output[idx].size(); i++){
			output[idx][i] = wt->output[i+offset].re; //sqrt(wt->output[i+offset].re*wt->output[i+offset].re+wt->output[i+offset].im*wt->output[i+offset].im);//wt->output[i+offset].re;
		}
	}
	
// 	double *oup;
// 	icwt(wt, oup);

	cwt_free(wt);
}

/** @scale0: level0 pyramind scale;  @dscale: scale_i = scale0*(2^{i*dsacle} ); @npyr: total number of pyramind*/
void g::proc::CWTAnalysisEsp(const std::vector<double>& raw, std::vector<std::vector<double> >& output, double scale0, double dscale, int npyr)
{
	const double* sigs = &raw[0];		//sst_nino3.dat
	cwt_object wt;

	size_t N = raw.size();
	double dt = 1;//2;		//sample rate	>  maybe we should use 2?
// 	npyr =  1; 			// Total Number of scales

	wt = cwt_init("dog", 2.0, N, dt, npyr);	//"morlet", "dog", "paul"
	setCWTScales(wt, scale0, dscale, "pow", 2.0);
	
// 	cwt_summary(wt);
	cwt(wt, sigs);

	output.resize(npyr);
	for(size_t k = npyr; k--;){
		int idx = npyr-k-1;
		
		output[idx].resize(raw.size());
		size_t offset = k*raw.size();
		for(size_t i = 0; i < output[idx].size(); i++){
			output[idx][i] = sqrt(wt->output[i+offset].re*wt->output[i+offset].re+wt->output[i+offset].im*wt->output[i+offset].im);//wt->output[i+offset].re;
		}
	}
	
// 	double *oup;
// 	icwt(wt, oup);

	cwt_free(wt);
}

void g::proc::DumpCWTSpectrogram(const std::vector< std::vector< double > >& waves, char* filename)
{
	std::ofstream o(filename);
	for(int i = 0; i < waves.size(); i++){
		for(int j = 0; j < waves[i].size(); j++){
			o<<waves[i][j]<<"\t";
		}
		o<<std::endl;
	}
	o.close();
}

double g::proc::DynamicTimeWarping(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment)
{	
	double* score[seq1.size()];
	double diff;
	
	for(int i = 0; i < seq1.size(); i++){
		score[i] = new double[seq2.size()];
	}
	
	for(int i = 0; i < seq1.size(); i++){
		for(int j = 0; j < seq2.size(); j++){
			score[i][j] = std::fabs(seq1[i]-seq2[j]);
		}
	}
	
	for(int i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}
	
	for(int j = 1; j < seq2.size(); j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1.size(); i++){
		for(int j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(std::min(score[i-1][j], score[i][j-1]), score[i-1][j-1]);
		}
	}
	
	diff = score[seq1.size()-1][seq2.size()-1];
	
	int i = seq1.size()-1, j = seq2.size()-1;

	alignment.clear();
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
		if(i == 0 && j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	
	std::reverse(alignment.begin(), alignment.end());
	
	for(int i = 0; i < seq1.size(); i++){
		delete [] score[i];
	}
	
	return diff;
}

double g::proc::DynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, std::vector<std::pair<int,int> >& alignment)
{	
	double* score[seq1.size()];
	double diff;
	
	for(int i = 0; i < seq1.size(); i++){
		score[i] = new double[seq2.size()];
	}
	
	for(int i = 0; i < seq1.size(); i++){
		for(int j = 0; j < seq2.size(); j++){
			score[i][j] = std::fabs(seq1[i]-seq2[j]);
		}
	}
	
	for(int i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}
	
	for(int j = 1; j < seq2.size(); j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1.size(); i++){
		for(int j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(score[i][j-1], score[i-1][j-1]);
		}
	}
	
	diff = score[seq1.size()-1][seq2.size()-1];
	
	int i = seq1.size()-1, j = seq2.size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		int ipre = i-1 < 0 ? 0 : i-1;
		int jpre = j-1 < 0 ? 0 : j-1;
		
		double premin = std::min(score[i][jpre], score[ipre][jpre]);
		
		if(premin == score[ipre][jpre]){
			i = ipre; j = jpre;
		}
		if(premin == score[i][jpre]){
			j = jpre;
		}
		if(i == 0 && j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	
	std::reverse(alignment.begin(), alignment.end());
	
	for(int i = 0; i < seq1.size(); i++){
		delete [] score[i];
	}
	
	return diff;
}

double g::proc::SubsequenceDynamicTimeWarping(const std::vector< double >& seq1, const std::vector< double >& seq2, 
											  std::vector<double>& local_score, std::vector<std::pair<int,int> >* alignment)
{
	//-- create score matrix --//
// 	std::cout<<seq1.size()<<" "<<seq2.size()<<std::endl;
	double* score[seq1.size()];
	for(int i = 0; i < seq1.size(); i++){
		score[i] = new double[seq2.size()];
	}

	for(int i = 0; i < seq1.size(); i++){
		for(int j = 0; j < seq2.size(); j++){
			score[i][j] = fabs(seq1[i]-seq2[j]);
		}
	}

	/*-- initialize X-axix --*/
	for(int i = 1; i < seq1.size(); i++){
		score[i][0] += score[i-1][0];
	}

	//-- initialize Y-axix --//
//	for(int j = 1; j < seq2.size(); j++){
//		score[0][j] += score[0][j-1];
//	}

	//-- fill-up score matrix --//
	for(int i = 1; i < seq1.size(); i++){
		for(int j = 1; j < seq2.size(); j++){
			score[i][j] += std::min(std::min(score[i-1][j], score[i][j-1]), score[i-1][j-1]);
		}
	}

	//-- obtain maximal score --//
	double minval = DBL_MAX;
	int min_j = -1;
	
	local_score.resize(seq2.size());
	memcpy(&(local_score[0]), &(score[seq1.size()-1][0]), sizeof(double)*seq2.size());
	
	for(int j = 0; j < seq2.size(); j++){
		if(local_score[j] < minval){
			minval = local_score[j];
			min_j = j;
		}
	}
	
	double diff = local_score[min_j];

	//-- generate alignment --//
	if(alignment){
		alignment->clear();
		int i = seq1.size()-1, j = min_j;
		while(true){
			alignment->push_back(std::make_pair(i,j));
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
				alignment->push_back(std::make_pair(i,j));
				break;
			}
		}
		std::reverse(alignment->begin(), alignment->end());
	}

	for(int i = 0; i < seq1.size(); i++){
		delete [] score[i];
	}

	return diff;
}

static inline double& SCORE(int i, int j, std::vector<double>* score, const std::vector<std::pair<int,int> >& bound){
	static double invalid = DBL_MAX;
	if(bound[i].first <= j && j <= bound[i].second){	
		return score[i][j-bound[i].first];
	}
	else{
		return invalid;
	}
}

void g::proc::DumpTimeSeries(const std::vector< double >& signal, const char* path)
{
	std::ofstream o;
	o.open(path);
	for(int i = 0; i < signal.size(); i++){
		o<<signal[i]<<std::endl;
	}
	o.close();
}

double g::proc::BoundDynamicTimeWarping(const std::vector<double>& sequence1, const std::vector<double>& sequence2, 
										const std::vector<std::pair<int,int> >& bound, std::vector<std::pair<int,int> >& alignment)
{
	/*check order*/
	bool firstorder = true;
	const std::vector<double>* seq1, * seq2;
	
	if(sequence1.size() > sequence2.size()){
		firstorder = false;
	}
	
	if(!firstorder){		//genome first
		seq1 = &sequence2;
		seq2 = &sequence1;
	}
	else{
		seq1 = &sequence1;
		seq2 = &sequence2;
	}
	
	std::vector<double> score[seq1->size()];
	double diff;
	
	for(int i = 0; i < seq1->size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			score[i].push_back(std::fabs((*seq1)[i]-(*seq2)[j]));
		}
	}
	
	for(int i = 1; i < seq1->size(); i++){
		if(SCORE(i, 0, score, bound) == DBL_MAX){
			break;
		}
		
		SCORE(i, 0, score, bound) += SCORE(i-1, 0, score, bound);
	}
	
	for(int j = 1; j <= bound[0].second; j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1->size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			double acc = std::min(std::min(SCORE(i-1, j, score, bound), SCORE(i, j-1, score, bound)), SCORE(i-1, j-1, score, bound));
			if(acc == DBL_MAX){
				score[i][j-bound[i].first] = acc;
			}
			else{
				score[i][j-bound[i].first] += acc;
			}
		}
	}
	
	diff = SCORE(seq1->size()-1, seq2->size()-1, score, bound);
	
	int i = seq1->size()-1, j = seq2->size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		int ipre = i-1 < 0 ? 0 : i-1;
		int jpre = j-1 < 0 ? 0 : j-1;
		
		double scipjp = SCORE(ipre, jpre, score, bound);
		double scipj = SCORE(ipre, j, score, bound);
		double scijp = SCORE(i, jpre, score, bound);
		
		double premin = std::min(std::min(scipj, scijp), scipjp);
		
		if(premin == scipjp){
			i = ipre; j = jpre;
		}
		if(premin == scipj){
			i = ipre;
		}
		if(premin == scijp){
			j = jpre;
		}
		if(i == 0 && j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	
	std::reverse(alignment.begin(), alignment.end());
	
	if(!firstorder){
		for(int i = alignment.size(); i--;){
			std::swap(alignment[i].first, alignment[i].second);
		}
	}
	
	return diff;
}

double g::proc::BoundDynamicTimeWarpingR(const std::vector<double>& seq1, const std::vector<double>& seq2, const std::vector<std::pair<int,int> >& bound, 
										 std::vector<std::pair<int,int> >& alignment)
{	
	std::vector<double> score[seq1.size()];
	double diff;
	
	for(int i = 0; i < seq1.size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			score[i].push_back(std::fabs(seq1[i]-seq2[j]));
		}
	}
	
	for(int i = 1; i < seq1.size(); i++){
		if(SCORE(i, 0, score, bound) == DBL_MAX){
			break;
		}
		
		SCORE(i, 0, score, bound) += SCORE(i-1, 0, score, bound);
	}
	
	for(int j = 1; j <= bound[0].second; j++){
		score[0][j] += score[0][j-1];
	}
	
	for(int i = 1; i < seq1.size(); i++){
		for(int j = bound[i].first; j <= bound[i].second; j++){
			double acc = std::min(SCORE(i, j-1, score, bound), SCORE(i-1, j-1, score, bound));//std::min(std::min(SCORE(i-1, j, score, bound), SCORE(i, j-1, score, bound)), SCORE(i-1, j-1, score, bound));
			if(acc == DBL_MAX){
				score[i][j-bound[i].first] = acc;
			}
			else{
				score[i][j-bound[i].first] += acc;
			}
		}
	}
	
	diff = SCORE(seq1.size()-1, seq2.size()-1, score, bound);
	
	
	int i = seq1.size()-1, j = seq2.size()-1;

	alignment.clear();
	while(true){
		alignment.push_back(std::make_pair(i,j));
		int ipre = i-1 < 0 ? 0 : i-1;
		int jpre = j-1 < 0 ? 0 : j-1;
		
		double scipjp = SCORE(ipre, jpre, score, bound);
// 		double scipj = SCORE(ipre, j, score, bound);
		double scijp = SCORE(i, jpre, score, bound);
		
		double premin = std::min(scijp, scipjp);//std::min(std::min(scipj, scijp), scipjp);
		
		if(premin == scipjp){
			i = ipre; j = jpre;
		}
// 		if(premin == scipj){
// 			i = ipre;
// 		}
		if(premin == scijp){
			j = jpre;
		}
		if(i == 0 || j == 0){
			alignment.push_back(std::make_pair(i,j));
			break;
		}
	}
	
	std::reverse(alignment.begin(), alignment.end());
	
	return diff;
}


static int EditDistance(std::vector<int> str1, std::vector<int> str2) 
{
	int m = str1.size();
	int n = str2.size();
	int dp[m+1][n+1]; 
  
    for(int i = 0; i <= m; i++){ 
		for(int j = 0; j <= n; j++){
            if(i == 0){
                dp[i][j] = j;
			}
            else if(j == 0){ 
                dp[i][j] = i;
			}
 
            else if(str1[i-1] == str2[j-1]){
                dp[i][j] = dp[i-1][j-1]; 
			}
			
            else{
                dp[i][j] = 1 + std::min(std::min(dp[i][j-1], dp[i-1][j]), dp[i-1][j-1]); 
			}
        } 
    } 
  
    return dp[m][n]; 
} 

double g::proc::EditDistanceError(const std::vector< std::pair< int, int > >& alignment, const std::vector< std::pair< int, int > >& ground_truth, double range)
{
	if(alignment.size() < 4){
		return DBL_MAX;
	}
	
	std::vector<std::vector<int> > idx_bag(alignment.back().first+1), ref_bag(ground_truth.back().first+1);
	
	for(int c = 0; c < alignment.size(); c++){
		idx_bag[alignment[c].first].push_back(alignment[c].second);
	}
	
	for(int c = 0; c < ground_truth.size(); c++){
		ref_bag[ground_truth[c].first].push_back(ground_truth[c].second);
	}
	
	double error = 0;
	int begin = idx_bag.size()*(1-range)*.5;
	int end = idx_bag.size()-begin;
	
	for(int i = begin; i < end; i++){
		error += EditDistance(idx_bag[i], ref_bag[i])/(double)ref_bag[i].size();
	}
	
	return error/(idx_bag.size()*range);
}


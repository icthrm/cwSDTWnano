#include <boost/math/special_functions/bessel.hpp>

static int GetGCD(int num1, int num2){
	int tmp = 0;
	while(num1 > 0){
		tmp = num1;
		num1 = num2 % num1;
		num2 = tmp;
	}
	
	return num2;
}

static int QuotientCeil(int num1, int num2){
	if(num1 % num2 != 0){
		return num1 / num2 + 1;
	}
	
	return num1 / num2;
}

static double sinc(double x){
	if(fabs(x - 0.0) < 0.000001){
		return 1;
	}
	
	return sin(M_PI * x) / (M_PI * x);
}

static void Firls(int length, std::vector<double> freq, const std::vector<double>& amplitude, std::vector<double>& result){
	std::vector<double> weight;
	int freq_size = freq.size();
	int weightSize = freq_size / 2;

	weight.reserve(weightSize);
	for(int i = 0; i < weightSize; i++){
		weight.push_back(1.0);
	}

	int filterLength = length + 1;

	for(int i = 0; i < freq_size; i++){
		freq[i] /= 2.0;
	}

	std::vector<double> dFreq;
	for(int i = 1; i < freq_size; i++){
		dFreq.push_back(freq[i] - freq[i - 1]);
	}

	length = (filterLength - 1) / 2;
	int Nodd = filterLength % 2;
	double b0 = 0.0;
	std::vector<double> k;
	if(Nodd == 0){
		for(int i = 0; i <= length; i++)
			k.push_back(i + 0.5);
	}
	else{
		for(int i = 0; i <= length; i++)
			k.push_back(i);
	}

	std::vector<double> b;
	int kSize = k.size();
	for(int i = 0; i < kSize; i++){
		b.push_back(0.0);
	}
	for(int i = 0; i < freq_size; i += 2){
		double slope = (amplitude[i + 1] - amplitude[i])
				/ (freq[i + 1] - freq[i]);
		double b1 = amplitude[i] - slope * freq[i];
		if(Nodd == 1){
			b0 += (b1 * (freq[i + 1] - freq[i])) + slope / 2.0
							* (freq[i + 1] * freq[i + 1] - freq[i] * freq[i])
							* fabs(weight[(i + 1) / 2] * weight[(i + 1) / 2]);
		}
		for(int j = 0; j < kSize; j++){
			b[j] += (slope / (4 * M_PI * M_PI) * (cos(2 * M_PI * k[j] * freq[i + 1])
							- cos(2 * M_PI * k[j] * freq[i])) / (k[j] * k[j]))
					* fabs(weight[(i + 1) / 2] * weight[(i + 1) / 2]);
			b[j] += (freq[i + 1] * (slope * freq[i + 1] + b1)
					* sinc(2 * k[j] * freq[i + 1])
					- freq[i] * (slope * freq[i] + b1)
							* sinc(2 * k[j] * freq[i]))
					* fabs(weight[(i + 1) / 2] * weight[(i + 1) / 2]);
		}
	}
	
	if(Nodd == 1){
		b[0] = b0;
	}
	
	std::vector<double> a;
	double w0 = weight[0];
	for(int i = 0; i < kSize; i++){
		a.push_back((w0 * w0) * 4 * b[i]);
	}
	
	if(Nodd == 1){
		a[0] /= 2;
		for(int i = length; i >= 1; i--){
			result.push_back(a[i] / 2.0);
		}
		result.push_back(a[0]);
		for(int i = 1; i <= length; i++){
			result.push_back(a[i] / 2.0);
		}
	}
	else{
		for(int i = length; i >= 0; i--){
			result.push_back(a[i] / 2.0);
		}
		for(int i = 0; i <= length; i++){
			result.push_back(a[i] / 2.0);
		}
	}
}

void Kaiser(const int order, const double bta, std::vector<double>& window){
	double bes = fabs(boost::math::cyl_bessel_i(0, bta));
	int odd = order % 2;
	double xind = (order - 1) * (order - 1);
	int n = (order + 1) / 2;
	std::vector<double> xi;
	xi.reserve(n);
	for(int i = 0; i < n; i++){
		double val = static_cast<double>(i) + 0.5 * (1 - static_cast<double>(odd));
		xi.push_back(4 * val * val);
	}
	std::vector<double> w;
	w.reserve(n);
	for(int i = 0; i < n; i++){
		w.push_back(boost::math::cyl_bessel_i(0, bta * sqrt(1 - xi[i] / xind)) / bes);
	}
	
	for(int i = n - 1; i >= odd; i--){
		window.push_back(fabs(w[i]));
	}
	for(int i = 0; i < n; i++){
		window.push_back(fabs(w[i]));
	}
}
	
#include <stdexcept>
#include <complex>
#include <vector>
#include <iostream>
#include <cmath>

template<class S1, class S2, class C>
class Resampler{
public:
    typedef    S1 inputType;
    typedef    S2 outputType;
    typedef    C coefType;

    Resampler(int upRate, int downRate, C *coefs, int coefCount);
    virtual ~Resampler();

    int        apply(S1* in, int inCount, S2* out, int outCount);
    int        neededOutCount(int inCount);
    int        coefsPerPhase(){ return _coefsPerPhase; }
    
private:
    int        _upRate;
    int        _downRate;

    coefType   *_transposedCoefs;
    inputType  *_state;
    inputType  *_stateEnd;
    
    int        _paddedCoefCount;  // ceil(len(coefs)/upRate)*upRate
    int        _coefsPerPhase;    // _paddedCoefCount / upRate
    
    int        _t;                // "time" (modulo upRate)
    int        _xOffset;
    
};

/*The coefficients are copied into local storage in a transposed, flipped
  arrangement.  For example, suppose upRate is 3, and the input number
  of coefficients coefCount = 10, represented as h[0], ..., h[9].
  Then the internal buffer will look like this:
        h[9], h[6], h[3], h[0],   // flipped phase 0 coefs
           0, h[7], h[4], h[1],   // flipped phase 1 coefs (zero-padded)
           0, h[8], h[5], h[2],   // flipped phase 2 coefs (zero-padded)   */
template<class S1, class S2, class C>
Resampler<S1, S2, C>::Resampler(int upRate, int downRate, C *coefs,
                                int coefCount): _upRate(upRate), _downRate(downRate), _t(0), _xOffset(0)
{
    _paddedCoefCount = coefCount;
    while (_paddedCoefCount % _upRate){
        _paddedCoefCount++;
    }
    _coefsPerPhase = _paddedCoefCount / _upRate;
    
    _transposedCoefs = new coefType[_paddedCoefCount];
    std::fill(_transposedCoefs, _transposedCoefs + _paddedCoefCount, 0.);

    _state = new inputType[_coefsPerPhase - 1];
    _stateEnd = _state + _coefsPerPhase - 1;
    std::fill(_state, _stateEnd, 0.);


    /* This both transposes, and "flips" each phase, while
     * copying the defined coefficients into local storage.
     * There is probably a faster way to do this
     */
    for (int i=0; i<_upRate; ++i){
        for (int j=0; j<_coefsPerPhase; ++j){
            if(j*_upRate + i  < coefCount)
                _transposedCoefs[(_coefsPerPhase-1-j) + i*_coefsPerPhase] =
                                                coefs[j*_upRate + i];
        }
    }
}

template<class S1, class S2, class C>
Resampler<S1, S2, C>::~Resampler(){
    delete [] _transposedCoefs;
    delete [] _state;
}

/* compute how many outputs will be generated for inCount inputs  */
template<class S1, class S2, class C>
int Resampler<S1, S2, C>::neededOutCount(int inCount)
{
    int np = inCount * _upRate;
    int need = np / _downRate;
    if((_t + _upRate * _xOffset) < (np % _downRate))
        need++;
    return need;
}

template<class S1, class S2, class C>
int Resampler<S1, S2, C>::apply(S1* in, int inCount, S2* out, int outCount)
{
    if(outCount < neededOutCount(inCount)) 
        throw std::invalid_argument("Not enough output samples");

    // x points to the latest processed input sample
    inputType *x = in + _xOffset;
    outputType *y = out;
    inputType *end = in + inCount;
    while (x < end){
        outputType acc = 0.;
        coefType *h = _transposedCoefs + _t*_coefsPerPhase;
        inputType *xPtr = x - _coefsPerPhase + 1;
        int offset = in - xPtr;
        if(offset > 0){
            // need to draw from the _state buffer
            inputType *statePtr = _stateEnd - offset;
            while (statePtr < _stateEnd){
                acc += *statePtr++ * *h++;
            }
            xPtr += offset;
        }
        while (xPtr <= x){
            acc += *xPtr++ * *h++;
        }
        *y++ = acc;
        _t += _downRate;

        int advanceAmount = _t / _upRate;

        x += advanceAmount;
        // which phase of the filter to use
        _t %= _upRate;
    }
    _xOffset = x - end;

    // manage _state buffer
    // find number of samples retained in buffer:
    int retain = (_coefsPerPhase - 1) - inCount;
    if(retain > 0){
        // for inCount smaller than state buffer, copy end of buffer
        // to beginning:
        std::copy(_stateEnd - retain, _stateEnd, _state);
        // Then, copy the entire (short) input to end of buffer
        std::copy(in, end, _stateEnd - inCount);
    }
    else{
        // just copy last input samples into state buffer
        std::copy(end - (_coefsPerPhase - 1), end, _state);
    }
    // number of samples computed
    return y - out;
}

/** This template function provides a one-shot resampling. Extra samples are padded to the end of the input in order to capture all of the non-zero 
output samples. The output is in the "results" vector which is modified by the function.

Note, I considered returning a vector instead of taking one on input, but then the C++ compiler has trouble with implicit template instantiation
(e.g. have to say upfirdn<float, float, float> every time - this way we can let the compiler infer the template types).

Thanks to Lewis Anderson (lkanders@ucsd.edu) at UCSD for the original version of this function. */
template<class S1, class S2, class C>
void upfirdn(int upRate, int downRate, const S1 *input, int inLength, C *filter, int filterLength, std::vector<S2> &results)
{
    // Create the Resampler
    Resampler<S1, S2, C> theResampler(upRate, downRate, filter, filterLength);

    // pad input by length of one polyphase of filter to flush all values out
    int padding = theResampler.coefsPerPhase() - 1;
    S1 *inputPadded = new S1[inLength + padding];
    for (int i = 0; i < inLength + padding; i++){
        if(i < inLength)
            inputPadded[i] = input[i];
        else
            inputPadded[i] = 0;
    }

    // calc size of output
    int resultsCount = theResampler.neededOutCount(inLength + padding); 

    results.resize(resultsCount);

    // run filtering
    int numSamplesComputed = theResampler.apply(inputPadded, inLength + padding, &results[0], resultsCount);
    delete[] inputPadded;
}

/** This template function provides a one-shot resampling. The output is in the "results" vector which is modified by the function.
In this version, the input and filter are vectors as opposed to pointer/count pairs. */
template<class S1, class S2, class C>
void upfirdn(int upRate, int downRate, const std::vector<S1> &input, std::vector<C> &filter, std::vector<S2> &results)
{
    upfirdn<S1, S2, C>(upRate, downRate, &input[0], input.size(), &filter[0], filter.size(), results);
}

/*
 * conv.h
 *
 *  Created on: May 1, 2013
 *      Author: Rafat Hussain
 */

#ifndef CONV_H_
#define CONV_H_

#include "real.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

typedef struct conv_set* conv_object;

conv_object conv_init(int N, int L);

struct conv_set{
	fft_real_object fobj;
	fft_real_object iobj;
	int ilen1;
	int ilen2;
	int clen;
};

int factorf(int M);

int findnext(int M);

int findnexte(int M);

void conv_direct(fft_type *inp1,int N, fft_type *inp2, int L,fft_type *oup);

void conv_directx(fft_type *inp1,int N, fft_type *inp2, int L,fft_type *oup);

//void conv_fft(const conv_object obj,fft_type *inp1,fft_type *inp2,fft_type *oup);

//void conv_fft(const conv_object obj,fft_type *inp1,fft_type *inp2,fft_type *oup);

void conv_fft(const conv_object obj,fft_type *inp1,fft_type *inp2,fft_type *oup);

//void free_conv(conv_object object);

void free_conv(conv_object object);


#ifdef __cplusplus
}
#endif

#endif /* CONV_H_ */

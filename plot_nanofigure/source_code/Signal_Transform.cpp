#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
using namespace std;

//-------- utility ------//
void getBaseName(string &in,string &out,char slash,char dot)
{
	int i,j;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	i++;
	for(j=len-1;j>=0;j--)
	{
		if(in[j]==dot)break;
	}
	if(j==-1)j=len;
	out=in.substr(i,j-i);
}
void getRootName(string &in,string &out,char slash)
{
	int i;
	int len=(int)in.length();
	for(i=len-1;i>=0;i--)
	{
		if(in[i]==slash)break;
	}
	if(i<=0)out=".";
	else out=in.substr(0,i);
}


//------------ load signal ------------//
//-> example
/*
478
571
531
542
536
541
534
549
*/
int Load_Signal(string &input, vector <double> &signal)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",input.c_str());
		exit(-1);
	}
	//load
	signal.clear();
	double raw_signal;
	int count=0;
	for(;;)
	{
		if(! (fin>>raw_signal) )break;
		signal.push_back(raw_signal);
		//count++
		count++;
	}
	//return
	return count;
}

//-------- load start,end,char ----------//
//-> example (1-base)
/*
1 2 C
2 6 C
6 55 C
55 64 T
64 68 G
68 69 G
*/
int Load_Event(string &input, vector <pair<int,int> > &range, vector <char> &code)
{
	ifstream fin;
	string buf,temp;
	//read
	fin.open(input.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"file %s not found!\n",input.c_str());
		exit(-1);
	}
	//load
	range.clear();
	code.clear();
	int count=0;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		int start,end;
		char c;
		www>>start>>end>>c;
		start--;
		end--;
		end--;
		range.push_back(pair<int,int>(start,end));
		code.push_back(c);
		//count++
		count++;
	}
	//return
	return count;
}


//------- signal transform (written by Renmin Han) --------------//
void Signal_Transform(vector<double>& input_signal, vector<double>& output_signal, int output_len = 50)
{
	//------ get input signal length -----//
	int input_len = input_signal.size();
	
	if(input_len == 1){
		output_signal.resize(output_len);
		for(int i = 0; i < output_len; i++){
			output_signal[i] = input_signal[0];
		}
		return;
	}
	
	int interval = output_len/(input_len-1);
	
	int intev[input_len];
	
	for(int i = 1; i < input_len; i++){
		intev[i] = interval;
	}
	
	int remains = output_len-interval*(input_len-1);
	for(int i = 1; remains > 0; i++){
		intev[i]++;
		remains--;
	}
	
	//------ output signal -------//
	output_signal.resize(output_len);
	
	int acc = 0;
	
	for(int i = 1; i < input_len; i++){
		double pre = input_signal[i-1];
		double cur = input_signal[i];
		
		double diff = (cur-pre)/intev[i];
		for(int k = 0; k < intev[i]; k++){
			output_signal[acc] = pre+diff*k;
			acc++;
			if(acc >= output_len){
				break;
			}
		}
		if(acc >= output_len){
			break;
		}
	}
	
	while(acc < output_len){
		output_signal[acc] = output_signal[acc-1];
		acc++;
	}
}

//================ main process ===================//
void Main_Process(string &signal, string &label, string &output, int len)
{
	//---- load signal -----//
	vector <double> sig;
	int sig_len=Load_Signal(signal, sig);
	//---- load label ------//
	vector <pair<int,int> > lab;
	vector <char> code;
	int lab_len=Load_Event(label, lab, code);
	//---- process -----//
	vector <double> out_sig;
	vector <pair<int,int> > out_lab;
	for(int i=0;i<lab_len;i++)
	{
		vector <double> insig;
		vector <double> outsig;
		//-> get current signal
		int start=lab[i].first;
		int termi=lab[i].second;
		insig.resize(termi-start+1);
		int cur=0; 
		for(int k=start;k<=termi;k++)
		{
			insig[cur]=sig[k];
			cur++;
		}
		//-> transform signal
		Signal_Transform(insig, outsig, len);
		//-> insert to output signal
		out_sig.insert(out_sig.end(),outsig.begin(),outsig.end());
		//-> insert to output label
		out_lab.push_back(pair<int,int>(i*len,i*len+len));
	}
 	//---- output -----//
 	FILE *fp;
 	string file;
 	//-> output transformed signal
 	file=output+".signal";
 	fp=fopen(file.c_str(),"wb");
 	for(int i=0;i<(int)out_sig.size()-1;i++)fprintf(fp,"%lf ",out_sig[i]);
	fprintf(fp,"%lf\n",out_sig[(int)out_sig.size()-1]);
 	fclose(fp);
 	//-> output transformed label
 	file=output+".label";
 	fp=fopen(file.c_str(),"wb");
 	for(int i=0;i<lab_len;i++)fprintf(fp,"%d %d %c \n",out_lab[i].first,out_lab[i].second,code[i]);
 	fclose(fp);
}


//---------- main ----------//
int main(int argc,char **argv)
{
	//---- Signal_Transform ----//
	{
		if(argc<5)
		{
			fprintf(stderr,"Signal_Transform <input_signal> <input_label> <transform_len> <out_name> \n");
			fprintf(stderr,"[note]: input_label should be 'start termi code' in 1-base \n");
			fprintf(stderr,"        transform_len should be set to 50 \n");
			exit(-1);
		}
		string input_signal=argv[1];
		string input_label=argv[2];
		int transform_len=atoi(argv[3]);
		string out_name=argv[4];
		//--- process ---//
		Main_Process(input_signal,input_label,out_name,transform_len);
		//--- exit ---//
		exit(0);
	}
}


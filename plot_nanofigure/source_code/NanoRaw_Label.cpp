#include <string>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <map>
#include <algorithm>
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

int Parse_Str(string &in,vector <string> &out)
{
	istringstream www(in);
	out.clear();
	int count=0;
	for(;;)
	{
		string buf;
		if(! (www>>buf) )break;
		out.push_back(buf);
		count++;
	}
	return count;
}

//----- Z-score ----//
void Z_score(vector <double> &in,double &mean,double &vari)
{
	//init
	mean=0;
	vari=1;
	//proc
	int i;
	int size=(int)in.size();
	if(size==0)return;
	//-> calculate mean
	for(i=0;i<size;i++)mean+=in[i];
	mean/=size;
	//-> calculate vari
	vari=0;
	for(i=0;i<size;i++)vari+=(in[i]-mean)*(in[i]-mean);
	vari=1.0*sqrt(vari/size);
	//-> calculate Z-score
	for(i=0;i<size;i++)in[i]=(in[i]-mean)/vari;
}



//----------- generate NanoRaw-style label output ----------//
//-> input
/*
  385   288 |         91          9 |       -0.915353        -1.17514          diff:        0.25979   CAGAA 9
  303   131 |         92         10 |        -2.08046        -1.82747          diff:        0.25299   AGAAT 10
  285   527 |         93         11 |        -2.33622        -1.86196          diff:       0.474257   GAATT 11
  278   527 |         94         11 |        -2.43568        -1.86196          diff:       0.573718   GAATT 11
  294   527 |         95         11 |        -2.20834        -1.86196          diff:        0.34638   GAATT 11
  369   527 |         96         11 |        -1.14269        -1.86196          diff:       0.719268   GAATT 11
  433    63 |         97         12 |       -0.233338        0.427542          diff:        0.66088   AATTT 12
  399    63 |         98         12 |       -0.716432        0.427542          diff:        1.14397   AATTT 12
  481    63 |         99         12 |        0.448676        0.427542          diff:      0.0211343   AATTT 12
....
*/
//-> output
/*
107 111 C
111 123 G
123 131 C
131 160 C
....
*/

void NanoRaw_Label(string &in)
{
	//read
	ifstream fin;
	string buf,temp;
	fin.open(in.c_str(), ios::in);
	if(fin.fail()!=0)
	{
		fprintf(stderr,"list %s not found!!\n",in.c_str());
		exit(-1);
	}
	//read
	int first=1;
	string prev="";
	string curr="";
	string prev_rec="";
	int start=0;
	int count=0;
	vector <double> signal;
	for(;;)
	{
		if(!getline(fin,buf,'\n'))break;
		istringstream www(buf);
		vector <string> out;
		int retv=Parse_Str(buf,out);
		curr=out[retv-1];
		if(first==1)
		{
			prev=curr;
			prev_rec=out[retv-2];
			first=0;
		}
		else
		{
			if(prev!=curr) //-> output
			{
				printf("%d %d %c\n",start+1,count+1,prev_rec[2]);
				start=count;
				prev=curr;
				prev_rec=out[retv-2];
			}
		}
		//fprintf(stderr,"%s\n",out[0].c_str());
		signal.push_back(atof(out[0].c_str()));
		count++;
	}
	//termi
	if(prev_rec!="") //-> output
	{
		printf("%d %d %c\n",start+1,count+1,prev_rec[2]);
	}
	//--- transfer to Z-score ---//
	double mean,vari;
	Z_score(signal,mean,vari);
	for(int i=0;i<(int)signal.size();i++)fprintf(stderr,"%lf\n",signal[i]);
}


//---------- main ----------//
int main(int argc,char **argv)
{
	//---- NanoRaw_Label ----//
	{
		if(argc<2)
		{
			fprintf(stderr,"NanoRaw_Label <nored_signal> \n");
			exit(-1);
		}
		string infile=argv[1];
		//proc
		NanoRaw_Label(infile);
		//exit
		exit(0);
	}
}


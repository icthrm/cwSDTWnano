import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


#----- function int2kmer --------#
def int2kmer(input_num,mer_k = 5):
	kmer = ''
	kmer_dict = ['A','C','G','T']
	for i in range(mer_k):
		kmer = kmer_dict[input_num%4]+kmer
		input_num = input_num/4
	return kmer



#=========== <signal> example ===========#
"""
 417 406 403 407 426 501 497 508 492 493 
"""

#----- load signal data ----#
def load_signal(infile):
	#-> read infile
	fh = open(infile)
	nanoraw_signal = []
	for line in list(fh):
		signal = line.split(' ')
		nanoraw_signal.append(signal)
	#-> return
	return nanoraw_signal


#=========== <start end char> example ===========#
"""
1 2 C
2 6 C
6 55 C
55 64 T
64 68 G
68 69 G
"""

#----- load range data ----#
def load_range(infile, base = 0):
	#-> read infile
	fh = open(infile)
	input_beg = list()
	input_end = list()
	input_lab = list()
	for line in fh:
		split_line = line.split()
		input_beg.append(int(split_line[0])+base)
		input_end.append(int(split_line[1])+base)
		input_lab.append(split_line[2])
	#-> return
	return input_beg,input_end,input_lab



#============= plot the figure =============#
def plot_figure(signal, beg,end,lab, output, start=None, termi=None):
	size=2000
	part=50
	color_dict = {'A':'red','C':'yellow','G':'green','T':'blue'}
	orilen=len(signal[0])
	#-> determine start and end if None
	if start is None:
		start=0
	if termi is None:
		termi=orilen-1
	if start < 0:
		start=0
	if termi >= orilen:
		termi=orilen-1
	#-> get subsequence of the input signal
	length=termi-start+1
	if length%size == 0:
		col=length/size
	else:
		col=int(length/size)+1
	#-> set the x-axis
	scalar=np.array([0, 5, 10, 15, 20, 25, 30, 35, 40])
	srange=scalar*50

	#-> plot the figure
	f,axs = plt.subplots(col,1,sharex = False,figsize=(20, 4*col), squeeze=False)
	for i in range(col):
		#--| get cur_start and cur_end
		cur_start=start+i*size
		cur_end=min(cur_start+size,termi)
		#--| plot cursignal
		cur=0
		for k in signal:
			cursignal=np.zeros(size)
			cursignal[0:cur_end-cur_start]=k[cur_start:cur_end]
			#---- plot 1st and 2nd signal with bold color ----#
			if cur == 0:
				axs[i,0].plot(cursignal, linewidth=5.0, color='red')
			elif cur == 1:
				axs[i,0].plot(cursignal, linewidth=2.5, color='green')
			else:
				axs[i,0].plot(cursignal, linewidth=0.2, color='black')
			cur=cur+1
		axs[i,0].set_ylim([-3,3])
		#--| plot curlabel
		for idx,val in enumerate(beg):
			if val>=cur_start and val<cur_end:
				axs[i,0].axvline(val-cur_start, color='grey', linewidth=0.05)
				axs[i,0].axvline(end[idx]-cur_start, color='grey', linewidth=0.05)
				axs[i,0].axvspan(val-cur_start,end[idx]-cur_start, ymax=0.05, alpha=0.5,color=color_dict[lab[idx]])
		#--| plot x-axis scale
		plt.sca(axs[i, 0])
		plt.xticks(srange, scalar+i*40)

	#-> save the figure
	f.savefig(output)



#------- usage -------#
def Usage():
	print 'python PlotSignal.py <raw_signals> <start_end_char> <output.png> <start> <termi>'
	print 'signal_label    : the file contains signal and label '
	print 'start_end_char  : the file contains start, end, and char '
	print 'output.png      : output PNG file'


#------- main --------#
def main(argv):
	if len(argv) != 5:
		Usage()
		sys.exit(-1)
	#-> get input arguments
	signal_label=argv[0]
	start_end_char=argv[1]
	output=argv[2]
	start=int(argv[3])
	termi=int(argv[4])
	#-> load input files
	nanoraw_signal=load_signal(signal_label)
	input_beg,input_end,input_lab=load_range(start_end_char)
	#-> process
	plot_figure(nanoraw_signal,input_beg,input_end,input_lab,output,start,termi)


#-------------- python main ---------------#
if __name__ == "__main__":
	## read <signal label> and <start end char>, plot the figure
	main(sys.argv[1:])




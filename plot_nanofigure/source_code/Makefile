

# make the subdir first
CXX = g++
#FLAG = -O3 -static
FLAG = -O3
SHELL :=/bin/bash

all:
	# ------ all files -------#
	$(CXX) NanoRaw_Label.cpp $(FLAG) -o NanoRaw_Label
	$(CXX) Signal_Transform.cpp $(FLAG) -o Signal_Transform
	mkdir -p ../util
	mv NanoRaw_Label Signal_Transform ../util

# clean all the temp files
clean:

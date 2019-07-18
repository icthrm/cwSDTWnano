project(SATnano) #Signal-level Analysis Toolkit for nanopore sequencing

cmake_minimum_required(VERSION 2.6)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no
C++11 support. Please use a different C++ compiler.")
endif()

#the store of data in clapack or cblas is column main sequence

set(WAVELET_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/wavelib/header)
set(FAST5_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/fast5/include)
#set(HDF5_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/hdf5_1.8.16/src)
set(CLAPACK_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/INCLUDE)
set(BLASWRAP_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/clapack_3.1/BLAS/WRAP)
set(CBLAS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/cblas/INCLUDE)
set(MATRIX_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/matrix)
# set(STATS_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/lib/stats/include)

include_directories(${PROJECT_SOURCE_DIR}/lib ${PROJECT_SOURCE_DIR}/src ${WAVELET_INCLUDE_DIR} ${HDF5_INCLUDE_DIR} 
	${FAST5_INCLUDE_DIR} ${CLAPACK_INCLUDE_DIR} ${BLASWRAP_INCLUDE_DIR} ${CBLAS_INCLUDE_DIR} ${MATRIX_INCLUDE_DIR})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_subdirectory(lib)
add_subdirectory(src)

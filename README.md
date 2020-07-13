# cwSDTWnano
Here we proposed two novel algorithms, the Direct Subsequence Dynamic Time Warping for nanopore raw signal search (DSDTWnano) and the continuous wavelet Subsequence DTW for nanopore raw signal search (cwSDTWnano), to enable the direct subsequence inquiry and exact mapping in the nanopore raw signal datasets. The proposed algorithms are based on the idea of Subsequence-extended Dynamic Time Warping (SDTW) and directly operates on the raw signals, without any loss of information. DSDTWnano could ensure an output of highly accurate query result and cwSDTWnano is the accelerated version of DSDTWnano, with the help of seeding and multi-scale coarsening of signals that based on continuous wavelet transform (CWT).

The source code of cwSDTWnano is within the project SATnano (sat-query), as well as cwdtw (sat-align). To compile the files, please use cmake. The updated code supports fast5 format, which require the installed hdf5 lib (here, for centos/fedora, please use the commend "dnf install hdf5-devel"; for ubuntu, please use the commend "apt-get install libhdf5-serial-dev"; otherwise, manually download from https://www.hdfgroup.org/downloads/hdf5/). To compile the source, please input the following code in terminal:

cd ./cwSDTWnano
mkdir build
cd build
cmake ..
make

# Example

sat-query -i genome_subseq -p nanopore_sig -o out -m 0

Here, -m 0 for direct query (DSDTWnano) and -m 1 for fast query (cwSDTWnano). The default mode is -m 0. More detailed information could be get by -h.

More testing data can be found at 
https://drive.google.com/drive/folders/1LuOxg9qE1l9AuDcfyUz9aF10X4cgmX5t?usp=sharing


# Output format:

         1          2                3               4             5               6                            7                 8
--------------------------------------------------------------------------------------------------------------------------------------
         0     56243 |         88.2554             514 |        -0.21721        0.227072          diff:       0.444281          AAACAA
         1     56244 |         106.833             602 |         1.26698         1.40997          diff:       0.142982          AACAAG
         1     56245 |         106.833             611 |         1.26698         1.53094          diff:        0.26396          AACAAG
         1     56246 |         106.833             614 |         1.26698         1.57127          diff:       0.304286          AACAAG
         1     56247 |         106.833             606 |         1.26698         1.46373          diff:        0.19675          AACAAG
         1     56248 |         106.833             610 |         1.26698          1.5175          diff:       0.250518          AACAAG
         1     56249 |         106.833             598 |         1.26698          1.3562          diff:      0.0892146          AACAAG
         1     56250 |         106.833             580 |         1.26698         1.11424          diff:       0.152741          AACAAG
         2     56251 |         90.1578             480 |      -0.0652257       -0.229956          diff:        0.16473          ACAAGG
         2     56252 |         90.1578             491 |      -0.0652257       -0.082094          diff:      0.0168683          ACAAGG
         2     56253 |         90.1578             495 |      -0.0652257      -0.0283261          diff:      0.0368997          ACAAGG
         2     56254 |         90.1578             514 |      -0.0652257        0.227072          diff:       0.292297          ACAAGG
         2     56255 |         90.1578             497 |      -0.0652257      -0.0014421          diff:      0.0637836          ACAAGG
         
---------
[Legend]:
---------

the 1st column shows the mapping of the first position (e.g., on query signal) starting from 0.
the 2nd column shows the mapping of the second position (e.g., on raw signal) starting from 0.

the 3rd column displays the original value of the first input signal (e.g., on query signal).
the 4th column displays the original value of the second input signal (e.g., on raw signal).

the 5th column indicates the Z-normalized value of the first input signal (e.g., on query signal).
the 6th column indicates the Z-normalized value of the second input signal (e.g., on raw signal).

the 7th column illustrates the absolute difference between the two Z-normalized values.
the 8th column displays the query sequence 6-mer at the index of 1st column.

# Example2

sat-query -i genome_subseq -p nanopore_sig -o out -m 0 -f

Here, -m 0 for direct query (DSDTWnano) and -m 1 for fast query (cwSDTWnano). The default mode is -m 0. More detailed information could be get by -h. -f is added to achieve a concise output.

# Output format:

0.161953

 query_idx      signal_idx      query_value    signal_value         6-mer              
--------------------------------------------------------------------------
         0          22272 |         75.9828             414         CTGATA
         0          22273 |         75.9828             409         CTGATA
         0          22274 |         75.9828             413         CTGATA
         0          22275 |         75.9828             423         CTGATA
         0          22276 |         75.9828             424         CTGATA
         0          22277 |         75.9828             411         CTGATA
         0          22278 |         75.9828             423         CTGATA
         0          22279 |         75.9828             412         CTGATA
         0          22280 |         75.9828             414         CTGATA
         0          22281 |         75.9828             406         CTGATA
         1          22282 |         77.0378             454         TGATAA
         1          22283 |         77.0378             495         TGATAA
         2          22284 |         109.129             585         GATAAG
         2          22285 |         109.129             556         GATAAG
         2          22286 |         109.129             585         GATAAG
         2          22287 |         109.129             573         GATAAG
         2          22288 |         109.129             578         GATAAG
         2          22289 |         109.129             552         GATAAG
         3          22290 |         90.9777             544         ATAAGT
         3          22291 |         90.9777             473         ATAAGT
         4          22292 |         79.1267             437         TAAGTA
         4          22293 |         79.1267             447         TAAGTA
         4          22294 |         79.1267             438         TAAGTA
         5          22295 |          77.238             427         AAGTAA
         5          22296 |          77.238             431         AAGTAA



         
---------
[Legend]:
---------

query_idx column shows the mapping of the first position (e.g., on query read).
signal_idx column shows the mapping of the second position (e.g., on raw signal) starting from 0.

query_value column displays the original value of the first input signal (e.g., the standard signal value corresponding to the time point on query read).
signal_value column displays the original value of the second input signal (e.g., on raw signal).

The format could be accepted by "plot_nanofigure", a small tool to display the alignment result between read and raw signal that produced by sat-query.

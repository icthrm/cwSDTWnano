# cwSDTWnano
Here we proposed two novel algorithms, the Direct Subsequence Dynamic Time Warping for nanopore raw signal search (DSDTWnano) and the continuous wavelet Subsequence DTW for nanopore raw signal search (cwSDTWnano), to enable the direct subsequence inquiry and exact mapping in the nanopore raw signal datasets. The proposed algorithms are based on the idea of Subsequence-extended Dynamic Time Warping (SDTW) and directly operates on the raw signals, without any loss of information. DSDTWnano could ensure an output of highly accurate query result and cwSDTWnano is the accelerated version of DSDTWnano, with the help of seeding and multi-scale coarsening of signals that based on continuous wavelet transform (CWT).

# example

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

the 1st column shows the mapping of the first position (e.g., on query signal) starting from 1,
the 2nd column shows the mapping of the second position (e.g., on raw signal) starting from 1,

the 3rd column displays the original value of the first input signal (e.g., on query signal),
the 4th column displays the original value of the second input signal (e.g., on raw signal),

the 5th column indicates the Z-normalized value of the first input signal (e.g., on query signal),
the 6th column indicates the Z-normalized value of the second input signal (e.g., on raw signal),

the 7th column illustrates the absolute difference between the two Z-normalized values.
the 8th column displays the query sequence 6-mer at the index of 1st column.



# :dna:ParallelRNA:dna:: Parallel RNA structure counting algorithm 
Secondary structures for nucleic acid sequences, DNA and RNA, are useful abstractions when examining certain behaviors of those molecules. We examine the problem of counting secondary structures compatible with an ordered multiset of sequences. In particular, we address the issue of accounting for indistinguishable secondary structures, for which no fast algorithm has been found. We provide a parallel algorithm for counting distinguishable secondary structures.

## 0. Prerequisites
g++ <br />
openmp <br />

## 0.1. How to compile and run
```
g++ -o max_pair_omp src/max_pair_omp.cpp -fopenmp 

./max_pair_omp input.txt num_threads
```

## 1. Diagonal Algorithm
The first type of parallel algorithm breaks up each diagonal into chunks that each OMP thread solves individually. Once a chunk is done, it waits for the remaining threads to finish to move on to the next diagonal section.
<br>
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/figure.jpg" width="300"> <br>
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/parallelRNA_vis.gif" width="600">


## 2. Triangle Algorithm
The second type of parallel algorithm breaks up each diagonal into chunks and computes a full triangle of the matrix. After it fills this triangle and all threads finish, it computes the inverted triangle such that after a triangle and inverted triangle are complete, the matrix has a perfect diagonal on which it can run the algorithm again.
<br>
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/tri-gif.gif" width="600">

## 3. Unsupervised Algorithm
Finally, we include an unsupervised algorithm. The two previous algorithms tell OMP the size of chunk each thread should work on. However, sometimes it is faster to let OMP assign threads tasks on its own. This algorithm is functionally similar to the Diagonal Algorithm, except it does not define a specific chucnk length or starting and ending position, it simply fills in the diagonals in whichever order OMP sees fit.

## Author Contributions
<!--## 2. Files-->
<!--The following files are included in this folder, in addition to this readme-->
<!--file, readme.md.-->
<!--<ul>-->
<!--<li>md.c: Main C program</li>-->
<!--<li>md.h: Header file for md.c</li>-->
<!--<li>md.in: Input parameter file (to be redirected to the standard input)</li>-->
<!--</ul>-->
<!--![Screen shot of MD simulation](ScreenShot.png)-->

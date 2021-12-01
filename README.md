# :dna:ParallelRNA:dna:: Parallel RNA structure counting algorithm 
Secondary structures for nucleic acid sequences, DNA and RNA, are useful abstractions when examining certain behaviors of those molecules. We examine the problem of counting secondary structures compatible with an ordered multiset of sequences. In particular, we address the issue of accounting for indistinguishable secondary structures, for which no fast algorithm has been found. We provide a parallel algorithm for counting distinguishable secondary structures.
<br>
<!-- ![alt text](https://github.com/masarunakajima/parallelRNA/blob/openMP/parallel%20rna.PNG) -->
<!-- ![alt text](https://github.com/masarunakajima/parallelRNA/blob/openMP/total%20fig.jpg) -->
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/figure.jpg" width="300">
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/parallelRNA_vis.gif" width="600">

## 0. Prerequisites
gcc/8.3.0 <br />
openmp <br />
Boost library

## 1. How to compile and run
```
gcc -o count_OMP -I/path/to/boost src/count_OMP.cpp -fopenmp \
-lstdc++ -L/path/to/boost/stage/lib/ -lboost_mpi -lboost_serialization

./count_OMP input.txt
```

## 2. How to compile and run for measurement
```
gcc -o count_OMP_adaptive -I/path/to/boost src/count_OMP_adaptive.cpp -fopenmp \
-lstdc++ -L/path/to/boost/stage/lib/ -lboost_mpi -lboost_serialization

./count_OMP_adaptive input.txt num_threads
```
<!--## 2. Files-->
<!--The following files are included in this folder, in addition to this readme-->
<!--file, readme.md.-->
<!--<ul>-->
<!--<li>md.c: Main C program</li>-->
<!--<li>md.h: Header file for md.c</li>-->
<!--<li>md.in: Input parameter file (to be redirected to the standard input)</li>-->
<!--</ul>-->
<!--![Screen shot of MD simulation](ScreenShot.png)-->

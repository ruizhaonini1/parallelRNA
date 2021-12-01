# :dna:ParallelRNA:dna:: Parallel RNA structure counting algorithm 
Secondary structures for nucleic acid sequences, DNA and RNA, are useful abstractions when examining certain behaviors of those molecules. We examine the problem of counting secondary structures compatible with an ordered multiset of sequences. In particular, we address the issue of accounting for indistinguishable secondary structures, for which no fast algorithm has been found. We provide a parallel algorithm for counting distinguishable secondary structures.
<br>
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/figure.jpg" width="300"> <br>
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/parallelRNA_vis.gif" width="600">

## 0. Prerequisites
g++ <br />
openmp <br />

## 1. How to compile and run
```
g++ -o max_pair_omp src/max_pair_omp.cpp -fopenmp 

./max_pair_omp input.txt num_threads
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

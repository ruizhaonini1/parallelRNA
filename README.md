# ParallelRNA: Parallel RNA structure counting algorithm
Secondary structures for nucleic acid sequences, DNA and RNA, are useful abstractions when examining certain behaviors of those molecules. We examine the problem of counting secondary structures compatible with an ordered multiset of sequences. In particular, we address the issue of accounting for indistinguishable secondary structures, for which no fast algorithm has been found. We provide a parallel algorithm for counting distinguishable secondary structures.

## 0. Prerequisites
mpic++ compiler
Boost library with MPI Implementation
gcc/11.2.0

## 1. How to compile and run
countMPI.cpp : Sequencial dynamic programming algorithm <br />
mpic++ -I/path/to/boost -o countMPI src/count_MPI.cpp -L/path/to/boost/stage/lib/ <br />
        &emsp; &emsp; -lboost_mpi -lboost_serialization<br />
mpirun -n numNodes countMPI input_file.txt

<!--## 2. Files-->
<!--The following files are included in this folder, in addition to this readme-->
<!--file, readme.md.-->
<!--<ul>-->
<!--<li>md.c: Main C program</li>-->
<!--<li>md.h: Header file for md.c</li>-->
<!--<li>md.in: Input parameter file (to be redirected to the standard input)</li>-->
<!--</ul>-->
<!--![Screen shot of MD simulation](ScreenShot.png)-->

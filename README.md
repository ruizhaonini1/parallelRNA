# ParallelRNA: Parallel RNA structure counting algorithm
This is a readme file for a parallel RNA structure counting algorithm.
## 0. Prerequisites
mpic++ compiler
Boost library with MPI Implementation
gcc/11.2.0

## 1. How to compile and run
countMPI.cpp : Sequencial dynamic programming algorithm <br />
mpic++ -I/path/to/boost -o countMPI src/countMPI.cpp -L/path/to/boost/stage/lib/ <br />
        &emsp; &emsp; -lboost_mpi -lboost_serialization
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

# ParallelRNA: Parallel RNA structure counting algorithm
This is a readme file for a parallel RNA structure counting algorithm.
## 0 Prerequisites
gcc/8.3.0 <br />
openmp <br />
Boost library

## 1. How to compile and run
count_OMP.cpp : Sequencial dynamic programming algorithm <br />
gcc -o count_OMP -I/path/to/boost src/count_OMP.cpp -fopenmp <br /> 
-lstdc++ -L/path/to/boost/stage/lib/ -lboost_mpi -lboost_serialization

<!--## 2. Files-->
<!--The following files are included in this folder, in addition to this readme-->
<!--file, readme.md.-->
<!--<ul>-->
<!--<li>md.c: Main C program</li>-->
<!--<li>md.h: Header file for md.c</li>-->
<!--<li>md.in: Input parameter file (to be redirected to the standard input)</li>-->
<!--</ul>-->
<!--![Screen shot of MD simulation](ScreenShot.png)-->

# ParallelRNA: Parallel RNA structure counting algorithm
<!-- ![alt text](https://github.com/masarunakajima/parallelRNA/blob/openMP/parallel%20rna.PNG) -->
<!-- ![alt text](https://github.com/masarunakajima/parallelRNA/blob/openMP/total%20fig.jpg) -->
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/figure.jpg" width="300">
<img src="https://github.com/masarunakajima/parallelRNA/blob/openMP/ezgif.com-gif-maker.gif" width="600">
## 0. Prerequisites
gcc/8.3.0 <br />
openmp <br />
Boost library

## 1. How to compile and run
count_OMP.cpp : Sequencial dynamic programming algorithm <br />
gcc -o count_OMP -I/path/to/boost src/count_OMP.cpp -fopenmp <br /> 
  &emsp; &emsp; &emsp;-lstdc++ -L/path/to/boost/stage/lib/ -lboost_mpi -lboost_serialization <br /> <br />
./count_OMP input.txt

<!--## 2. Files-->
<!--The following files are included in this folder, in addition to this readme-->
<!--file, readme.md.-->
<!--<ul>-->
<!--<li>md.c: Main C program</li>-->
<!--<li>md.h: Header file for md.c</li>-->
<!--<li>md.in: Input parameter file (to be redirected to the standard input)</li>-->
<!--</ul>-->
<!--![Screen shot of MD simulation](ScreenShot.png)-->

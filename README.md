# ParallelRNA: Parallel RNA structure counting algorithm
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

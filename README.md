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
count_OMP.cpp : Sequencial dynamic programming algorithm <br />
gcc  &emsp;-o count_OMP  &emsp;-I/<i>path</i>/<i>to</i>/<i>boost</i>  &emsp; src/count_OMP.cpp  &emsp;<b>-fopenmp</b> <br /> 
  &emsp; &emsp; &emsp;<b>-lstdc++</b>  &emsp;-L/<i>path</i>/<i>to</i>/<i>boost</i>/stage/lib/  &emsp;<b>-lboost_mpi</b>  &emsp;<b>-lboost_serialization</b> <br /> <br />
./count_OMP input.txt

## 2. How to compile and run for measurement
count_OMP_adaptive.cpp : Sequencial dynamic programming algorithm <br />
gcc  &emsp;-o count_OMP_adaptive  &emsp;-I/<i>path</i>/<i>to</i>/<i>boost</i>  &emsp; src/count_OMP_adaptive.cpp  &emsp;<b>-fopenmp</b> <br /> 
  &emsp; &emsp; &emsp;<b>-lstdc++</b>  &emsp;-L/<i>path</i>/<i>to</i>/<i>boost</i>/stage/lib/  &emsp;<b>-lboost_mpi</b>  &emsp;<b>-lboost_serialization</b> <br /> <br />
./count_OMP_adaptive input.txt <i>num_threads</i>

<!--## 2. Files-->
<!--The following files are included in this folder, in addition to this readme-->
<!--file, readme.md.-->
<!--<ul>-->
<!--<li>md.c: Main C program</li>-->
<!--<li>md.h: Header file for md.c</li>-->
<!--<li>md.in: Input parameter file (to be redirected to the standard input)</li>-->
<!--</ul>-->
<!--![Screen shot of MD simulation](ScreenShot.png)-->

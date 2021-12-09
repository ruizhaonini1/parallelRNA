#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=run_5000.out
#SBATCH --account=anakano_429
##SBATCH -p qcb

N=5000
MAX_T=8

gcc -fopenmp -lstdc++ -lm -o max_pair src/max_pair.cpp
gcc -fopenmp -lstdc++ -lm -o max_pair_omp src/max_pair_omp.cpp
gcc -fopenmp -lstdc++ -lm -o max_pair_unsup_omp src/max_pair_unsup_omp.cpp
gcc -fopenmp -lstdc++ -lm -o max_pair_omp_tri src/max_pair_omp_tri.cpp


echo "***omp algorithm***"
for t in $( eval echo {1..$MAX_T} )
do
  echo "nthread = ${t}"
  ./max_pair_omp sample_inputs/input_${N}.txt $t
  echo ""
done
echo "**************************"
echo ""

echo "***omp_tri algorithm***"
for t in $( eval echo {1..$MAX_T} )
do
  echo "nthread = ${t}"
  ./max_pair_omp_tri sample_inputs/input_${N}.txt $t
  echo ""
done
echo "**************************"
echo ""


echo "***unsup_omp algorithm***"
for t in $( eval echo {1..$MAX_T} )
do
  echo "nthread = ${t}"
  ./max_pair_unsup_omp sample_inputs/input_${N}.txt $t
  echo ""
done
echo "**************************"

#sequential method
echo "***Sequential algorithm***"
./max_pair sample_inputs/input_${N}.txt
echo "**************************"



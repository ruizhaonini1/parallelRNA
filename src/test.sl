#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=2
#SBATCH --time=00:01:59
#SBATCH --output=test.out
#SBATCH -Aanakano_429

mpirun -np 2 ./count_MPI 

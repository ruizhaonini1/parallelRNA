#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --output=max_pair.out
#SBATCH -A anakano_429

echo '##### input with 10 bases'
time ./src/max_pair ./sample_inputs/input_10.txt
echo ''

echo '##### input with 50 bases'
time ./src/max_pair ./sample_inputs/input_50.txt
echo ''

echo '##### input with 100 bases'
time ./src/max_pair ./sample_inputs/input_100.txt
echo ''

echo '##### input with 500 bases'
time ./src/max_pair ./sample_inputs/input_500.txt
echo ''


echo '##### input with 1000 bases'
time ./src/max_pair ./sample_inputs/input_1000.txt
echo ''

echo '##### input with 5000 bases'
time ./src/max_pair ./sample_inputs/input_5000.txt
echo ''

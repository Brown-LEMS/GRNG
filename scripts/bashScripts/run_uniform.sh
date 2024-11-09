#!/bin/bash
#SBATCH -J u-6D-3L-204k
#SBATCH -o slurms/uniform/slurm-%j.out 
#SBATCH --account=carney-bkimia-condo
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

cd /users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/
module load gcc/8.3

export d="uniform"
export D=6
export N=204800

# run script
/usr/bin/time -v ./testResults               \
                    -D $D                    \
                    -N $N                    \
                    -n 8                     \
                    -r 0.66610 0.38681
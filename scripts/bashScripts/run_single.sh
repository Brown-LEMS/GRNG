#!/bin/bash
#SBATCH -J LA-L7-Test
#SBATCH -o slurms/uniform/slurm-%j.out 
#SBATCH --account=carney-bkimia-condo
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH --time=72:00:00

cd /users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/
module load gcc/8.3

export d="LA"

# run script
/usr/bin/time -v ./testResults               \
                    -d $d                    \
                    -n 16                    \
                    -C false                 \
                    -r 726.71411 270.63030 100.78345 37.53203 13.97703 5.20508
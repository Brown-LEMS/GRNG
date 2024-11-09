#!/bin/bash
#SBATCH -J rayar-LA
#SBATCH -o slurms/uniform/slurm-%j.out 
#SBATCH --account=carney-bkimia-condo
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=16G
#SBATCH --time=48:00:00

cd /users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/
module load gcc/8.3

export d="LA"
# export D=6

# #               0    1    2     3     4     5      6      7      8      9       10      11      12       13       14       15
# export N_Array=(3200 6400 12800 25600 51200 102400 204800 409600 819200 1638400 3276800 6553600 13107200 26214400 52428800 104857600)
# export N=${N_Array[$SLURM_ARRAY_TASK_ID]}
# echo $N

# run script
/usr/bin/time -v ./rayar                     \
                    -d $d                    \
                    -n 16                     \
                    -G true                   

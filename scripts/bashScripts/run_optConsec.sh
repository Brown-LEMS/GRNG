#!/bin/bash
#SBATCH -J OPT-LA-5L
#SBATCH -o slurms/slurm-%j.out 
#SBATCH --account=carney-bkimia-condo
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH --time=48:00:00

cd /users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/
module load gcc/8.3

# /usr/bin/time -v ./optimize_consecutive -d $d -L $begin -N $end -n 8 -r $r1 -s $s1
# /usr/bin/time -v ./optimize_consecutive -d LA       \
#                                         -n 16       \
#                                         -L 1073727  \
#                                         -N 1073727  \
#                                         -r 616.40440 194.70661 61.50291 19.42722 6.13657 \
#                                         -s 50 25 12 8 4
/usr/bin/time -v ./optimize_consecutive -d LA       \
                                        -n 8       \
                                        -L 3200  \
                                        -N 409600  \
                                        -r 700 395.05553 147.78847 55.28699 \
                                        -s 100 50 25 15
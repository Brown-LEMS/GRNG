#!/bin/bash
#SBATCH -J Corel3_OPT
#SBATCH -o slurms/slurm-%j.out 
#SBATCH --account=carney-bkimia-condo
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=8G
#SBATCH --time=48:00:00

cd /users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/
module load gcc/8.3

export d="Corel3"
export begin=4800
export end=19200
export r1=1.49625
export s1=0.005

/usr/bin/time -v ./optimize_consecutive -d $d -L $begin -N $end -n 8 -r $r1 -s $s1

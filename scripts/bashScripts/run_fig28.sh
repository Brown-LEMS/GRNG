#!/bin/bash
#SBATCH -J fig28
#SBATCH -o slurms/fig28/slurm-%j.out 
##SBATCH --account=carney-bkimia-condo
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=48:00:00
##SBATCH --array=0-13

cd /users/cfoste18/data/cfoste18/SISAP/HierarchicalRNG/bin/
module load gcc/8.3

export d="uniform"
export D=2
export N=102400

# export r1_Array=(0.09841 0.0955 0.08128 0.06443 0.05329 0.04744 0.03544 0.03385 0.02795 0.02273 0.01905 0.01596 0.01337 0.01173 0.0099 0.00835 0.00704 0.005947)
export r1_Array=(0.09841 0.0955 0.08128 0.06443 0.05329 0.04744 0.03544 0.03385 0.02795 0.02273 0.01905 0.01596 0.01337 0.01173)
# export r1=${r1_Array[$SLURM_ARRAY_TASK_ID]}
export r1=0.02273

echo $N
echo $r1

# run script
/usr/bin/time -v ./build_search               \
                    -D $D                    \
                    -N $N                    \
                    -n 1                     \
                    -v true                 \
                    -i true                 \
                    -C false                \
                    -r $r1               

#!/bin/bash

#SBATCH --job-name=conway
#SBATCH --output=conway_%j.out
#SBATCH --error=conway_%j.err
#SBATCH --nodelist=node[1755-1757]
#SBATCH --ntasks-per-node=3
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=100

module load openmpi
mpirun -n 9 conways_base life.pgm 1 10 0 0 0 0

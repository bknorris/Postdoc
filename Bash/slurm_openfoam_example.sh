#!/bin/bash

##########################
# Set up submit to ...
#SBATCH -J Test
#SBATCH --time=12:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 24
#SBATCH -o test-%j.out
##########################

nodes=$SLURM_JOB_NUM_NODES            # Number of nodes
cores=24                      # number of cores

# LOAD THE OPENFOAM MODULE

ml load gcc openmpi/2.0.1

ml load openfoam/4.1

sed -i "/numberOfSubdomains/c\\numberOfSubdomains    $(($nodes*$cores));" system/decomposeParDict

# Decompose solution (serial)
echo "DECOMPOSE MESH WITH decomposePar..."
decomposePar -cellDist > log.decomposeParDict 2>&1

# Run the solution (parallel)
echo "RUNNING THE SOLVER WITH fireFoam..."
mpirun -np $(($nodes*$cores)) fireFoam -parallel > log.fireFoam 2>&1
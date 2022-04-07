#!/bin/bash
#SBATCH --job-name=bknOF  # name your job
#SBATCH --nodes=1  # number of nodes
#SBATCH --ntasks=2  # number of CPUs
#SBATCH --ntasks-per-node=2  # tasks per node
#SBATCH --account=pcmsc  # account code
#SBATCH --time=0-12:00:00  # job time D-HH:MM:SS
#SBATCH --mail-type=ALL  # BEGIN, END, FAIL, REQUEUE, and ALL
#SBATCH --mail-user=bknorris@usgs.gov  # email address
#SBATCH -o testrun-%j.out  # name of output file

source /home/bknorris/OpenFOAM/OpenFOAM-v1912/etc/bashrc
#------------------------------------------------------------------------------
cp -r 0.orig/ 0/

echo "Decomposing mesh..."
blockMesh > log.blockMesh 2>&1

decomposePar -force > log.decomposeParDict 2>&1

echo "Setting fields..."
mpirun -np 2 setFields -parallel > log.setFields 2>&1

echo "Running interFoam solver..."
mpirun -np 2 interFoam -parallel > log.interFoam 2>&1

#------------------------------------------------------------------------------
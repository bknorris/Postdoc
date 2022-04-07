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
MODEL_PATH='/caldera/projects/usgs/hazards/pcmsc/openfoam/'
SCENARIO='Scenario_4'

cd Model
rm -rf log.* 0/* constant/polyMesh/*
cp -r 0.org/* 0/
cp -r $MODEL_PATH$SCENARIO/Extrude/constant/polyMesh/* constant/polyMesh #copy mesh
cp -r $MODEL_PATH/3D_Models/Domains_sampleVolume_Length_4_V2.stl constant/triSurface/bathymetry_sample.stl
#------------------------------------------------------------------------------
echo "Decomposing mesh..."
decomposePar -force > log.decomposeParDict 2>&1

echo "Setting fields..."
mpirun -np 2 setFields -parallel > log.setFields 2>&1

echo "Running interFoam solver..."
mpirun -np 2 interFoam -parallel > log.interFoam 2>&1
#------------------------------------------------------------------------------
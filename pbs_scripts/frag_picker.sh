#!/bin/bash
#SBATCH --account 
#SBATCH --partition 
#SBATCH --nodes 1
#SBATCH --ntasks 16
##SBATCH --job-name fragments
#SBATCH --mem=16g
#SBATCH --time 2:00:00
##SBATCH --output ./log/%x-%j.log

module load PE-gnu/2.0

$ROSETTA_BIN/fragment_picker.default.linuxgccrelease @$1
#mpirun -n 16 $ROSETTA_BIN/fragment_picker.mpi.linuxgccrelease @$1

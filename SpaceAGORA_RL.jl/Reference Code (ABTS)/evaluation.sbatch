#!/bin/bash
#SBATCH --nodes=1
#SBATCH --tasks=40
#SBATCH --time=4:00:00
#SBATCH --job-name="Evaluation"
#SBATCH --partition=secondary
#SBATCH --output=dqn.o%j
#SBATCH --mail-user=gfalcon2@illinois.edu
#SBATCH --error=dqn.e%j
#SBATCH --mail-type=ALL

cd ~/scratch/DRLAerobraking/

module load openmpi/4.1.0-gcc-7.2.0
module load python/3
module load cuda

python3 Evaluation.py

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=40
#SBATCH --time=72:00:00
#SBATCH --partition=putnamrg
#SBATCH --job-name="DQN"
#SBATCH --output=dqn.o%j
#SBATCH --mail-user=gfalcon2@illinois.edu
#SBATCH --error=dqn.e%j
#SBATCH --mail-type=ALL
#SBATCH --gres=gpu:V100:1

cd ~/scratch/DRLAerobraking/

module load openmpi/4.1.0-gcc-7.2.0
module load python/3
module load cuda

srun python3 Main_DRL.py

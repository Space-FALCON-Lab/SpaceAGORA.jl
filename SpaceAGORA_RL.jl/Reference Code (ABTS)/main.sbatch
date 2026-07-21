#!/bin/bash


cd ~/scratch/DRLAerobraking/

JOB_01=`sbatch run_dqn.sbatc |cut -f 4 -d " "`
JOB_02=`sbatch --dependency=afterany:$JOB_01 evaluation.sbatch |cut -f 4 -d " "`

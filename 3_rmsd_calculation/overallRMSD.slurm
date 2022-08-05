#!/bin/bash
#SBATCH --job-name=overallRMSD_1
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH -n 40
#SBATCH --ntasks-per-node=40
#SBATCH --output=task_file/%x_%j.out
#SBATCH --error=task_file/%x_%j.err

module load miniconda3
source activate pymol
cd /lustre/home/acct-stu/stu/deepocean2


python ${SLURM_JOB_NAME}.py


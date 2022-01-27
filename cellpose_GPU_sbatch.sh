#!/bin/bash
#SBATCH --partition gpu
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --gres gpu:1
#SBATCH --mem 40G
#SBATCH --time 23:00:00
#SBATCH --job-name 2101
#SBATCH --array=1-24

python Register_segment_for_cellpose_tracking.py AU02101 $SLURM_ARRAY_TASK_ID

#!/bin/bash
#SBATCH --partition gpu
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --gres gpu:1
#SBATCH --mem 40G
#SBATCH --time 23:00:00
#SBATCH --job-name 2001110
#SBATCH --array=1-10

python /home/exacloud/gscratch/HeiserLab/software/image_analysis_pipelines/AU565/python/Process_LCI_nuclear.py AU02001 $SLURM_ARRAY_TASK_ID

#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --mem 40G
#SBATCH --time 4:00:00
#SBATCH --job-name 2100801
#SBATCH --array=1-24

python /home/exacloud/gscratch/HeiserLab/software/image_analysis_pipelines/AU565/python/Process_LCI_nuclear.py $SLURM_JOB_NAME $SLURM_ARRAY_TASK_ID

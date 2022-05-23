#!/bin/bash
#SBATCH --partition gpu
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --gres gpu:1
#SBATCH --mem 40G
#SBATCH --time 4:00:00
#SBATCH --job-name HC00801
#SBATCH --array=2-6

python /home/exacloud/gscratch/HeiserLab/software/image_analysis_pipelines/AU565/python/Process_LCI_cyto2_LC2.py HC00801 $SLURM_ARRAY_TASK_ID

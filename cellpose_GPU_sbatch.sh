#!/bin/bash
#SBATCH --partition gpu
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
# SBATCH --gpus-per-task 1
#SBATCH --gres gpu:v100:1
#SBATCH --mem 10G
#SBATCH --time 23:00:00
#SBATCH --job-name CP_AU_B

srun python AU565_Live_cell_image_segmentation.py AU01401 B1
srun python AU565_Live_cell_image_segmentation.py AU01401 B2
srun python AU565_Live_cell_image_segmentation.py AU01401 B3
srun python AU565_Live_cell_image_segmentation.py AU01401 B4
srun python AU565_Live_cell_image_segmentation.py AU01401 B5
srun python AU565_Live_cell_image_segmentation.py AU01401 B6

#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --mem 10G
#SBATCH --time 23:00:00
#SBATCH --job-name OME_Com

srun python AU565_create_companions.py AU000601
srun python AU565_create_companions.py AU000602
srun python AU565_create_companions.py AU000701
srun python AU565_create_companions.py AU000702
srun python AU565_create_companions.py AU000801
srun python AU565_create_companions.py AU000802
srun python AU565_create_companions.py AU000901
srun python AU565_create_companions.py AU000902
srun python AU565_create_companions.py AU001001
srun python AU565_create_companions.py AU001002
srun python AU565_create_companions.py AU001101
srun python AU565_create_companions.py AU001102


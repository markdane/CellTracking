#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --mem 10G
#SBATCH --time 23:00:00
#SBATCH --job-name cp151

srun python Live_cell_imaging_file_copy.py AU_I_L_015_01_1 AU01501

"AU01401","AU01402","AU01501","AU01502","AU01601","AU01602","AU01701","AU01702","AU01801","AU01802","AU01901","AU01902","AU02001","AU02002","AU02101"

srun -c 8 -J AI1401 -o AI1401_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01401 &
srun -c 8 -J AI1402 -o AI1402_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01402 &
srun -c 8 -J AI1501 -o AI1501_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI1502 -o AI1502_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI1601 -o AI1601_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01601 &
srun -c 8 -J AI1602 -o AI1602_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01602 &
srun -c 8 -J AI1701 -o AI0701_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI1702 -o AI0702_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI1801 -o AI001_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI1802 -o AI1802_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI1901 -o AI1901_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI1902 -o AI1902_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI2001 -o AI2001_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI2002 -o AI2002_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI2101 -o AI2101_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01601 &

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
srun -c 8 -J AI1501 -o AI1501_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI1802 -o AI1802_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI1901 -o AI1901_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI1902 -o AI1902_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI2001 -o AI2001_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 &
srun -c 8 -J AI2002 -o AI2002_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01502 &
srun -c 8 -J AI2101 -o AI2101_out -t 23:00:00 python Live_cell_imaging_file_copy.py AU01601 &
srun -c 8 -J C2301 -o C2301_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU02301 &
srun -c 8 -J C2401 -o C2401_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU02401 AU02401 &
srun -c 8 -J C2501 -o C2501_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU02501 AU02501 &
srun -c 8 -J C1501 -o C1501_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU01501 AU01501 &
srun -c 8 -J C2801 -o C2801_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU02801 AU02801 &
srun -c 8 -J C2901 -o C2901_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU02901 AU02901 &
srun -c 8 -J C3001 -o C3001_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU03001 AU03001 &
srun -c 8 -J C3101 -o C3101_out.txt -t 23:00:00 python Live_cell_imaging_file_copy.py AU03101 AU03101 &




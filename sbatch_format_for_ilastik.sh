#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 2
#SBATCH --mem 10G
#SBATCH --time 23:00:00
#SBATCH --job-name cp151


srun -c 8 -J FI1401 -o FI1401_out -t 23:00:00 python Format_for_ilastik.py AU01401 &
srun -c 8 -J FI1402 -o FI1402_out -t 23:00:00 python Format_for_ilastik.py AU01402 &
srun -c 8 -J FI1501 -o FI1501_out -t 23:00:00 python Format_for_ilastik.py AU01501 &
srun -c 8 -J FI1502 -o FI1502_out -t 23:00:00 python Format_for_ilastik.py AU01502 &
srun -c 8 -J FI1601 -o FI1601_out -t 23:00:00 python Format_for_ilastik.py AU01601 &
srun -c 8 -J FI1602 -o FI1602_out -t 23:00:00 python Format_for_ilastik.py AU01602 &
srun -c 8 -J FI1701 -o FI1701_out -t 23:00:00 python Format_for_ilastik.py AU01701 &
srun -c 8 -J FI1702 -o FI1702_out -t 23:00:00 python Format_for_ilastik.py AU01702 &
srun -c 8 -J FI1801 -o FI1801_out -t 23:00:00 python Format_for_ilastik.py AU01801 &
srun -c 8 -J FI1802 -o FI1802_out -t 23:00:00 python Format_for_ilastik.py AU01802 &
srun -c 8 -J FI1901 -o FI1901_out -t 23:00:00 python Format_for_ilastik.py AU01901 &
srun -c 8 -J FI1902 -o FI1902_out -t 23:00:00 python Format_for_ilastik.py AU01902 &
srun -c 8 -J FI2001 -o FI2001_out -t 23:00:00 python Format_for_ilastik.py AU02001 &
srun -c 8 -J FI2002 -o FI2002_out -t 23:00:00 python Format_for_ilastik.py AU02002 &
srun -c 8 -J FI2101 -o FI2101_out -t 23:00:00 python Format_for_ilastik.py AU02101 &

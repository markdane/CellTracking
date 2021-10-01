#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G
#SBATCH --time 23:00:00

srun -c 4 -o FI0601.txt --job-name FI0601 --time 23:00:00 python Format_for_ilastik.py AU00601 &
srun -c 4 -o FI0602.txt --job-name FI0602 --time 23:00:00 python Format_for_ilastik.py AU00602 &
srun -c 4 -o FI0701.txt --job-name FI0701 --time 23:00:00 python Format_for_ilastik.py AU00701 &
srun -c 4 -o FI0702.txt --job-name FI0702 --time 23:00:00 python Format_for_ilastik.py AU00702 &
srun -c 4 -o FI0801.txt --job-name FI0801 --time 23:00:00 python Format_for_ilastik.py AU00801 &
srun -c 4 -o FI0802.txt --job-name FI0802 --time 23:00:00 python Format_for_ilastik.py AU00802 &
srun -c 4 -o FI0901.txt --job-name FI0901 --time 23:00:00 python Format_for_ilastik.py AU00901 &
srun -c 4 -o FI0902.txt --job-name FI0902 --time 23:00:00 python Format_for_ilastik.py AU00902 &
srun -c 4 -o FI1001.txt --job-name FI1001 --time 23:00:00 python Format_for_ilastik.py AU01001 &
srun -c 4 -o FI1002.txt --job-name FI1002 --time 23:00:00 python Format_for_ilastik.py AU01002 &
srun -c 4 -o FI1101.txt --job-name FI1101 --time 23:00:00 python Format_for_ilastik.py AU01101 &
srun -c 4 -o FI1102.txt --job-name FI1102 --time 23:00:00 python Format_for_ilastik.py AU01102 &
srun -c 4 -o FI1102.txt --job-name FI1102 --time 23:00:00 python Format_for_ilastik.py AU01102 &
srun -c 4 -o FI2301.txt --job-name FI2301 --time 23:00:00 python Format_for_ilastik.py AU02301 &

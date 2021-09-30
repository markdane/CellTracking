#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 8
#SBATCH --mem 10G
#SBATCH --time 23:00:00

python generate_apply_ilastik_jobs.py AU00601 &
python generate_apply_ilastik_jobs.py AU00602 &
python generate_apply_ilastik_jobs.py AU00701 &
python generate_apply_ilastik_jobs.py AU00702 &
python generate_apply_ilastik_jobs.py AU00801 &
python generate_apply_ilastik_jobs.py AU00802 &
python generate_apply_ilastik_jobs.py AU00901 &
python generate_apply_ilastik_jobs.py AU00902 &
python generate_apply_ilastik_jobs.py AU01001 &
python generate_apply_ilastik_jobs.py AU01002 &
python generate_apply_ilastik_jobs.py AU01101 &
python generate_apply_ilastik_jobs.py AU01102 &


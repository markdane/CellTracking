#!/bin/bash
#SBATCH --partition exacloud
#SBATCH --account HeiserLab
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G
#SBATCH --time 23:00:00

python generate_ilastik_jobs.py AU00601 &
python generate_ilastik_jobs.py AU00602 &
python generate_ilastik_jobs.py AU00701 &
python generate_ilastik_jobs.py AU00702 &
python generate_ilastik_jobs.py AU00801 &
python generate_ilastik_jobs.py AU00802 &
python generate_ilastik_jobs.py AU00901 &
python generate_ilastik_jobs.py AU00902 &
python generate_ilastik_jobs.py AU01001 &
python generate_ilastik_jobs.py AU01002 &
python generate_ilastik_jobs.py AU01101 &
python generate_ilastik_jobs.py AU01102 &
python generate_ilastik_jobs.py AU01401 &
python generate_ilastik_jobs.py AU01501 &
python generate_ilastik_jobs.py AU01502 &
python generate_ilastik_jobs.py AU99999 &
python generate_ilastik_jobs.py AU02301 &
python generate_ilastik_jobs.py AU02401 &
python generate_ilastik_jobs.py AU02501 &
python generate_ilastik_jobs.py AU02901 &
python generate_ilastik_jobs.py AU03001 &


#!/bin/bash
#SBATCH --job-name=lsmethod-scan
#SBATCH --output=results/numeric_checks/slurm-%j.out
#SBATCH --error=results/numeric_checks/slurm-%j.err
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1

source ~/.bashrc
mamba activate lsmethod

python scripts/scan_random_points.py --npoints 1000 --eps 0.08 --out results/numeric_checks/scan_${SLURM_JOB_ID}.csv

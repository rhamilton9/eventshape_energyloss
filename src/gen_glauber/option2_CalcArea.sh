#!/bin/bash
#SBATCH --job-name=calculate_area
#SBATCH --time=4:30:00
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=100G
#SBATCH --output="log/slurm%A_%a.out"
#SBATCH --mail-type=ALL

cd /home/rjh78/project/glauber

[ ! -d compiled ] && mkdir compiled

module load ROOT/6.26.10-foss-2022b

root -l -q "calc_area.c(${SLURM_ARRAY_TASK_ID},${SLURM_ARRAY_TASK_COUNT}, true)"

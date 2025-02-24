#!/bin/bash
#SBATCH --job-name=glauber_gen
#SBATCH --time=3:30:00
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=90G
#SBATCH --output="log/slurm%A_%a.out"
#SBATCH --mail-type=ALL

nevent=100000

cd /home/rjh78/project/glauber

[ ! -d out ] && mkdir out
[ ! -d log ] && mkdir log
[ ! -d compiled ] && mkdir compiled

module load ROOT/6.26.10-foss-2022b

root -l -q "runGlauber_local.c($nevent, ${SLURM_ARRAY_TASK_ID}, ${SLURM_ARRAY_TASK_COUNT})"
root -l -q "calc_area.c(${SLURM_ARRAY_TASK_ID},${SLURM_ARRAY_TASK_COUNT}, true)"

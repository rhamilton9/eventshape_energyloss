#!/bin/bash
#SBATCH --job-name=glauber_gen
#SBATCH --time=1:00:00
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=10G
#SBATCH --output="log/slurm%A_%a.out"
#SBATCH --mail-type=ALL

nevent=5000000

cd /home/rjh78/project/glauber

[ ! -d out ] && mkdir out
[ ! -d log ] && mkdir log

module load ROOT/6.26.10-foss-2022b

root -l -q "runGlauber_local.c($nevent, ${SLURM_ARRAY_TASK_ID}, ${SLURM_ARRAY_TASK_COUNT})"

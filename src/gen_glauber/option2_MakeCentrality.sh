#!/bin/bash
#SBATCH --job-name=make_centrality
#SBATCH --time=0:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --output="log/slurm%A.out"
#SBATCH --mail-type=ALL

module load ROOT/6.26.10-foss-2022b

root -l -q "mergetree.c(1)"
root -l -q make_centrality.c
#rm out/glauber_*thread*.root

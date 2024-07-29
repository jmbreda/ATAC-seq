#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 10G
#SBATCH --time 01:00:00
#SBATCH --output output.txt
#SBATCH --error error.txt
#SBATCH --qos serial

module load gcc/11.3.0
module load py-deeptools
bamCoverage -b results/mapping/SRR8119853_sorted.bam -o results/mapping/SRR8119853_coverage.bw

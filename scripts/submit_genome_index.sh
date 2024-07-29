#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1
#SBATCH --mem 200G
#SBATCH --time 12:00:00
#SBATCH --output output.txt
#SBATCH --error error.txt
#SBATCH --qos serial

bwa-mem2 index /scratch/jbreda/genome/mm10/Mus_musculus.GRCm38.dna.primary_assembly.fa


#!/bin/bash
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 8
#SBATCH --mem 50G
#SBATCH --time 24:00:00
#SBATCH --output stdout/download_fastq
#SBATCH --error stderr/download_fastq
#SBATCH --qos serial

module load sratoolkit

infile='resources/SRR_Acc_List.txt'
mapfile -t Accession_list < $infile
outdir='resources/fastq'

prefetch --progress "${Accession_list[@]}"
fasterq-dump  --split-files --outdir $outdir --split-files --mem 50GB --progress "${Accession_list[@]}"



#!/bin/bash

module load sratoolkit
Accession_list='resources/SRR_Acc_List.txt'

while IFS= read -r line
do
  echo "$line"
  fasterq-dump "$line" --outdir /scratch/jbreda/M-F_Liver_ATAC-seq/resources/reads -p
done < $Accession_list
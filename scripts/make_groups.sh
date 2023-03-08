#!/bin/bash
infile=../resources/SraRunTable.txt
outfile=../config/tissue_sex_groups.txt
SEX=(male female)
printf "TISSUE_SEX_SAMPLES:\n" > $outfile
for tissue in $(cut -f 37 -d',' $infile |sort -u)
do
    printf "\t- %s:\n" $tissue >> $outfile
    for sex in ${SEX[@]}
    do
        printf "\t\t- %s: " $sex >> $outfile
        awk -F, -v s=$sex -v t=$tissue '($34==s && $37==t) {print $1}' ../resources/SraRunTable.txt | tr '\n' ' ' >> $outfile
        printf "\n" >> $outfile
    done
done

outfile=../config/tissue_groups.txt
printf "TISSUE_SAMPLES:\n" > $outfile
for tissue in $(cut -f 37 -d',' $infile |sort -u)
do
    printf "\t- %s: " $tissue >> $outfile
    awk -F, -v t=$tissue '($37==t) {print $1}' ../resources/SraRunTable.txt | tr '\n' ' ' >> $outfile
    printf "\n" >> $outfile
done
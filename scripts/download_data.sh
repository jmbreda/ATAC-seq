#!/bin/bash

Accession_list='../resources/SRR_Acc_List.txt'
for dataset in $(cat $Accession_list)
do
	echo $dataset
	fasterq-dump $dataset --outdir /bigdata/jbreda/ATACseq/resources -p
done
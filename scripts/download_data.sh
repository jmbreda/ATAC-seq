#!/bin/bash

DATASETS=(SRR8119838 SRR8119839 SRR8119852 SRR8119853)
for dataset in ${DATASETS[@]}
do
	echo $dataset
	fasterq-dump $dataset --outdir resources/ -p
done
	 	

#!/bin/bash

if [ $# -lt 2 ]
then
    echo not enough arguments: \($#\)
    echo samtools_genome_info infile outfile
else
    infile=$1
    outfile=$2
    samtools view -H $infile| perl -ne 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print $1,"\t",$2,"\n"}' > $outfile
fi
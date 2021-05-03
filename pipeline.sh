#!/bin/bash
n=$#
if [ $n -lt 5 ]; then
    echo "Following pars have to be specified: reference genome, sequences 1, sequences 2, plot title and output file."
    echo "Optional args: -ni - the reference will not be indexed."
    exit
fi

ref=$1
pieces="$2 $3"
title=$4
output=$5

echo "REF: $ref"
echo "PIE: $pieces"
echo "TIT: $title"
echo "OUT: $output"

if [[ $@ != *"-ni"* ]]; then
    echo "Indexing..."
    bwa index $ref
else
    echo "Indexing skipped..."
fi

echo "Aligning..."
bwa mem $ref $2 $3 > results-wt/alignment.bam

echo "Sorting..."
samtools sort results-wt/alignment.bam -o results-wt/sorted.bam

echo "Calculating depth..."
samtools depth results-wt/sorted.bam > results-wt/output.txt

echo "Generating graph..."
python3 plot.py results-wt/output.txt "$title" "$output" -s

echo "DONE! :-)"

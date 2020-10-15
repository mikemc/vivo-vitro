#!/bin/bash
#SBATCH -o map-sample-%j.out
set -euo pipefail

## Args
# Path with burst database
burst_dir=$1
# Fastq file with reads to map
fastq_file=$2
# Output directory
out_dir=$3
# Path to blast6-to-tibble.R
compress_script=$4
# Number of cpus to use
nproc=$5

## Paths
# Fasta file to serve as BURST input
fasta_file=$out_dir/$(basename $fastq_file).fasta
# BURST output file
burst_out=$out_dir/$(basename $fastq_file)-burst.tsv

## Convert Fastq to Fasta; replace spaces in sequence headers
seqtk seq -A $fastq_file | sed "s/ /_/g" > $fasta_file

## Map with burst
burst_linux_DB12 -r ${burst_dir}/db.edx -a ${burst_dir}/db.acx \
    -q ${fasta_file} -o ${burst_out} \
    --mode ALLPATHS --threads $nproc --id 0.97 --forwardreverse
rm $fasta_file

## Compress burst out into a tibble Rds
Rscript --no-init-file $compress_script $burst_out

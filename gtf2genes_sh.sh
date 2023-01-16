#!/bin/bash
gtf=$1
fasta=$2
name=$3
gff2bed < "$gtf" > "$name".bed
bedtools getfasta -fi "$fasta" -bed "$name".bed  -fo "$name"_gene.fa

#!/bin/bash
# requirements needed: bedtools https://bedtools.readthedocs.io/en/latest/

# Motif enrichment relative to crosslink site
# This pipeline will flank 300bps upstream and downstream each crosslink site from the input BED file and plot the motif enrichment from the second input file in txt format.

crosslinks_rep1=$1	# RBP corrlinks sites in BED format of replciate 1
crosslinks_rep2=$2	# .. replicate 2
crosslinks_rep3=$3	# .. replicate 3
motif=$4		# list of motifs in TXT format, where each motif is listen in seperate line
output=$5		# output file bane
GENOME_FASTA="../hg38.fa" # input to the reference genome fasta file

# flank the region 300 bps upstream and downstream
python flankBEDpositionsCustom.py ${crosslinks_rep1}-3SS.bed ${crosslinks_rep1}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${crosslinks_rep2}-3SS.bed ${crosslinks_rep2}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${crosslinks_rep3}-3SS.bed ${crosslinks_rep3}-3SS-flanked300.bed 300 300

# get fasta 
bedtools getfasta -s -name -fi ${GENOME_FASTA} -bed ${crosslinks_rep1}-3SS-flanked300.bed -fo ${crosslinks_rep1}-3SS-flanked300.fasta
bedtools getfasta -s -name -fi ${GENOME_FASTA} -bed ${crosslinks_rep2}-3SS-flanked300.bed -fo ${crosslinks_rep2}-3SS-flanked300.fasta
bedtools getfasta -s -name -fi ${GENOME_FASTA} -bed ${crosslinks_rep3}-3SS-flanked300.bed -fo ${crosslinks_rep3}-3SS-flanked300.fasta

# get motif density
python k-mer_coverage-norm.py ${crosslinks_rep1}-3SS-flanked300.fasta ${motifs} ${crosslinks_rep1}-3SS-flanked300.motif.density.txt
python k-mer_coverage-norm.py ${crosslinks_rep2}-3SS-flanked300.fasta ${motifs} ${crosslinks_rep2}-3SS-flanked300.motif.density.txt
python k-mer_coverage-norm.py ${crosslinks_rep3}-3SS-flanked300.fasta ${motifs} ${crosslinks_rep3}-3SS-flanked300.motif.density.txt


# draw the motif density relative to crosslink site
Rscript RNAmap-motif-enrichment.R ${crosslinks_rep1}-3SS-flanked300.motif.density.txt ${crosslinks_rep2}-3SS-flanked300.motif.density.txt ${crosslinks_rep3}-3SS-flanked300.motif.density.txt ${name}

#!/bin/bash

# RNA splicing map 
# requirements needed: bedtools https://bedtools.readthedocs.io/en/latest/

crosslinks=$1	# RBP corrlinks sites in BED format
enhanced=$2		# coordinates of enhanced exons in BED format
silenced=$3		# coordinates of silenced exons in BED format
control=$4		# coordinates of control exons in BED format
name=$5			# outdput name

mkdir ${name}

# get splice sites
python getSpliceSites-exonBED.py ${enhanced} ${enhanced}
python getSpliceSites-exonBED.py ${silenced} ${silenced}
python getSpliceSites-exonBED.py ${control} ${control}

# flank region
python flankBEDpositionsCustom.py ${enhanced}-3SS.bed ${enhanced}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${silenced}-3SS.bed ${silenced}-3SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${control}-3SS.bed ${control}-3SS-flanked300.bed 300 300

python flankBEDpositionsCustom.py ${enhanced}-5SS.bed ${enhanced}-5SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${silenced}-5SS.bed ${silenced}-5SS-flanked300.bed 300 300
python flankBEDpositionsCustom.py ${control}-5SS.bed ${control}-5SS-flanked300.bed 300 300

# sort
sort -k1,1 -k2,2n -k6,6 ${enhanced}-3SS-flanked300.bed > ${enhanced}-3SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${silenced}-3SS-flanked300.bed > ${silenced}-3SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${control}-3SS-flanked300.bed > ${control}-3SS-flanked300-sorted.bed

sort -k1,1 -k2,2n -k6,6 ${enhanced}-5SS-flanked300.bed > ${enhanced}-5SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${silenced}-5SS-flanked300.bed > ${silenced}-5SS-flanked300-sorted.bed
sort -k1,1 -k2,2n -k6,6 ${control}-5SS-flanked300.bed > ${control}-5SS-flanked300-sorted.bed

# get crosslink coverage 
bedtools coverage -sorted -s -b ${crosslinks} -a ${enhanced}-3SS-flanked300-sorted.bed -d > ${enhanced}-3SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${crosslinks} -a ${silenced}-3SS-flanked300-sorted.bed -d > ${silenced}-3SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${crosslinks} -a ${control}-3SS-flanked300-sorted.bed -d > ${control}-3SS-flanked300-sorted-crosslink-coverage.bed

bedtools coverage -sorted -s -b ${crosslinks} -a ${enhanced}-5SS-flanked300-sorted.bed -d > ${enhanced}-5SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${crosslinks} -a ${silenced}-5SS-flanked300-sorted.bed -d > ${silenced}-5SS-flanked300-sorted-crosslink-coverage.bed
bedtools coverage -sorted -s -b ${crosslinks} -a ${control}-5SS-flanked300-sorted.bed -d > ${control}-5SS-flanked300-sorted-crosslink-coverage.bed

# plot RNA map 
Rscript RNAmaps-cDNAstarts-enrichment.R ${enhanced}-3SS-flanked300-sorted-crosslink-coverage.bed ${enhanced}-5SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-3SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-5SS-flanked300-sorted-crosslink-coverage.bed ${control}-3SS-flanked300-sorted-crosslink-coverage.bed ${control}-5SS-flanked300-sorted-crosslink-coverage.bed ${name}

# move results into folder
mv ${name}.pdf ${name}
mv ${name}.enhanced.3SS.tab ${name}
mv ${name}.enhanced.5SS.tab ${name}
mv ${name}.silenced.3SS.tab ${name}
mv ${name}.silenced.5SS.tab ${name}

# clean
rm ${enhanced}-3SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-3SS-flanked300-sorted-crosslink-coverage.bed ${control}-3SS-flanked300-sorted-crosslink-coverage.bed
rm ${enhanced}-5SS-flanked300-sorted-crosslink-coverage.bed ${silenced}-5SS-flanked300-sorted-crosslink-coverage.bed ${control}-5SS-flanked300-sorted-crosslink-coverage.bed
rm ${enhanced}-3SS-flanked300-sorted.bed ${silenced}-3SS-flanked300-sorted.bed ${control}-3SS-flanked300-sorted.bed
rm ${enhanced}-5SS-flanked300-sorted.bed ${silenced}-5SS-flanked300-sorted.bed ${control}-5SS-flanked300-sorted.bed
rm ${enhanced}-3SS-flanked300.bed ${silenced}-3SS-flanked300.bed ${control}-3SS-flanked300.bed
rm ${enhanced}-5SS-flanked300.bed ${silenced}-5SS-flanked300.bed ${control}-5SS-flanked300.bed 
rm ${enhanced}-3SS.bed ${silenced}-3SS.bed ${control}-3SS.bed
rm ${enhanced}-5SS.bed ${silenced}-5SS.bed ${control}-5SS.bed


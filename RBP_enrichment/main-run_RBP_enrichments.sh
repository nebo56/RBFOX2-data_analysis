#!/bin/bash

# script will return percentage of crosslinks per target region.
# INPUT folder = FILE_PATH: must contain sorted crosslink BED file, where each line is one crosslink position.

FILE_PATH=/*.bed  # path to a folder with crosslink positions in BED format
INPUT_TARGET=$1   # target regions in BED format
INPUT_CONTROL=$2  # control regions in BED format
OUTPUT=$3 #text file for results

# create headers 
echo  "RBP  norm.enrichment" > $OUTPUT.region.txt

for CL in $FILE_PATH
do
  filename=${CL##*/}     #remove the path from the filename
  RBP=$(echo $filename| cut -d'_' -f 1) #get RBP name out of the file name vy removing everything after "_"

  # controls
  bedtools coverage -s -b $CL -a gencode.v27.primary_assembly.annotation.exon.number.filtered.random10000-flanked-both_SS.bed -d > control.${RBP}
  
  # selected region - regulated exons
  bedtools coverage -s -b $CL -a deltaPSI_tumor_normal_lich_analysis_AS.fdr_0.05.DPSI_5.strand-allSS.flanked.bed -d > region.${RBP}
  
  # get crosslink enrichment for each peak - normalised average crosslink number from surrounding control region
  Rscript get_enrichment-mean.R control.${RBP} region.${RBP} ${RBP} $OUTPUT.region.txt

  # clean
  rm control.${RBP} 
  rm region.${RBP}
done

# plot histograms
Rscript plot_RBP_enrichment.R $OUTPUT.region.txt $OUTPUT.region.pdf

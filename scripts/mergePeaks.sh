#!/usr/bin/env bash
##replace with your own merged output path
awk -v OFS='\t' ' { print $1, $2, $3} ' /nfs/lab/katha/multiomics/scripts/call_peaks/*.merged* > allPeaksAnno.bed
sort -k1,1 -k2,2n allPeaksAnno.bed > allPeaksAnno_sorted.bed
bedtools merge -i  allPeaksAnno_sorted.bed > mergedPeak.txt

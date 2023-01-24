#!/usr/bin/env bash

while getopts c:t:b:o: flag
do
    case "${flag}" in
        c) cells_temp=${OPTARG};;
        t) tagAlign=${OPTARG};;
        b) barcodes_temp=${OPTARG};;
        o) peakCallDir=${OPTARG};;
    esac
done

#create array of cells/leiden
mapfile -t cells < $cells_temp

#print out all the cells/leiden
echo "${cells[*]}"

#create array of barcode paths
mapfile -t barcodes < $barcodes_temp

#getting length of $cells
num_of_cells=${#cells[@]}

N=10

for (( i=0; i<${num_of_cells}; i++ )); do
	((j=j%N)); ((j++==0)) && wait
    (	
    #setting current cell type/leiden   
    cell=${cells[$i]}
	
    echo "working on leiden ${cell}"
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "$dt"

    barcodeFile=${barcodes[$i]}
    outfile="${peakCallDir}split.${cell}.tagAlign"

    echo "zgrep -F -f ${barcodeFile} ${tagAlign} > ${outfile}"

    zgrep -F -f ${barcodeFile} ${tagAlign} > ${outfile}
    echo "done working on ${cell}"
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "done splitting tagAlign for ${cell}      $dt"
    ) &
done


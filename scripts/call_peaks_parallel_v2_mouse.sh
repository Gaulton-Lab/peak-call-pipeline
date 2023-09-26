#!/usr/bin/env bash
while getopts c:t:b:o:g: flag
do
    case "${flag}" in
        c) cells_temp=${OPTARG};;
        t) splitTag_temp=${OPTARG};;
        b) barcodes_temp=${OPTARG};;
        o) peakCallDir=${OPTARG};;
	g) GSIZE_FILE=${OPTARG};;
    esac
done

#create array of cells/leiden
mapfile -t cells < $cells_temp

#print out all the cells/leiden
echo "${cells[*]}"

#create array of splitTag paths
mapfile -t splitTags < $splitTag_temp

#create array of barcode paths
mapfile -t barcodes < $barcodes_temp

#setting number for parallelization
N=6

#getting length of $cells
num_of_cells=${#cells[@]}
echo "number of cell types is $num_of_cells"

#make dir for log files
log_dir=${peakCallDir}"log_files"
mkdir $log_dir

for (( i=0; i<${num_of_cells}; i++ )); do
   (
	((j=j%N)); ((j++==0)) && wait
    #setting current cell type/leiden
    cell=${cells[$i]}

    #make log_file
    log_file=${log_dir}"/"${cell}".log_file.txt"
    touch $log_file

    echo "current cell type is $cell" >> $log_file
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "working on    ${cell}       $dt"  >> $log_file
    echo "current split tagAlign file is ${splitTags[$i]}.gz"  >> $log_file

    ## 1. Call Peaks
    #input: zipped split tagAlign file
    #       path to output directory
    #       current cell
    #output: .narrowPeak
    #        .xls
    #        .summits.bed

    tagalignFilegz="${splitTags[$i]}.gz"
    echo "${cell} macs2 call"  >> $log_file
    macs2 callpeak -t ${tagalignFilegz} --outdir ${peakCallDir} -n ${cell} -q 0.05 --nomodel --keep-dup all -g mm

    ## 2 . count reads
    #input: un-zipped split tagAlign
    #output: fileSize file containing amount of reads in split tagAlign

    tagalignFile=${splitTags[$i]}
    fileSizeFile="${peakCallDir}fileSize.${cell}.totalReadCount.txt"

    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "${cell} wc -l       $dt"  >> $log_file
    sizeA=$(wc -l ${tagalignFile})
    sizearray=($sizeA)
    size=${sizearray[0]}
    echo ${size} > ${fileSizeFile}

    ## 2.5 sort tagalign file
    #input: un-zipped tagAlign
    #output: un-zipped tagAlign file

    sortedTagalignFile="${peakCallDir}${cell}.filt.rmdup.sort.tagAlign"
    echo "${cell} sorting tagalign"  >> $log_file
    sort -k1,1 -k2,2n ${tagalignFile} > ${sortedTagalignFile}

    ## 3. make bdg file
    #input: scaleFactor
    #       sorted un-zipped tagAlign file
    #       genome sizes file
    #output: bedGraph file

    SCALE_FACTOR=`awk -v nreads=${size} 'BEGIN {printf("%.20f", 1e6/nreads)}'`
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "${cell} total reads = ${size}      $dt"  >> $log_file
    echo "${cell} scale factor = ${SCALE_FACTOR}      $dt"  >> $log_file

    bdgFile="${peakCallDir}${cell}.scale_1e6.bdg"

    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "${cell} calling bedtools genomecov     $dt"  >> $log_file
    bedtools genomecov -i ${sortedTagalignFile} -bg -scale ${SCALE_FACTOR} -g ${GSIZE_FILE} > ${bdgFile}

    ## 4. sort bdg file
    #input: bedGraph file
    #output: sorted bedGraph file

    bdgSortedFile="${peakCallDir}${cell}.scale_1e6.sorted.bdg"
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "$cell sorting bdg file    $dt"  >> $log_file
    sort -k1,1 -k2,2n ${bdgFile} > ${bdgSortedFile}

    ## 5. make bigwig
    #input: sorted bedgraph file
    #output: bigWig file

    bwFile="${peakCallDir}${cell}.scale_1e6.bw"
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "${cell} making bigwig   $dt"  >> $log_file
    bedGraphToBigWig ${bdgSortedFile} ${GSIZE_FILE} ${bwFile}

    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "done making bw for ${cell}      $dt"  >> $log_file

    ## 6 make annotated bed file
    #input: narrowPeaks file created in step 1
    #output: allPeaksAnnoBed file

    narrowPeakFile="${peakCallDir}${cell}_peaks.narrowPeak"
    allPeaksAnnoBed="${peakCallDir}${cell}.all_peaks.anno.bed"
    awk -v  cellType=${cell} -v OFS='\t' ' { print $1, $2, $3, cellType } ' ${narrowPeakFile} > ${allPeaksAnnoBed}

    ## 7. sort annotated bed file
    #input: allPeaksAnnoBed file
    #output: sorted allPeaksAnnoBed file

    allPeaksAnnoSortBed="${peakCallDir}${cell}.all_peaks.anno.sort.bed"
    sort -k1,1 -k2,2n ${allPeaksAnnoBed} > ${allPeaksAnnoSortBed}

    ## 8. bedtools merge
    #input: sorted allPeaksAnnoBed file
    #output: mergedPeak file

    mergedPeakFile="${peakCallDir}${cell}.merged_peaks.anno.bed"
    bedtools merge -i  ${allPeaksAnnoSortBed} > ${mergedPeakFile}

    ## 9. bedtools intersect
    #input: mergedPeak file
    #       allPeaksAnnoSortBed file
    #output: intersect1 File

    intersect1File="${peakCallDir}${cell}.all_peaks.anno.INTERSECT1.bed"
    bedtools intersect -a ${mergedPeakFile} -b ${allPeaksAnnoSortBed} -wa -wb > ${intersect1File}

    ## 10. bedtools merge again
    #input: intersect 1 file
    #output: merged2PeakFile
    #NOTE: .merged_peaks.anno.bed is basically overwritten

    merged2PeakFile="${peakCallDir}${cell}.merged_peaks.anno.bed"
    bedtools  merge -i ${intersect1File} -c 7 -o distinct >  ${merged2PeakFile}
    dt=$(date '+%d/%m/%Y %H:%M:%S');
    echo "done merging peaks for ${cell}      $dt" >> $log_file

    ## 11. Remove all unnecessary/repetitive intermediates
    rm $allPeaksAnnoBed
    rm $allPeaksAnnoSortBed
    rm $intersect1File
    rm $bdgFile
    rm $tagalignFilegz

    ) &
done

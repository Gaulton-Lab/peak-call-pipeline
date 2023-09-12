# Example of Peak Call Pipeline

1. Make TagAligns for each sample 
    
    ```r
    ### run in parallel on all samples
    #!/bin/bash
    tagAlign_Pyscript_fp=/nfs/lab/projects/multiomic_islet/code/peak_calling/make_tagAligns.py
    samples=(R207 R217 R218 R221 R223 R226 R228 R234 R237 R238 R246 R247 R275 R277 R284 R290 R292 R316 R317 R319 R325 R326 R327 R332 R343 R349 R353 R354 R362 R363 R364)
    echo "${samples[@]}"
    N=8
    for sample in ${samples[@]}; do
    	((i=i%N)); ((i++==0)) && wait
    	(
    		sample_dir=/nfs/lab/projects/multiomic_islet/data/multiomics/cellranger/deep-shallow/${sample}/outs #path to directory with bam file (cellranger outputs)
    		cd $sample_dir
    		outs_dir=/nfs/lab/projects/multiomic_islet/outputs/multiome/indv_sample_processing/${sample}/
    		/usr/bin/python3 $tagAlign_Pyscript_fp -b "atac_possorted_bam.bam" -k ${outs_dir}filtered_barcodes.txt -o $outs_dir >> ${outs_dir}/log_make_tagAlign.txt
    	) &
    done
    
    ```
    
2. Merge the TagAligns into one file
    1. Modify the `merge_tagAligns.R` script in a few places: samples list, general path to `atac_keep_reads.tagAlign.gz`, and output path it writes to
3. Make barcode files for each celltype/leiden (needs to be a \n delimited .txt file)
    1. Doing this in a notebook: `/nfs/lab/hmummey/multiomic_islet/notebooks/221120_FINAL_combined_clustering_characterization_v2`. Basically uses the metadata of a Seurat object to write files of all barcodes per cell type. 
    2. Also make a file with the filepaths to each barcodes files on separate lines
4. Split merged tagAlign file by cell type (normally this only takes a few hours, but this time I started it at 1:21pm, then checked back at 4:46pm and it still wasn’t done running… eventually log said it was done at 7:31pm, this was a day with lots of high CPU processes on the server so was likely slower as a result.)
    1. Script: `/nfs/lab/katha/multiomics/islet-multiome-processing/call_peaks/splitTagAlign_parallelized.sh`
    2. Inputs:
        1. -c (file with all cell types in order)
        2. -t (merged tagAlign from step 2)
        3. -b (paths to all cell type specific barcodes files, same order as -c)
        4. -o outputs
    3. Command: 

    ```bash
    screen -r peak_calling
    
    cd /nfs/lab/katha/multiomics/islet-multiome-processing/call_peaks
    cells_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/celltypes.txt
    merged_tagAlign_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/31merged.atac_keep_reads.tagAlign2.gz
    ct_barcodes_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/barcodes_byCT.txt
    
    bash splitTagAlign_parallelized.sh -c $cells_fp -t $merged_tagAlign_fp -b $ct_barcodes_fp -o /nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/
    ```

5. Subsample cell types with more than 100,000,000 reads in the split tagAlign file
    1. Check the number of reads with `wc -l *tagAlign > tagAlign_counts.txt`
        1. Results
        
        ```bash
           672336397 split.acinar.tagAlign
          1779955636 split.alpha.tagAlign
          3276605905 split.beta.tagAlign
           352962627 split.delta.tagAlign
           140260923 split.ductal.tagAlign
             6573896 split.endothelial.tagAlign
           216682544 split.gamma.tagAlign
            15140738 split.immune.tagAlign
             1581251 split.schwann.tagAlign
            21924584 split.stellate.tagAlign
        ```
        
        1. Celltypes with > 100 mil: acinar, alpha, beta, delta, ductal, gamma, 
        2. Celltypes with < 100 mil: endothelial, immune, schwann, stellate
    2. Subsampling code
    
    ```bash
    screen -r peak_calling
    cd /nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2
    
    for ct in acinar alpha beta delta ductal gamma 
    do
    	orig_tagAlign="split."$ct".tagAlign"
    	new_tagAlign="split.subsample."$ct".tagAlign"
    	subsample --reservoir -n 100000000 $orig_tagAlign > $new_tagAlign
    	rm $orig_tagAlign
    done
    
    #fixing names
    for ct in endothelial immune schwann stellate
    do
    	orig_tagAlign="split."$ct".tagAlign"
    	new_tagAlign="split.subsample."$ct".tagAlign"
    	mv $orig_tagAlign $new_tagAlign
    done
    ```
    
6. Gzip the split tagAligns (as a separate file)

```bash
for ct in acinar alpha beta delta ductal endothelial gamma immune schwann stellate
do
	new_tagAlign="split.subsample."$ct".tagAlign"
	gzip -cvf $new_tagAlign > $new_tagAlign".gz"
done
```

7. Activate the conda environment (clone the environ from katha’s dir)
    1. `conda create --name call_peaks --clone /home/kakorgao/.conda/envs/kat_py_37/` 
    2. In screen: `conda activate /home/hmummey/.conda/envs/call_peaks`
    3. Alternative option: create the environment from `call_peaks.yml` (basically this environment just has `macs2 v2.2.7.1` and some other packages that I’m not sure are necessary)
8. Run the call peaks script!
    1. Inputs:
        1. -c (file with all cell types in order)
        2. -t (un-gzipped split tagAligns from step 5, same order as -c)
        3. -b (paths to all cell type specific barcodes files, same order as -c)
        4. -o output dir
        5. -g genome sizes file
    2. Command
    
    ```bash
    screen -r peak_calling
    cd /nfs/lab/katha/multiomics/islet-multiome-processing/call_peaks
    cells_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/celltypes.txt
    split_tagAlign_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/split.subsample.tagAligns_byCT.txt
    ct_barcodes_fp=/nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/barcodes_byCT.txt
    GSIZE_FILE=/nfs/lab/elisha/nPOD_output/scripts/references/hg38.chrom.sizes
    
    bash call_peaks_parallel.sh -c $cells_fp -t $split_tagAlign_fp -b $ct_barcodes_fp -o /nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/ -g $GSIZE_FILE    
    ```
    
9. Merge the peak files 

```bash
cd /nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/
echo `date`
awk -v OFS='\t' ' { print $1, $2, $3} ' /nfs/lab/projects/multiomic_islet/outputs/multiome/call_peaks/recluster_final_majorCTs_v2/*.merged* > allPeaksAnno.bed
sort -k1,1 -k2,2n allPeaksAnno.bed > allPeaksAnno_sorted.bed
bedtools merge -i  allPeaksAnno_sorted.bed > mergedPeak.txt
echo `date`
```

10. Overlap all individual celltype peak sets with the merged peak set (useful for downstream analyses so just doing this now)

```bash
for ct in acinar alpha beta delta ductal endothelial gamma immune schwann stellate
do 
	bedtools intersect -wa -a mergedPeak.txt -b ${ct}'.merged_peaks.anno.bed' > ${ct}'.merged_peaks.anno.mergedOverlap.bed'
done
```

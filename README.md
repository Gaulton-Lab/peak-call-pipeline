# Call Peaks Pipeline
This repository contains scripts and sample files designed to call peaks on multiomics seurat clustering outputs. Scripts and README.md were initially designed by Katha Korgaonkar, and recently adapted by Hannah Mummey and Rebecca Melton. See [this file](https://github.com/Gaulton-Lab/peak-call-pipeline/blob/main/example_peak_call_pipeline.md) for an example implementation of all steps of the pipeline with commands. Note: this pipeline supports calling peaks for any genome build, you just need to update the chrom.sizes file that is input into the pipeline script.
## Step 1 - Make TagAligns for each sample
- [Make tagAligns script](https://github.com/Gaulton-Lab/peak-call-pipeline/blob/main/scripts/make_tagAligns.py): `make_tagAligns.py`
    - Inputs: 
    - Outputs: .tagAlign.gz for each sample

## Step 2 - Merge Sample TagAligns into one file
- [Merge tagAligns script](https://github.com/Gaulton-Lab/peak-call-pipeline/blob/main/scripts/merge_tagAligns.R): `merge_tagAligns.R`
    - Inputs: list of samples and the location of their tagAligns
    - Outputs: merged .tagAlign.gz file

- **NOTE:** this script is hard coded - so be sure to adapt it to your needs

## Step 3 - Making barcode files for each celltype / leiden
- Pull out the barcodes for each cell type / leiden value into a \n delimted .txt file. You should have the same amount of .txt files as you do cell types/ leiden values.

## Step 4 - Split merged tagAlign file by cell type/leiden
- [Split tagAlign script](https://github.com/Gaulton-Lab/peak-call-pipeline/blob/main/scripts/splitTagAlign_parallelized.sh): `splitTagAlign_parallelized.sh`
    - Inputs: -c `cells.txt`, -t `tagAlign.tsv.gz`,  -b `barcodes.txt`, -o `/path/to/output/directory/`
        - `cells.txt` is a \n delimited .txt file with the names of all cell types you are calling peaks on 
		- `tagAlign.tsv.gz` is the output from step 2 
		- `barcodes.txt` is a \n delimited .txt file with the paths to the celltype/leiden barcode .txt files
	- **NOTE:** for the -b .txt file, the paths should be in the same order as the cell/leiden in the -c .txt file
    - Outputs: split .tagAlign.gz file for each celltype / leiden
- **SAMPLE COMMAND**: `bash splitTagAlign_parallelized.sh -c cells.txt -t /path/to/atac_reads.tagAlign.gz -b barcodes.txt -o /path/to/outputs/`

## Step 5 -  Subsample cell types/ leiden with more than 100 Million reads in split tagAligns
- Check the number of reads by using this command: `wc -l * > tagAlign_counts.txt`
    - **NOTE:** for the samples that have less than 100million -> rename files with suffix .subsample.tagAlign 
- `subsample --reservoir -n 100000000 path/to/.splitTagAlign > path/to/.subsample.tagAlign`

## Step 6 -  GZIP splitTagAligns
- **NOTE:** the call peaks script needs both gzip and unzipped versions of the tagAlign file
- `gzip -cvf path/to/.subsampletagAlign > path/to/.subsampletagAlign.gz`

## Step 7 - Activating conda Environment
- To get all the packages needed to run the call peaks script, you must run it in the correct conda environment
    - Create a new conda env from .yml file using cmd: `conda env create -f call_peaks_environment.yml`
    - Once created, run this command: `conda activate call_peaks`
## Step 8 - Running the Call Peaks Script
- [Call_peaks script](https://github.com/Gaulton-Lab/peak-call-pipeline/blob/main/scripts/call_peaks_parallel_v2.sh): `call_peaks_parallel_v2.sh`
    - Inputs: -c `cells.txt`, -t `tagAligns.txt`,  -b `barcodes.txt`, -o `/path/to/output/directory/`, -g `genome.chrom.sizes`
        - `cells.txt` is a \n delimited .txt file with the names of all cell types you are calling peaks on 
		- `tagAlign.txt` is a \n delimited .txt file with the paths to un-gzipped split tagAligns 
		- `barcodes.txt` contains the \n delimited paths to celltype/leiden barcode .txt files
		- `genome.chrom.sizes` is a \t delimited file with the length (in base pairs) of each chromosome in your genome build (including all random contigs)
	- **NOTE:** for the -t and -b .txt files, the paths should be in the same order as the cell/leiden in the -c .txt file
	- **NOTE:** gzipped versions of the splitTagAligns must be in the same directory as the unzipped ones
 - **SAMPLE COMMAND**: `bash call_peaks_parallel_v2.sh -c cell.txt -t /nfs/lab/katha/multiomics/scripts/cleaned_again/tag_list.txt -b barcodes.txt -o /path/to/outputs/ -g genome.chrom.sizes`
## Step 9 - Merge Peak Files
- [Merge peaks script](https://github.com/Gaulton-Lab/peak-call-pipeline/blob/main/scripts/mergePeaks.sh): `mergePeaks.sh`
    - **NOTE:** this script is hard coded - so be sure to adapt it to your needs

## Step 10 - Overlap cell type peak sets with merged peak set

Use bedtools intersect `-wa -a mergedPeak.txt -b celltype peaks` to get files of merged peaks found in each cell type. 

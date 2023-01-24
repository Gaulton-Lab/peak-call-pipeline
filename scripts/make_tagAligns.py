#!/usr/bin/env python3

import os
import sys
import gzip
import argparse
import subprocess
import pysam
import numpy as np
import pandas as pd
import scipy.sparse
from multiprocessing import Pool
from datetime import datetime


def main(args):
    print('Reading in the keep barcodes:',datetime.now())
    #read in the keep barcodes and store them as a set
    keep = open(args.keep).read().splitlines()
    keep_set = set(keep)

    #define necessary files
    tagalign_file = args.outdir + 'atac_keep_reads.tagAlign.gz'

    print('Reading through the bamfile:',datetime.now())
    #read through the bamfile and add lines to tagalign file if in keep_set barcodes
    bamfile = pysam.AlignmentFile(args.bam, 'rb')
    with gzip.open(tagalign_file, 'wt') as tagout:
        genome_size = {item['SN']:item['LN'] for item in bamfile.header['SQ']}
        for read in bamfile:
            if not read.is_proper_pair or read.is_duplicate or read.is_secondary or read.is_qcfail or read.is_supplementary:
                continue
            if read.mapping_quality < 30:
                continue
            if not read.has_tag('CB'):
                continue
            #get barcode and see if it's in the keep set, if not don't proceed with this read
            barcode = read.get_tag('CB')
            if not barcode in keep_set:
                continue
            read_chr = read.reference_name
            if not read_chr in ['chr{}'.format(i) for i in range(1,23)] + ['chrX','chrY','chrM']:
                continue
            read_start = max(1, read.reference_end - args.shift - args.extsize - 5 if read.is_reverse else read.reference_start + args.shift + 4)
            read_end = min(genome_size[read_chr], read.reference_end - args.shift - 5 if read.is_reverse else read.reference_start + args.shift + args.extsize + 4)
            read_qual = read.mapping_quality
            if read.is_reverse:
                read_orient = '-'
            else:
                read_orient = '+'
            print(read_chr, read_start, read_end, barcode, read_qual, read_orient, sep='\t', file=tagout)
        bamfile.close()

    #in the ATAC lfm script we sort the bed file we make... not sure if that's necessary with a tagalign
    #I don't think it is (not sorted in Josh's pipeline script), and we've already ignored any reads from
    #barcodes not in the keep_set so I think we're done here!
    print('Done:',datetime.now())


def process_args():
    parser = argparse.ArgumentParser(description='Filter fragments file')
    io_group = parser.add_argument_group('I/O arguments')
    io_group.add_argument('-b', '--bam', required=True, type=str, default='atac_possorted_bam.bam', help='Path to 10X output ATAC bam file')
    io_group.add_argument('-k', '--keep', required=True, type=str, help='List of barcodes to keep')
    io_group.add_argument('-o', '--outdir', required=True, type=str, help='Desired output directory')

    tagalign_group = parser.add_argument_group('tagAlign generation arguments')
    tagalign_group.add_argument('-s', '--shift', required=False, type=int, default=-100, help='Read shift length')
    tagalign_group.add_argument('-e', '--extsize', required=False, type=int, default=200, help='Read extension size')
    return parser.parse_args()


if __name__ == '__main__':
        #logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s'>
        args = process_args()
        main(args)

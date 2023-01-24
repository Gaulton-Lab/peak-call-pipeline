samples <- c('R207','R217','R218','R221','R223',
             'R226','R228','R234','R237','R238',
             'R246', 'R247', 'R275', 'R277', 'R284', 
             'R290', 'R292', 'R316', 'R317', 'R319',
             'R325','R326','R327','R332','R343',
             'R349','R353','R354','R362','R363','R364')
system(sprintf('for SAMPLE in %s; do zcat /nfs/lab/projects/multiomic_islet/outputs/multiome/indv_sample_processing/$SAMPLE/atac_keep_reads.tagAlign.gz | awk -v SAMPLE=$SAMPLE \'BEGIN{FS=OFS="\t"} {print $1,$2,$3,SAMPLE"_"$4,$5,$6}\'; done | sort -k1,1 -k2,2n -S 64G | bgzip -c -@ 16 > /nfs/lab/hmummey/multiomic_islet/intermediates/31merged.atac_keep_reads.tagAlign2.gz', paste0(samples, collapse=' '), samples[[1]]))

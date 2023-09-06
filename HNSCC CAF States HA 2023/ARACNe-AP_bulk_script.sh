#!/bin/sh

java -Xmx5G -jar aracne.jar -e 230710_caf_bulk/combined_cafs_counts_adj_filt.txt -o 230710_caf_bulk/output --tfs lambert2018.txt --pvalue 1E-8 --seed 123 --calculateThreshold

for i in `seq 1 1 100`
do
    java -Xmx5G -jar aracne.jar -e 230710_caf_bulk/combined_cafs_counts_adj_filt.txt -o 230710_caf_bulk/output --tfs lambert2018.txt --pvalue 1E-8 --seed $i
done

java -Xmx5G -jar aracne.jar -o 230710_caf_bulk/output --consolidate


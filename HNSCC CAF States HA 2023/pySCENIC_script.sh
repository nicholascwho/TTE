#!/bin/bash

cd /hpctmp/a0155826x/pyscenic

source /app1/ebenv && module load pySCENIC/0.12.1-foss-2022a

pyscenic grn --num_workers 24 --output output/adj.tsv --method grnboost2 pri_caf_counts.loom hs_hgnc_curated_tfs.txt
pyscenic ctx output/adj.tsv hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather --annotations_fname motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname pri_caf_counts.loom --mode "dask_multiprocessing" --output output/reg.csv --num_workers 24 --mask_dropouts
pyscenic aucell pri_caf_counts.loom output/reg.csv --output output/pri_caf_aucell.txt --num_workers 24

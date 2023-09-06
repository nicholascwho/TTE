#!/bin/sh

java -Xmx5G -jar aracne.jar -e 230608_caf/pri_caf_log1p_counts_sct.txt -o 230608_caf/output --tfs lambert2018.txt --pvalue 1E-8 --seed 123 --calculateThreshold

for i in `seq 1 10 100`
do
    java -Xmx5G -jar aracne.jar -e 230608_caf/pri_caf_log1p_counts_sct.txt -o 230608_caf/output --tfs lambert2018.txt --pvalue 1E-8 --seed $i
done

java -Xmx5G -jar aracne.jar -o 230608_caf/output --consolidate


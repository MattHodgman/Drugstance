# append individual output files into complete tsv
out=$1 # dir that output files are in
ls ${out}/* | xargs -I {} cat {} >> all_drug_distances.tsv # merge tsvs
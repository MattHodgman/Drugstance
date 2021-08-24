# append individual output files into complete tsv
out=$1 # dir that output files are in
file=$2 # name of output file
ls ${out}/* | xargs -I {} cat {} >> $file # merge tsvs
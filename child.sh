#!/bin/bash
#SBATCH -c 20                              # Request cores
#SBATCH -t 0-01:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem-per-cpu=1G                   # Memory total in GiB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID (%j)
                                           # You can change the filenames given with -o and -e to any filenames you'd like

# args
drugs_list=$1 # the list of drugs to be used in this child process
out=$2 # dir to write output file to

# get id
id=$(echo $drugs_list | grep -oP '(?<=_)\w+')

# source bashrc to activate conda in shell
source ~/.bashrc

# start env with python 3.8
conda activate py38env

# Run computeDistances.py using multiple processes
python3 computeSemanticDistances.py -n 20 -d $drugs_list -a data/drugs -g data/chembl.gpkl -p data/drug_node_dict.pkl -i $id -o $out

# close conda env
conda deactivate
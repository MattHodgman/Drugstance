# args
mesh=$1 # MeSH data in ASCII format
indications=$2 # TSV with the columns 'pref_name' and 'mesh_heading'

# prep data
mkdir data
python3 mapMesh.py -i $mesh -o data # map MeSH headings to their tree numbers and visa-versa
python3 makeGraph.py -i $indications -h mesh_headings.pkl -n mesh_numbers.pkl -o data # create a DAG of the MeSH tree
mv $indications data/

# compare all drugs
bash control.sh 10 $indications # compute both the semantic distance and overlap between all drugs, using SLURM and multiprocessing across 20 cores

# append intermediary output files together into final results
bash append.sh semantic_distances
bash append.sh overlaps
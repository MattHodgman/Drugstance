mesh=$1
indications=$2
mkdir data
python3 mapMesh.py -i $mesh -o data
python3 makeGraph.py -i $indications -h mesh_headings.pkl -n mesh_numbers.pkl -o data
mv $indications data/
bash control.sh 10
bash append.sh output
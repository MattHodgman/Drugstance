# Drugstance
Calculate the [semantic distance](https://academic.oup.com/bioinformatics/article/29/13/i53/195366) and [overlap coefficient](https://en.wikipedia.org/wiki/Overlap_coefficient) between drugs based on their indications.

## Example usage (Slurm)
When running this pipeline on a large number of drugs, it is recommended to use SLURM to split computations across multiple jobs, each multiprocessing across 20 cores. To do this simple run the following command:
```
bash pipeline d2021.bin indications.tsv
```
Note: Files create and used in the pipeline will be written to a directory called `data/` and output files will be written to a directory called `output/`.

## Example usage (Docker)
This is a standalone version of the algorithm that runs within Docker. There is no multiprocessing nor job usage. Therefore it will run very slowly and require a lot of memory. We recommend running it with a smaller set of drugs and indications.
```
docker run --rm -v "$PWD":/data labsyspharm/drugstance:latest python3 /app/drugstance -i d2021.bin -m indications.tsv -o /data/
``` 

## Input Files
`d2021.bin` is all [MeSH data](https://www.nlm.nih.gov/databases/download/mesh.html) downloaded in ASCII format. `indications.tsv` is a TSV file from [ChEMBL](https://www.ebi.ac.uk/chembl/) that contains in the column `pref_name` the name of the drug and in the column `mesh_heading` a valid MeSH heading that is an indication of that drug. For example:

pref_name | mesh_heading
--------- | ------------
TOFACITINIB | Immune System Diseases
TOFACITINIB | Arthritis, Rheumatoid
ASPIRIN | Pain
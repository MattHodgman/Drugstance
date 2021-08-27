import csv
import argparse
import math
import pandas as pd
from multiprocessing import Pool


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Compute the overlap between drugs in ChEMBL using the MeSH headings of their indications.')
    parser.add_argument('-n', '--num-cpus', help="The number of cpus to use. Default is 4.", default=4, type=int, required=True)
    parser.add_argument('-c', '--chembl', help='A TSV file containing ChEMBL drug indication information. Must contain the drug name under the column \'pref_name\' and the indication name under \'mesh_heading\'.', type=str, required=True)
    parser.add_argument('-d', '--drugs', help='A file containing a list of drugs for computation', type=str, required=True)
    parser.add_argument('-a', '--all-drugs', help='A file containing a list of all drugs in ChEMBL.', type=str, required=True)
    parser.add_argument('-i', '--id', help='A unique id to use for this scripts output file.', type=str, required=True)
    parser.add_argument('-o', '--output', help='Output directory.', default='.', type=str, required=True)
    args = parser.parse_args()
    return args


'''
Get a drug's indications
'''
def getIndications(drug):
    indications = set(chembl[chembl['pref_name'] == drug]['mesh_heading'].tolist())
    return indications


'''
Calculate the overlap coefficient for two drugs.
'''
def computeOverlap(drug1, drug2):
    # get the indications for each drug
    indications1 = getIndications(drug1)
    indications2 = getIndications(drug2)

    intersection = len(list(indications1 & indications2)) # get the size of the intersection of both sets
    min_size = min(len(indications1),len(indications2)) # get the size of the smaller set
    overlap = intersection / min_size # compute overlap 

    return overlap


'''
Calculate the overlap between every pairwise combination of drugs, no repeats.
'''
def run_comparisons(drugs, all_drugs):

    results = []
    for drug1 in drugs:
        overlaps = [drug1]
        for drug2 in all_drugs:
            overlap = computeOverlap(drug1, drug2) # compute overlap
            overlaps.append(overlap)
        results.append(overlaps)

    return results
    

'''
Initializer for multiprocessing to generate global variables to use in each proces.
'''
def initializer():
    global all_drugs
    global chembl

    # load in all drugs
    f = open(args.all_drugs, 'r')
    all_drugs = [line.rstrip() for line in f]

    # load in ChEMBL
    chembl = pd.read_csv(args.chembl, sep='\t')


'''
Main function for each process. Computes all comparisons.
'''
def main(drugs):
    results = run_comparisons(drugs, all_drugs)
    return results


'''
Write results to table.
'''
def write_results_new(overlaps):
    with open(f'{args.output}/drug_overlaps_{args.id}.tsv', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        if args.id == 'aa':
            all_drugs.insert(0,'Drug')
            writer.writerow(all_drugs)
        writer.writerows(overlaps)

    
'''
Main.
'''
if __name__ == "__main__":

    global args # hopefully this allows the subprocesses to access them

    # get arguments
    args = parseArgs() # parse arguments
    n = args.num_cpus # number of processes/cpus to use

    print('loading data...')

    # load in small drug list
    f = open(args.drugs, 'r')
    drugs = [line.rstrip() for line in f]

    # load in all drugs
    f = open(args.all_drugs, 'r')
    all_drugs = [line.rstrip() for line in f]

    # divide drugs for processes
    num_drugs = math.ceil(len(drugs) / n)
    lists = [drugs[i:i + num_drugs] for i in range(0, len(drugs), num_drugs)]
  
    print('running processes...')

    # compute overlaps across n threads
    with Pool(n, initializer, ()) as p:
        overlaps = p.map(main, lists)

    # merge results
    overlaps = [j for i in overlaps for j in i]

    print('writing results...')
    write_results_new(overlaps)
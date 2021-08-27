import argparse
import pandas as pd


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Average kernels for an improved metric of drug distance. Input matrices should be similarity scores between 0 and 1, where 1 indicates equivalence.')
    parser.add_argument('-i', '--input', help='A list of similarity matrices (TSV files) with the different kernels to be averaged.', nargs='*', action='store', dest='input', type=str, required=True)
    parser.add_argument('-o', '--output', help='Output directory.', default='.', type=str, required=False)
    args = parser.parse_args()
    return args


'''
Main.
'''
if __name__ == "__main__":
    args = parseArgs() # parse arguments
    DRUG = 'Drug' # column header to use as index
    matrices = []

    # read, sort, and convert to distances, all input matrices
    for matrix in args.input:
        X = pd.read_csv(matrix, delimiter='\t', index_col=DRUG) # read matrix into pandas
        
        X = X.sort_index() # sort rows
        X = X.reindex(sorted(X.columns), axis = 1) # sort columns
        X = 1 - X # convert to distance
        matrices.append(X)

    # add all the kernels together
    sum = matrices[0] # init df
    for X in matrices:
        sum += X
    
    # get the average
    sum = sum / len(matrices)

    # write average distance kernel to file
    sum.to_csv(f'{args.output}/drug_average_kernel.tsv', sep='\t')
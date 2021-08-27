import argparse
import pandas as pd
from sklearn.metrics.pairwise import rbf_kernel


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Compute the Radial Basis Function Kernel for a distance matrix')
    parser.add_argument('-i', '--input', help='A distance matrix in TSV format.', type=str, required=True)
    parser.add_argument('-s', '--sigma', help='The value of sigma to use in the RBF kernel.', default=10, type=int, required=False)
    parser.add_argument('-o', '--output', help='Output directory.', default='.', type=str, required=False)
    args = parser.parse_args()
    return args


'''
Main.
'''
if __name__ == "__main__":
    args = parseArgs() # parse arguments
    DRUG = 'Drug' # column header to use as index

    X = pd.read_csv(args.input, delimiter='\t', index_col=DRUG) # read matrix into df
    K = rbf_kernel(X, gamma = args.sigma) # do the RBG kernel magic!
    Y = pd.DataFrame(K, columns=X.columns).set_index(X.index)

    filename = args.input.split('/')[-1].split('.tsv')[0] # get input file name without extension
    Y.to_csv(f'{args.output}/{filename}_rbf_kernel.tsv', sep='\t') # write RBF kernel to file
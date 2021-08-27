'''
This script is adapted from https://code.tutsplus.com/tutorials/working-with-mesh-files-in-python-linking-terms-and-numbers--cms-28587
it takes as input the file from MeSH https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/asciimesh/d2021.bin but with the current year, this is the MeSH database in ascii format.
'''

import re
import argparse
import pickle


'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Map MeSH headings to their tree numbers and visa-versa. Write as dictionaries to pickle files.')
    parser.add_argument('-i', '--input', help='The MeSH database in ascii format', type=str, required=True)
    parser.add_argument('-o', '--output', help='Path to the output directory', default='.', type=str, required=False)
    args = parser.parse_args()
    return args


'''
Main.
'''
if __name__ == '__main__':

    # parse arguments
    args = parseArgs()
 
    terms = {}
    numbers = {}
    
    # read input
    meshFile = args.input
    with open(meshFile, mode='rb') as file:
        mesh = file.readlines()

    # parse input
    for line in mesh:
        meshTerm = re.search(b'MH = (.+)$', line)
        if meshTerm:
            term = meshTerm.group(1)
            term = term.decode()
        meshNumber = re.search(b'MN = (.+)$', line)
        if meshNumber:
            number = meshNumber.group(1)
            numbers[number.decode()] = term
            if term in terms:
                terms[term].append(number.decode('utf-8'))
            else:
                terms[term] = list()
                terms[term].append(number.decode('utf-8'))

    # write dicts to pickle files
    f = open(f'{args.output}/mesh_numbers.pkl', 'wb')
    pickle.dump(numbers, f)

    f = open(f'{args.output}/mesh_headings.pkl', 'wb')
    pickle.dump(terms, f)
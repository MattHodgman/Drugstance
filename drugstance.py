import re
import argparse
import networkx as nx
import math
import pandas as pd
import numpy as np
import csv

'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Compute the semantic distance and overlap between drugs in ChEMBL using the MeSH headings of their indications.')
    parser.add_argument('-i', '--input', help='A TSV file containing ChEMBL drug indication information. Must contain the drug name under the column \'pref_name\' and the valid MeSH heading (indication) under \'mesh_heading\'.', type=str, required=True)
    parser.add_argument('-m', '--mesh', help='A file all MeSH data in ASCII format', type=str, required=True)
    parser.add_argument('-o', '--output', help='Output directory.', default='.' type=str, required=False)
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
Calculate mis-information between two drug graphs.
'''
def computeMI(drug1, drug2):
    nodes = np.setdiff1d(list(drug_node_dict[drug1]), list(drug_node_dict[drug2]))

    mi = 0
    for n in nodes:
        mi += G.nodes[n]['ia']

    return mi


'''
Compute remaining uncertainty between two drug graphs.
'''
def computeRU(drug1, drug2):
    nodes = np.setdiff1d(list(drug_node_dict[drug2]), list(drug_node_dict[drug1]))

    ru = 0
    for n in nodes:
        ru += G.nodes[n]['ia']

    return ru


'''
Calculate te semantic distance (sd) between two drugs by summing the mis-information and remaining uncertainty values.
'''
def semanticDistance(drug1, drug2):
    sd = computeMI(drug1, drug2) + computeRU(drug1, drug2)
    return sd


'''
Calculate the semantic distance between every pairwise combination of drugs, no repeats.
'''
def runComparisons(drugs):
    all_sd = []
    all_o = []
    for drug1 in drugs:
        s_distances = [drug1]
        overlaps = [drug1]
        for drug2 in drugs:

            # compute semantic distance
            sd = semanticDistance(drug1, drug2)
            s_distances.append(sd)

            # compute overlap
            overlap = computeOverlap(drug1, drug2)
            overlaps.append(overlap)

        all_sd.append(s_distances)
        all_o.append(overlaps)

    return all_sd, all_o


'''
Write results to a TSV.
'''
def writeResults(results, f):
    with open(f'{args.output}/{f}', 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        drugs.insert(0,'Drug')
        writer.writerow(drugs)
        writer.writerows(results)


'''
Get the parent number (move UP the tree one level/generation).
'''
def up(n):
    sep = '.'
    n = n.split(sep)
    if len(n) > 1:
        n.pop()
        n = sep.join(n)
    else:
        return n[0][0]
    return n


'''
Get indications headings.
'''
def getHeadings(drug):
    indications = chembl[chembl['pref_name'] == drug.upper()]
    headings = sorted(list(set(indications.mesh_heading)))
    return headings


'''
Make a graph where the nodes are MeSH headings and the directed edges form the heirarchy. Each node has as an attribute the list of drugs that include that node.
'''
def makeGraph(drugs):

    node_drug_dict = {}

    # initialize DAG
    G = nx.DiGraph()

    for drug in drugs:
        drug_node_dict[drug] = set()
        for heading in getHeadings(drug):
            # add heading node
            if heading not in list(G.nodes):
                G.add_node(heading)
                node_drug_dict[heading] = set()
                
            drug_node_dict[drug].add(heading)
            node_drug_dict[heading].add(drug)

            # add all parents
            for n in mesh_headings[heading]:
                c_heading = heading
                for i in range(n.count('.') + 1):
                    p = up(n) # get parent number
                    p_heading = mesh_numbers[p] # get parent heading

                    # add parent heading node
                    if p_heading not in list(G.nodes):
                        G.add_node(p_heading)
                        node_drug_dict[p_heading] = set()
                    
                    drug_node_dict[drug].add(p_heading)
                    node_drug_dict[p_heading].add(drug)
                    
                    if not G.has_edge(p_heading, c_heading):
                        G.add_edge(p_heading, c_heading) # add directed edge from parent to child

                    n = p # move up path
                    c_heading = p_heading
    
    nx.set_node_attributes(G, node_drug_dict, 'drugs')

    return G


'''
Compute information accretion for each node in a graph.
'''
def computeIA(G):
    # annotate nodes with information accretion value from probability
    node_ia_dict = {}

    for node in G.nodes:
        n_drugs = len(G.nodes[node]['drugs']) # get the number of drugs with this heading 

        ## get the number of drugs with all parent headings

        parents = G.predecessors(node) # get parent nodes

        # get drugs for each parent
        drug_sets = []
        for parent in parents:
            drug_sets.append(G.nodes[parent]['drugs'])

        n_p_drugs = len(set().union(*drug_sets)) # count the number of drugs that have all parent headings

        # compute probability
        if n_p_drugs == 0:
            prob = 1
        else:
            prob = n_drugs / n_p_drugs

        # compute information accretion
        ia = -math.log2(prob) # does this work for single values or does it have to be an array?

        node_ia_dict[node] = ia

    nx.set_node_attributes(G, node_ia_dict, 'ia')

    return G


'''
Create two dictionaries:
    mesh_headings: (key) MeSH heading as string (value) list of MeSH tree numbers as strings
    mesh_numbers: (key) MeSH tree number as string (value) MeSH heading as string
'''
def mapMeSH(meshFile):
    mesh_headings = {}
    mesh_numbers = {}
    
    # read input
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
            mesh_numbers[number.decode()] = term
            if term in mesh_headings:
                mesh_headings[term].append(number.decode('utf-8'))
            else:
                mesh_headings[term] = list()
                mesh_headings[term].append(number.decode('utf-8'))
    
    return mesh_headings, mesh_numbers


'''
Main.
'''
if __name__ == "__main__":
    args = parseArgs() # parse arguments

    print('Loading input data...')
    chembl = pd.read_csv(args.input, sep='\t') # load in ChEMBL
    drugs = list(set(chembl['pref_name'])) # get drug list
    drug_node_dict = {} # init dict

    print('Mapping MeSH headings and numbers...')
    mesh_headings, mesh_numbers = mapMeSH(args.mesh) # map MeSH headings to numbers and visa-versa

    print('Making graph...')
    G = makeGraph(drugs) # build graph

    print('Computing information accretion...')
    G = computeIA(G) # compute and add ia values

    print(f'Computing semantic distance and overlap between all {len(drugs)} drugs...')
    semantic_distances, overlaps = runComparisons(drugs)

    print('Writing results...')
    writeResults(semantic_distances, 'semantic_distances.tsv')
    writeResults(overlaps, 'overlaps.tsv')
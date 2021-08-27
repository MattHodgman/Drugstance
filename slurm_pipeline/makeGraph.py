import pickle
import networkx as nx
import math
import pandas as pd
import argparse

'''
Parse arguments. None are required.
'''
def parseArgs():
    parser = argparse.ArgumentParser(description='Construct a graph using the MeSH headings of drugs in ChEMBL. Write graph and a drug:node dict to pickle files.')
    parser.add_argument('-i', '--input', help='A TSV of drugs and their indications.', type=str, required=True)
    parser.add_argument('-o', '--output', help='Path to output directory.', default='.', type=str, required=False)
    parser.add_argument('-h', '--headings', help='A pickle file of a dictionary that has MeSH headings as keys and a list their numbers as values.', type=str, required=True)
    parser.add_argument('-n', '--numbers', help='A pickle file of a dictionary that has MeSH numbers as keys and their headings as values.', type=str, required=True)
    args = parser.parse_args()
    return args


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
def get_headings(drug):
    indications = chembl[chembl['pref_name'] == drug.upper()]
    headings = sorted(list(set(indications.mesh_heading)))
    return headings


'''
Make a graph where the nodes are MeSH headings and the directed edges form the heirarchy. Each node has as an attribute the list of drugs that include that node.
'''
def make_graph(drugs):

    node_drug_dict = {}

    # initialize DAG
    G = nx.DiGraph()

    for drug in drugs:
        drug_node_dict[drug] = set()
        for heading in get_headings(drug):
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
def compute_ia(G):
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
Main.
'''
if __name__ == '__main__':

    print('loading data...')

    # parse arguments
    args = parseArgs()

    # load in MeSH pickle files
    f = open(args.headings, 'rb')
    mesh_headings = pickle.load(f)
    f = open(args.numbers, 'rb')
    mesh_numbers= pickle.load(f)

    # load in ChEMBL
    chembl = pd.read_csv(args.input, sep='\t')

    drugs = list(set(chembl['pref_name'])) # get drug list

    drug_node_dict = {}

    print('making graph...')
    G = make_graph(drugs) # build graph

    print('computing information accretion...')
    G = compute_ia(G) # compute and add ia values

    print('writing output to pickle files...')
    # save graph to pickle file
    nx.write_gpickle(G, f'{args.out}/chembl.gpkl')

    # save drug_node_dict to pickle
    dict_file = open(f'{args.out}/drug_node_dict.pkl', 'wb')
    pickle.dump(drug_node_dict, dict_file)
#Useful functions

import networkx as nx
import hypernetx as hnx

'''
Input:  hypernetx hypergraph
        h_id : hyperedge id

Output: set of hyperedge nodes
'''
def getHyperedgeNodeSet(hypergraph, h_id):
    output = set()
    for node in hypergraph.edges[h_id]:
        output.add(node)
    return output

'''
Input: two sets s1, s2
Output: jaccard similarity between s1 and s2
'''
def jaccard_similarity(s1, s2):
    intersection = s1.intersection(s2)
    union = s1.union(s2)
    return len(intersection)/len(union)

'''
Input: two sets s1, s2
Output: intersection between s1 and s2
'''
def common_neighbors_similarity(s1, s2):
    intersection = s1.intersection(s2)
    return len(intersection)


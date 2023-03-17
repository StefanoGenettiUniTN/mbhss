import networkx as nx
import hypernetx as hnx
import summary
import motif

def mbhss(component_adj, summary, k):
    '''
    Given an input hypergraph H(X,E) and integer k, find a summary hypergraph S
    for H with at most k supernodes X (|X| <= k)
    '''
    print("mbhss summarization algorithm start")

    #edges
    #{(87, 791), (491, 615), (39, 265, 720), ..., }
    edges = set()    
    for edge_id in summary.summaryHypergraph.edges:
        edge = []
        for supernode_id in summary.summaryHypergraph.edges[edge_id]:
            edge.append(supernode_id)

        if len(edge)<=4:    #chiedere a LOTITO se cambia qualcosa mettere 4 nel caso di ricerca motivi di ordine 3
            edges.add(tuple(edge))
    
    print(edges)

    output = {}
    
    print("====================motif 3 analysis ========================")
    output['motifs_3'] = motif.motifs_order_3(edges)
    print("============================================")
    print("")
    print("====================motif 4 analysis ========================")
    output['motifs_4'] = motif.motifs_order_4(edges)
    print("===========================================")
    print("")

    print("")
    print("motifs_3 counter")
    print(output['motifs_3'])
    print("")
    print("motifs_4 counter")
    for to_be_printed in output['motifs_4']:
        if to_be_printed[1] >= 1:
            print(to_be_printed)
import networkx as nx
import hypernetx as hnx
import summary
import motif_lib

def mbhss(component_adj, summary, k):
    '''
    Given an input hypergraph H(X,E) and integer k, find a summary hypergraph S
    for H with at most k supernodes X (|X| <= k)
    '''
    print("mbhss summarization algorithm start")

    #edges (we need to transform the set of edges according the requirements of the motif analysis algorithm)
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
    output['motifs_3'], h3_motif_set = motif_lib.motifs_order_3(edges)

    for m3 in h3_motif_set:
        print(m3)
    
    for to_be_printed in output['motifs_3']:
        if to_be_printed[1]>0:
            print(to_be_printed)

    print("============================================")
    print("")
    print("====================motif 4 analysis ========================")
    output['motifs_4'], h4_motif_set = motif_lib.motifs_order_4(edges)

    for m4 in h4_motif_set:
        print(m4)

    for to_be_printed in output['motifs_4']:
        if to_be_printed[1]>0:
            print(to_be_printed)
    print("===========================================")
    print("")

    #######################################################################################################

    #for the moment (march the 28th) we merge order three motifs
    #in the order they have been visisted by the algorithm
    for m3 in h3_motif_set:
        print(f"visiting motif: {m3.components}")

        #before the merge of the supernodes which populate the motif
        #we check wether these supernodes currently exist
        #indeed it is possible that the supernodes at hand have already been
        #merged in a previous iteration
        if check_motif(summary, m3):
            print("OK the motif can be merged in a new supernode")
            motif_lib.motif_condense(summary, m3)
            print("condense complete")
        else:
            print("unable to merge the motif")
        

'''
Check that the supernodes which populate the input_motif
are currrently in the input_summary.
Iff all the supernodes are in the input_summary return true.
'''
def check_motif(input_summary, input_motif):
    for c in input_motif.components:
        if input_summary.getSupernode(c)==None:
            return False
    return True

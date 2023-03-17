#this file contains the motif analysis algorithm. The code
#has been taken from the following repository made available
#by Quintino Francesco Lotito et al.

#https://github.com/FraLotito/higher-order-motifs.git

import math
import itertools
import numpy as np

#Motif order 3 analysis
def motifs_order_3(edges):
    N = 3
    full, visited = motifs_ho_full(edges, N)   #motifs with higher order interaction
    standard = motifs_standard(edges, N, visited)  #motifs composed only by pairwise interactions

    #print("print full")
    #print(full)
    #print("end print full")

    #print("print standard")
    #print(standard)
    #print("end print standard")

    res = []
    for i in range(len(full)):
        res.append((full[i][0], max(full[i][1], standard[i][1])))   #non capisco perché venga fatto il massimo. Infatti teoricamente full dovrebbe contemplare i motivi con high order motifs, mentre standard i motivi composti solo da low order interactions

    #print("print res")
    #print(res)
    #print("end print res")

    return res


#Motif order 4 analysis
def motifs_order_4(edges):
    N = 4
    full, visited = motifs_ho_full(edges, N, )
    not_full, visited = motifs_ho_not_full(edges, N, visited)
    standard = motifs_standard(edges, N, visited)

    res = []
    for i in range(len(full)):
        res.append((full[i][0], max([full[i][1], not_full[i][1], standard[i][1]])))

    return res

############################################################
############################################################

#count motifs that contains the higher order interaction
#of cardinality N
def motifs_ho_full(edges, N):
    mapping, labeling = generate_motifs(N)

    #print("print mapping")
    #print(mapping)
    #print("end mapping")

    #print mapping for N=3
    #{
    # ((1, 2, 3),): {((1, 2, 3),)},
    # ((1, 2), (1, 2, 3)): {((1, 2), (1, 2, 3)), ((1, 2, 3), (2, 3)), ((1, 2, 3), (1, 3))},
    # ((1, 2), (1, 3)): {((1, 2), (2, 3)), ((1, 2), (1, 3)), ((1, 3), (2, 3))},
    # ((1, 2), (1, 2, 3), (1, 3)): {((1, 2), (1, 2, 3), (2, 3)), ((1, 2), (1, 2, 3), (1, 3)), ((1, 2, 3), (1, 3), (2, 3))},
    # ((1, 2), (1, 3), (2, 3)): {((1, 2), (1, 3), (2, 3))},
    # ((1, 2), (1, 2, 3), (1, 3), (2, 3)): {((1, 2), (1, 2, 3), (1, 3), (2, 3))}
    # }
    #end mapping

    #print("print labeling")
    #print(labeling)
    #print("end labeling")

    #print labeling for N=3
    #{
    # ((1, 2, 3),): 0,
    # ((1, 2), (1, 2, 3)): 0,
    # ((1, 2, 3), (1, 3)): 0,
    # ((1, 2, 3), (2, 3)): 0,
    # ((1, 2), (1, 3)): 0,
    # ((1, 2), (2, 3)): 0,
    # ((1, 3), (2, 3)): 0,
    # ((1, 2), (1, 2, 3), (1, 3)): 0,
    # ((1, 2), (1, 2, 3), (2, 3)): 0,
    # ((1, 2, 3), (1, 3), (2, 3)): 0,
    # ((1, 2), (1, 3), (2, 3)): 0,
    # ((1, 2), (1, 2, 3), (1, 3), (2, 3)): 0
    # }
    #end labeling

    T = {}
    for e in edges:
        T[tuple(sorted(e))] = 1

    visited = {}

    def count_motif(nodes):

        #nodes is a list of exactly length N
        #so if we are looking for order three motifs, then
        #nodes have always length 3

        print("start count motif")

        nodes = tuple(sorted(tuple(nodes)))
        p_nodes = power_set(nodes)
        
        #nodes
        #(156, 201, 869)
        
        #p_nodes (all possible sets)
        #[[], [156], [201], [156, 201], [869], [156, 869], [201, 869], [156, 201, 869]]

        motif = []
        for edge in p_nodes:
            if len(edge) >= 2:
                edge = tuple(sorted(list(edge)))
                if edge in T:   #l'insieme di nodi edge esiste in T, ovvero esiste nell'ipergrafo in input
                    motif.append(edge)

        print("motif")
        print(motif)
        print("end print motif")

        #start count motif
            #motif
            #[(177, 489), (177, 637), (489, 637), (177, 489, 637)]
            #end print motif
        #...end count motif

        #start count motif
            #motif
            #[(185, 254), (185, 258), (254, 258), (185, 254, 258)]
            #end print motif
        #...end count motif

        #start count motif
            #motif
            #[(295, 674), (295, 954), (674, 954), (295, 674, 954)]
            #end print motif
        #...end count motif


        #Assign a unique ID to each node
        #Example:
        #nodes = (38, 219, 866)
        #m = {38: 1, 219: 2, 866: 3}
        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        #update counter
        labeled_motif = []
        for e in motif:
            new_e = []
            for node in e:
                new_e.append(m[node])
            new_e = tuple(sorted(new_e))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))

        if labeled_motif in labeling:
            labeling[labeled_motif] += 1

        print("...end count motif")

    for e in edges:
        if len(e) == N: #consider iff the size of the hyperedge is N
            #print(e)
            visited[e] = 1
            nodes = list(e)

            #edge: (38, 219, 866)
            #nodes: [38, 219, 866]

            count_motif(nodes)  #search motifs and update counter

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]
            
        out.append((motif, count))

    out = list(sorted(out))

    #counter of each higher order motif of order N found in the hypergraph
    #[(((1, 2), (1, 2, 3)), 58), (((1, 2), (1, 2, 3), (1, 3)), 231), (((1, 2), (1, 2, 3), (1, 3), (2, 3)), 1802), (((1, 2), (1, 3)), 0), (((1, 2), (1, 3), (2, 3)), 0), (((1, 2, 3),), 0)]

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    return out, visited

#used in motifs_order_4: counts order motifs of order 4
#which are composed of order 3 interactions
def motifs_ho_not_full(edges, N, visited):
    mapping, labeling = generate_motifs(N)

    T = {}
    graph = {}
    for e in edges:
        if len(e) >= N:
            continue

        T[tuple(sorted(e))] = 1

        for e_i in e:
            if e_i in graph:
                graph[e_i].append(e)
            else:
                graph[e_i] = [e]

    def count_motif(nodes):
        print("start count motif")
        nodes = tuple(sorted(tuple(nodes)))
        p_nodes = power_set(nodes)
        
        motif = []
        for edge in p_nodes:
            if len(edge) >= 2:
                edge = tuple(sorted(list(edge)))
                if edge in T:
                    motif.append(edge)

        print("motif")
        print(motif)
        print("end print motif")

        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e = []
            for node in e:
                new_e.append(m[node])
            new_e = tuple(sorted(new_e))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))

        if labeled_motif in labeling:
            labeling[labeled_motif] += 1
        
        print("...end count motif")

    for e in edges:
        if len(e) == N - 1: #in this function we consider only order N-1 interactions
            nodes = list(e)
            
            for n in nodes:
                for e_i in graph[n]:
                    tmp = list(nodes)
                    tmp.extend(e_i)
                    tmp = list(set(tmp))
                    if len(tmp) == N and not (tuple(sorted(tmp)) in visited):
                        visited[tuple(sorted(tmp))] = 1
                        count_motif(tmp)

    out = []

    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]
            
        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    return out, visited

#in visisted there are order N hyperedges which have been already visited
#in this function we consider only pairwise interactions
def motifs_standard(edges, N, visited):
    mapping, labeling = generate_motifs(N)

    graph = {}  #dictionary: keys are the nodes, the values are the adjacency list of the node
    T = {}

    z = set()   #z contains all the nodes of the hypergraph
    for e in edges:
        for n in e:
            z.add(n)

    for e in edges:
        #in this function only pairwise interactions are considered
        if len(e) == 2:
            T[tuple(sorted(e))] = 1 #annotate that the edge exists
            a, b = e
            if a in graph:
                graph[a].append(b)
            else:
                graph[a] = [b]

            if b in graph:
                graph[b].append(a)
            else:
                graph[b] = [a]

    def count_motif(nodes):
        print("inizia count_motif")
        nodes = tuple(sorted(tuple(nodes)))

        if nodes in visited:
            return

        p_nodes = power_set(nodes)
        
        motif = []
        for edge in p_nodes:
            edge = tuple(sorted(list(edge)))
            if edge in T:
                motif.append(edge)

        print("print motif")
        print(motif)
        print("end print motif")

        #inizia count_motif
            #print motif
            #[(252, 494), (252, 720), (494, 720)]
            #end print motif
        #finisce count_motif
        
        #inizia count_motif
            #print motif
            #[(252, 720), (366, 720)]
            #end print motif
        #finisce count_motif
        
        #inizia count_motif
            #print motif
            #[(252, 626), (252, 720)]
            #end print motif
        #finisce count_motif

        #herafter aggiorna contatore dei motivi
        m = {}
        idx = 1
        for i in nodes:
            m[i] = idx
            idx += 1

        labeled_motif = []
        for e in motif:
            new_e = []
            for node in e:
                new_e.append(m[node])
            new_e = tuple(sorted(new_e))
            labeled_motif.append(new_e)
        labeled_motif = tuple(sorted(labeled_motif))

        if labeled_motif in labeling:
            labeling[labeled_motif] += 1
        else:
            #sembrerebbe non accadere mai
            print("this is not in labeling")
            print(labeled_motif)
            print("-")

        print("finisce count_motif")

    #count_motif is called here
    def graph_extend(sub, ext, v, n_sub):

        if len(sub) == N:
            count_motif(sub)
            return

        while len(ext) > 0:
            w = ext.pop()
            tmp = set(ext)

            for u in graph[w]:
                if u not in sub and u not in n_sub and u > v:
                    tmp.add(u)

            new_sub = set(sub)
            new_sub.add(w)
            new_n_sub = set(n_sub).union(set(graph[w]))
            graph_extend(new_sub, tmp, v, new_n_sub)

    c = 0
    
    k = 0
    for v in graph.keys():
        v_ext = set()
        for u in graph[v]:
            if u > v:
                v_ext.add(u)
        k += 1

        graph_extend(set([v]), v_ext, v, set(graph[v]))
        c += 1

    out = []

    #update motif counter
    for motif in mapping.keys():
        count = 0
        for label in mapping[motif]:
            count += labeling[label]
            
        out.append((motif, count))

    out = list(sorted(out))

    D = {}
    for i in range(len(out)):
        D[i] = out[i][0]

    return out

#############################################################
#############################################################

def generate_motifs(N):
    n = N
    assert n >= 2

    h = [i for i in range(1, n + 1)]
    A = []

    for r in range(n, 1, -1):
        A.extend(list(itertools.combinations(h, r)))

    B = power_set(A)

    C = []
    for i in range(len(B)):
        if is_connected(B[i], N):
            C.append(B[i])

    isom_classes = {}

    for i in C:
        edges = sorted(i)
        relabeling_list = list(itertools.permutations([j for j in range(1, n + 1)]))
        found = False
        for relabeling in relabeling_list:
            relabeling_i = relabel(edges, relabeling)
            #print(relabeling_i)
            if tuple(relabeling_i) in isom_classes:
                found = True
                break
        if not found:
            isom_classes[tuple(edges)] = 1

    mapping = {}
    labeling = {}

    for k in isom_classes.keys():
        mapping[k] = set()
        relabeling_list = list(itertools.permutations([j for j in range(1, n + 1)]))
        for relabeling in relabeling_list:
            relabeling_i = relabel(k, relabeling)
            labeling[tuple(sorted(relabeling_i))] = 0
            mapping[k].add(tuple(sorted(relabeling_i)))
    
    return mapping, labeling


############################################################
############################################################

def power_set(A): 
    subsets = []
    N = len(A)

    for mask in range(1<<N):
        subset = []

        for n in range(N):
            if ((mask>>n)&1) == 1:
                subset.append(A[n])

        subsets.append(subset)

    return subsets

def is_connected(edges, N):
    nodes = set()
    for e in edges:
        for n in e:
            nodes.add(n)

    if len(nodes) != N:
        return False

    visited = {}
    for i in nodes:
        visited[i] = False
    graph = {}
    for i in nodes:
        graph[i] = []
    
    for edge in edges:
        for i in range(len(edge)):
            for j in range(len(edge)):
                if edge[i] != edge[j]:
                    graph[edge[i]].append(edge[j])
                    graph[edge[j]].append(edge[i])
    
    q = []
    nodes = list(nodes)
    q.append(nodes[0])
    while len(q) != 0:
        v = q.pop(len(q) - 1)
        if not visited[v]:
            visited[v] = True
            for i in graph[v]:
                q.append(i)
    conn = True
    for i in nodes:
        if not visited[i]:
            conn = False
            break
    return conn

def relabel(edges, relabeling):
    res = []
    for edge in edges:
        new_edge = []
        for v in edge:
            new_edge.append(relabeling[v - 1])
        res.append(tuple(sorted(new_edge)))
    return sorted(res)
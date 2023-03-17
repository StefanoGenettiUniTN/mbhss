import supernode
import matplotlib.pyplot as plt
import networkx as nx
import hypernetx as hnx
from itertools import combinations

class Summary:
    def __init__(self):
        """
        superList: dictionary{key:=supernode unique identifier ; value:=corresponding supernode object}
        autoincrementId: progressive identifier of the supernodes and hyperedges which populate the data structure. This counter is used to assign auto-increment unique integer identifiers. DO NOT use this attribute to count the number of supernodes which populate the data structure.
        components: dictionary{key:=id of the node of the original graph ; value:=supernode id where the node is in the summary}
        numComponents: number of components of the original graph represented by the summary
        """
        self.superList = {}
        self.autoincrementId = 0
        self.components = {}
        self.numComponents = 0
        self.summaryHypergraph = hnx.Hypergraph()

        #Remark:    supernodes and hyperepdges share the same autoincrementId because according to
        #           hypernetx specifications, it is not possible to have an hyperedge with the same
        #           id of a node

    def addComponent(self, key):
        """
        Add a new component to the summary.
        Consequently a new supernode is created to host the newcomer.
        """
        
        #add 1 to the number of components
        self.numComponents += 1

        #create a supernode with one component
        newSupernode = supernode.Supernode(self.autoincrementId)
        newSupernode.addComponent(key)
        self.superList[self.autoincrementId] = newSupernode

        #link components "key" to the corresponding supernode
        self.components[key] = self.autoincrementId

        #increment the number of supernodes
        self.autoincrementId += 1

    def getSupernode(self, key):
        """
        If supernode with key is in the summary then return
        the Supernode.
        """

        #use the get method to return the Supernode if it
        #exists otherwise it will return None
        return self.superList.get(key)

    def getComponentSupernode(self, componentId):
        """
        Find the supernode which containes the component with id componentId.
        If there it does not exist, return None.
        """
        if componentId in self.components:
            return self.superList.get(self.components[componentId])
        
        return None

    def getComponentSupernodeId(self, componentId):
        """
        Find the id of thesupernode which containes the component with id componentId.
        If there it does not exist, return None.
        """
        if componentId in self.components:
            return self.components[componentId]
        
        return None

    def getVertices(self):
        """
        Return all the supernodes in the summary
        """

        return self.superList.keys()

    def getComponents(self):
        """
        Return all the components represented by the summary
        """

        return self.components.keys()
    
    def insertHyperedge(self, elements):
        newHyperedge = hnx.Entity(self.autoincrementId, elements)
        self.autoincrementId += 1
        self.summaryHypergraph.add_edge(newHyperedge)

        #Update the neighbourhoods of each component of the new hyperedge
        possibleElementsCouples = list(combinations(elements, 2))
        for couple in possibleElementsCouples:
            c1 = couple[0]
            c2 = couple[1]
            
            #Get the corresponding supernodes
            s1 = self.getSupernode(c1)
            s2 = self.getSupernode(c2)

            s1.updateNeighbor(s2.getId(), 1)
            s2.updateNeighbor(s1.getId(), 1)

            s1.setIntersectionProfile(newHyperedge.uid, 1)
            s2.setIntersectionProfile(newHyperedge.uid, 1)

    def merge(self, s1, s2):
        """
        Merge supernode identified by id s1 with supernode identified by id s2.
        Returns the identifier of the new supernode.
        """
        #Check if supernode s1 != s2
        if s1==s2:
            print(f"Error function merge({s1},{s2}): cannot merge {s1} with itself")
            return

        #Get supernode s1
        supernodeS1 = self.getSupernode(s1)
        if supernodeS1==None:
            print(f"Error function merge({s1},{s2}): supernode {s1} does not exist")
            return

        #Get supernode s2
        supernodeS2 = self.getSupernode(s2)
        if supernodeS2==None:
            print(f"Error function merge({s1},{s2}): supernode {s2} does not exist")
            return
        
        #Create new supernode
        supernodeS3 = supernode.Supernode(self.autoincrementId)
        self.superList[self.autoincrementId] = supernodeS3
        self.autoincrementId += 1

        #Add to supernodeS3 all the components of supernodeS1
        for c in supernodeS1.components:
            supernodeS3.addComponent(c)
            self.components[c] = supernodeS3.id

        #Add to supernodeS3 all the components of supernodeS2
        for c in supernodeS2.components:
            supernodeS3.addComponent(c)
            self.components[c] = supernodeS3.id
        
        ##Update S3 partecipation (intersection profile) with summary hyperedges
        for s1_h_partecipation in supernodeS1.getAdjHyperedges():
            s1_h_intersection_profile = supernodeS1.getIntersectionProfile(s1_h_partecipation)
            supernodeS3.setIntersectionProfile(s1_h_partecipation, s1_h_intersection_profile)

            self.summaryHypergraph.add_node_to_edge(supernodeS3.getId(), s1_h_partecipation)

        for s2_h_partecipation in supernodeS2.getAdjHyperedges():
            s2_h_intersection_profile = supernodeS2.getIntersectionProfile(s2_h_partecipation)
            supernodeS3.updateIntersectionProfile(s2_h_partecipation, s2_h_intersection_profile)

            self.summaryHypergraph.add_node_to_edge(supernodeS3.getId(), s2_h_partecipation)

        self.summaryHypergraph.remove_node(supernodeS1.getId())
        self.summaryHypergraph.remove_node(supernodeS2.getId())

        removed_hyperedges_id = []
        for s3_h_partecipation in supernodeS3.getAdjHyperedges():
            hyperedge_object = self.summaryHypergraph.edges[s3_h_partecipation]
            if hyperedge_object.size()==1:
                self.summaryHypergraph.remove_edge(s3_h_partecipation)
                removed_hyperedges_id.append(s3_h_partecipation)
        
        for s3_h_partecipation in removed_hyperedges_id:
            supernodeS3.removeAdjHyperedge(s3_h_partecipation)

        ##...end update S3 partecipation (intersection profile) with summary hyperedges

        #Set internal_edges[supernodeS3] = internal_edges[supernodeS1] + internal_edges[supernodeS3] + edges between supernodeS1 and supernodeS3
        supernodeS3.incrementInternalEdges(supernodeS1.getInternalEdges())
        supernodeS3.incrementInternalEdges(supernodeS2.getInternalEdges())
        supernodeS3.incrementInternalEdges(supernodeS1.getWeight(supernodeS2.id))

        #Copy all the neighbours of S1 (different from S2) to S3
        for n in supernodeS1.getConnections():
            if n != supernodeS2.id:
                w = supernodeS1.getWeight(n)
                supernodeS3.addNeighbor(n, w)

                #Delete and update node n neighbourhood with respect to supernodeS1
                neighbourSupernode = self.getSupernode(n)
                neighbourSupernode.removeNeighbor(supernodeS1.id)
                neighbourSupernode.addNeighbor(supernodeS3.id, w)

        #Copy all the neighbours of S2 (different from S2) to S3: if a neighbour n has been already added in the previous step, increment the weight by 1
        for n in supernodeS2.getConnections():
            if n != supernodeS1.id:
                w = supernodeS2.getWeight(n)
                supernodeS3.updateNeighbor(n, w)

                #Delete and update node n neighbourhood with respect to supernodeS1
                neighbourSupernode = self.getSupernode(n)
                neighbourSupernode.removeNeighbor(supernodeS2.id)
                neighbourSupernode.updateNeighbor(supernodeS3.id, w)

        #Delete S1 and S2 from superList
        del self.superList[supernodeS1.id]
        del self.superList[supernodeS2.id]        

        #Return S3 id
        return supernodeS3.id
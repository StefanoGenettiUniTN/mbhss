class Supernode:
    def __init__(self, key):
        """
        id: unique identifier of the supernode
        connectedTo: dictionary{key:=id of the adjacent supernode ; value:= number of edges between the two connected supernodes}
        cardinality: number of components which populate the supernode
        internalEdges: number of connections among the components which populate the supernode
        components: set of components ids which populate the supernode
        intersection_profile: dictionary{key:= id of the relationship where the supernode partecipate ; value:= intersection profile between the supernode and the key hyperedge}
        """
        self.id = key
        self.connectedTo = {}
        self.cardinality = 0
        self.internalEdges = 0
        self.components = set()
        self.intersection_profile = {}

    def __str__(self):
        return f"supernode [{str(self.id)}] \n cardinality: {str(self.cardinality)} \n connected to (id, number_of_edges): {str([(x, self.connectedTo[x]) for x in self.connectedTo])} \n components: {str(self.components)} \n internal edges: {str(self.internalEdges)} \n intersection profile: {str([(x, self.intersection_profile[x]) for x in self.intersection_profile])} \n components: {str(self.components)}"

    def addNeighbor(self, nbr, weight=1):
        self.connectedTo[nbr] = weight

    def updateNeighbor(self, nbr, weight):
        if nbr in self.connectedTo.keys():
            self.connectedTo[nbr] += weight
        else:
            self.connectedTo[nbr] = weight
    
    def removeNeighbor(self, nbr):
        if nbr in self.connectedTo.keys():
            del self.connectedTo[nbr]
    
    def getConnections(self):
        return self.connectedTo.keys()

    def getAdjHyperedges(self):
        return self.intersection_profile.keys()

    def setIntersectionProfile(self, hyperedge, _intersection_profile):
        self.intersection_profile[hyperedge] = _intersection_profile

    def updateIntersectionProfile(self, hyperedge, _intersection_profile):
        if hyperedge in self.intersection_profile:
            self.intersection_profile[hyperedge] += _intersection_profile
        else:
            self.intersection_profile[hyperedge] = _intersection_profile

    def getIntersectionProfile(self, hyperedge):
        if hyperedge in self.intersection_profile:
            return self.intersection_profile.get(hyperedge)
        return 0

    def removeAdjHyperedge(self, hyperedge):
        if hyperedge in self.intersection_profile.keys():
            del self.intersection_profile[hyperedge]

    def getId(self):
        return self.id

    def getWeight(self, nbr):
        if nbr in self.connectedTo:
            return self.connectedTo.get(nbr)
        return 0
    
    def addComponent(self, key):
        self.cardinality += 1
        self.components.add(key)

    def getComponents(self):
        return self.components

    def getInternalEdges(self):
        return self.internalEdges

    def incrementInternalEdges(self, numEdges):
        self.internalEdges += numEdges
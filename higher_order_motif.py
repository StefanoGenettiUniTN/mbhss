#This class is used to represent a higher order motif which is found in the input hypergraph

class Higher_order_motif:
    def __init__(self, order, components, motifClass):
        """
        order : order of the motif, i.e. number of components. In our domain this can be 3 or 4
        components : tuple of supernode id which are in the motif
        motifClass   :   encodes the motif tipology. For instance ((1,2), (1,3), (2,3))
        """
        self.order = order
        self.components = components
        self.motifClass = motifClass

    def __str__(self):
        return f"[{str(self.order)}][{str(self.motifClass)}] \n supernode components: {str(self.components)}"

    
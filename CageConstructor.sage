class XEdge:
    """The end points of an edge, as an ordered pair, together with
    its "age".
    """
    def __init__(self):
        """Initialize XEdge (xtra edge)
        """
        self.ends = (0,1)
        self.age = 0

    def increment_age(self,increment=1):
        """Increment age of the XEdge. Default: increment by 1.
        """
        self.age = self.age + increment

class XGraph:
    """Xtra Graph. It is determined by the list of its XEdges, and the
    number of vertices.
    """
    def __init__(self):
        """Initialize Xtra Graph
        """
        self.edgelist=[]  # list of XEdges
        self.verts=0      # an integer. Vertices are labeled
                          # 0,1,..,self.verts
        self.untouched=[] # list of vertices such that all its incident
                          # edges are 'permanent'.
        self.g=3          # girth sougth
        self.k=3          # regularity sougth

    def graph(self):
        """Returns the underlying graph of the XGraph
        """
        g = Graph(self.verts)
        for e in self.edgelist:
            g.add_edge(e.ends)
        return g
        
def TreeForCage(n,g,k):
    """Returns an XGraph with underlying graph a k-regular tree plus
    some isolated vertices. The tree will be the base for a
    (k,g)-cage.
    
    Arguments:
    - `n`: total number of vertices
    - `g`: girth sought
    - `k`: regularity sought
    """
    def vlso(s):
        """Number of vertices up to level s. Odd case
        
        Arguments:
        - `s`: level
        """
        return (k*(k-1)^s-2)/(k-2)
    def vlse(s):
        """Number of vertices up to level s. Even case
        
        Arguments:
        - `s`: level
        """
        return 2*(((k-1)^(s+1)-1)/(k-2))

    T = XGraph()
    T.verts = n
    T.g = g
    T.k = k
    listedges = []

    if is_odd(g):
        l = (g-3)/2
        listedges = listedges + [(0,i+1) for i in range(k)]
        for i in range(l):
            for j in range(k*(k-1)^i):
                listedges = listedges +\
                    [(vlso(i)+j,vlso(i+1)+j*(k-1)+t) for t in range(k-1)]
        T.untouched = range(vlso(l))

    if is_even(g):
        l = (g-2)/2
        listedges = listedges + [(0,1)]
        for i in [-1]+range(l-1):
            for j in range(2*(k-1)^(i+1)):
                listedges = listedges +\
                    [(vlse(i)+j,vlse(i+1)+j*(k-1)+t) for t in range(k-1)]
        T.untouched = range(vlse(l-1))

    for e in listedges:
        edge = XEdge()
        edge.ends = (e[0],e[1])
        edge.age = 0
        T.edgelist.append(edge)
    
    return T

def PossibleNewEdges(G,problem='cage'):
    """List of edges that can be added to the graph.

    This is not a method of an XGraph so that we can work in different
    problems.
    
    Arguments:
    - `G`: An XGraph that should be extended
    - `problem`: name of the problem
    """
    auxg = G.graph()
    elegible_vertices = [x for x in range(G.verts) if not x in G.untouched]
    auxh = auxg.subgraph(elegible_vertices)
    edges_a_priori_elegible = auxh.complement().edges(labels=False)
    good_edges=[]
    for e in edges_a_priori_elegible:
        if problem == 'cage':
            v0 = e[0]
            v1 = e[1]
            if auxg.distance(v0,v1)>=G.g-1 and\
                    max([auxg.degree(v0),auxg.degree(v1)])<G.k:
                good_edges.append(e)
    return good_edges

        
    

        

        
        



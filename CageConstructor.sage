import random

class XEdge:
    """The end points of an edge, as an ordered pair, together with
    its "age". 

    The idea is that the "edge age" gives a measure of how long the
    edge has been present in the graph.

    whendeleted will record the run when an edge was deleted, so as to
    not include it again so soon. By now, edges never deleted will
    have this equal to -1.
    """
    def __init__(self):
        """Initialize XEdge (xtra edge)
        """
        self.ends = (0,1)
        self.age = 1
        self.whendeleted = -1

    def increment_age(self,increment=1):
        """Increment age of the XEdge. Default: increment by 1.
        """
        self.age = self.age + increment

    def pretty_print(self):
        return [self.ends,self.age,self.whendeleted]

class XGraph:
    """Xtra Graph. It is determined by the list of its XEdges, and the
    number of vertices.
    """
    def __init__(self):
        """Initialize Xtra Graph
        """
        self.edgelist=[]  # list of XEdges that could be added
        self.edgeperm=[]  # list of XEdges that have to be added
        self.verts=0      # an integer. Vertices are labeled
                          # 0,1,..,self.verts
        self.g=3          # girth sought
        self.k=3          # regularity sought

    def graph(self):
        """Returns the underlying graph of the XGraph
        """
        g = Graph(self.verts)
        for e in self.edgeperm:
            g.add_edge(e.ends)
        for e in self.edgelist:
            if e.age>0:
                g.add_edge(e.ends)
        return g

    def pretty_print_edgelist(self):
        return [e.pretty_print() for e in self.edgelist]

    def degree(self):
        """Calculates degree on the underlying graph
        """
        return self.graph().degree()

    def girth(self):
        """Calculates girth on the underlying graph
        """
        return self.graph().girth()

    def show(self):
        """Applies show to the underlying graph
        """
        return self.graph().show()       

    def show3d(self):
        """Applies show3d to the underlying graph
        """
        return self.graph().show3d(vertex_colors=\
                                       {(0.8,0.8,0.8):self.graph().vertices()},\
                                       color_by_label=True)
    
def EdgeValidInCage(G,e,g,k):
    """Checks if the edge e could be added to G and still have a
    (k,g)-graph.
    
    Arguments:
    - `G`: a Sage graph
    - `e`: a tuple
    - `g`: girth
    - `k`: regularity
    """
    return G.distance(e[0],e[1]) >= g-1 and\
        max([G.degree(e[0]),G.degree(e[1])]) < k

def IncrementAgeOfEdges(X,increment=1):
    """Returns a new XGraph where we increment the age of all XEdges
    with positive age in X.
    """
    for e in X.edgelist:
        if e.age>0:
            e.increment_age(increment)

def ResetAgeOfEdges(X):
    """Returns a new XGraph where we the age of all XEdges with
    positive age is now equal to 1.
    """
    for e in X.edgelist:
        if e.age>0:
            e.age=1

def TreeForCage(n,g,k):
    """Returns an XGraph with n vertices and underlying graph a
    k-regular tree plus some isolated vertices. The tree will be the
    base for a (k,g)-cage.
    
    Arguments:
    - `n`: total number of vertices
    - `g`: girth sought, has to be >=4.
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

    if is_odd(g):
        l = (g-3)/2
        for i in range(k):
            edge = XEdge()
            edge.ends = (0,i+1)
            edge.age = '-1'
            T.edgeperm.append(edge)
        for i in range(l):
            for j in range(k*(k-1)^i):
                for t in range(k-1):
                    edge = XEdge()
                    edge.ends = (vlso(i)+j,vlso(i+1)+j*(k-1)+t)
                    edge.age = '-1'
                    T.edgeperm.append(edge)

    if is_even(g):
        l = (g-2)/2
        edge = XEdge()
        edge.ends = (0,1)
        edge.age = '-1'
        T.edgeperm.append(edge)
        for i in [-1]+range(l-1):
            for j in range(2*(k-1)^(i+1)):
                for t in range(k-1):
                    edge = XEdge()
                    edge.ends = (vlse(i)+j,vlse(i+1)+j*(k-1)+t)
                    edge.age = '-1'
                    T.edgeperm.append(edge)

    auxg = Graph(n) # auxg will be the graph of the tree. We need it
                    # so that we can determine the edges that could
                    # never be added to the graph. Note that only
                    # possible edges are added to edgelist, and they
                    # will be the only considered in the future.
    auxg.add_edges([e.ends for e in T.edgeperm])
    nonedges = auxg.complement().edges(labels=False)

    for e in nonedges:
        if EdgeValidInCage(auxg,e,g,k):
            edge = XEdge()
            edge.ends = e
            edge.age = 0
            T.edgelist.append(edge)
    
    return T

def PossibleNewEdges(X,problem='cage'):
    """List of edges that can be added to the graph.

    This is not a method of an XGraph so that we can work in different
    problems.
    
    Arguments:
    - `X`: An XGraph that should be extended
    - `problem`: name of the problem
    """
    good_edges=[]
    edges_a_priori_elegible = \
        [edge.ends for edge in filter(lambda e:e.age==0,X.edgelist)]
    for e in edges_a_priori_elegible:
        if problem == 'cage':
            if EdgeValidInCage(X.graph(),e,X.g,X.k):
                good_edges.append(e)
    return good_edges

def FirstEdgeAvailable(X,edgelist):
    """Returns the first edge available (the first in edgelist).
    """
    return edgelist[0]

def EdgeWithDegreeSumMax(X,edgelist):
    """Returns the edge with maximum degree sum of its extremes.
    """
    def degree_sum(e):
        g = X.graph()
        e0 = e.ends[0]
        e1 = e.ends[1]
        return g.degree(e0)+g.degree(e1)

    edgewithsums = sorted(edgelist,key=degree_sum,reverse=True)
    return edgewithsums[0]    

def EdgeWithDegreeSumMaxNotRecent(X,edgelist):
    """Returns the edge with maximum degree sum of its extremes,
    taking into account the last time it was deleted.
    """
    def degree_sum(e):
        g = X.graph()
        e0 = e.ends[0]
        e1 = e.ends[1]
        return g.degree(e0)+g.degree(e1)

    edgewithsums = sorted(edgelist,key=lambda e:e.whendeleted)
    edgewithsums = sorted(edgewithsums,key=degree_sum,reverse=True)
    return edgewithsums[0]    

def ChooseDelRandomEdges(X,edgelist):
    """Choose edges randomly for deletion.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        edgs = random.sample(edgelist,i)
    else:
        edgs = edgelist
    return edgs

def ChooseDelOldRandomEdges(X,edgelist):
    """Choose edges randomly for deletion, starting with the oldest.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        oedgs = sorted(edgelist,key=lambda e:e.age,reverse=True)
        edgs = oedgs[:i]
    else:
        edgs = edgelist
    return edgs

def ChooseDelOldRandomEdgesNotRecentlyDeleted(X,edgelist):
    """Choose edges randomly for deletion, starting with the oldest,
    and considering whendeleted.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        oedgs = sorted(edgelist,key=lambda e:e.whendeleted)
        oedgs = sorted(edgelist,key=lambda e:e.age,reverse=True)
        edgs = oedgs[:i]
    else:
        edgs = edgelist
    return edgs

def EdgesCageProblem(X,edgelist):
    """Returns the edges that can be added to X in the cage
    problem ('Elegible edges').
    """
    return filter(lambda e:EdgeValidInCage(X.graph(),e.ends,X.g,X.k),\
                      edgelist)

def XGraphWithEdgeAdded(X,selectf=EdgesCageProblem,\
                            addf=EdgeWithDegreeSumMaxNotRecent):
    """Given an XGraph, returns an XGraph with one more XEdge,
    according to some method.
    
    Arguments:
    - `X`: an X graph
    - `method`: method used to add the edge
    """
    edges_a_priori_elegible = filter(lambda e:e.age==0,X.edgelist)
    found = False
    if len(edges_a_priori_elegible) > 0:
        elegible_edges = selectf(X,edges_a_priori_elegible)
        if len(elegible_edges) > 0:
            new_edge = addf(X,elegible_edges)
            found = True
        if found:
            print "Adding ",new_edge.ends
            IncrementAgeOfEdges(X)
            new_edge.age = 1

def ExtendXGraph(X,selectf=EdgesCageProblem,addf=FirstEdgeAvailable):
    """Extend an XGraph according to some method.
     
    Arguments:
    - `X`: the XGraph
    - `method`: the method used
    """
    m = X.graph().size()
    XGraphWithEdgeAdded(X,selectf,addf)
    while X.graph().size() > m:
        m = m+1
        XGraphWithEdgeAdded(X,selectf,addf)

def IsNotCageYet(X):
    """Returns True whenever the graph X with girth X.g and max degree
    less than k is not a cage.
    """
    return set(X.graph().degree()) <> set([X.k])

def SearchForGraph(X,limit=200,\
                       notdonef = IsNotCageYet,\
                       selectf = EdgesCageProblem,\
                       addf = EdgeWithDegreeSumMaxNotRecent,\
                       delf = ChooseDelOldRandomEdges):
    """Continuously look for a graph with certain properties, deleting
    edges if necessary.
    
    Arguments:
    - `X`: an XGraph
 
    - `limit`: maximum number of tries

    - `notdonef`: a boolean function that returns True if the graph
      still does not satisfy the requirement.

    - `selectf`: a function that chooses edges that can theoretically
      be added.

    - `addf`: a function that chooses which edges to actually add.

    - `delf`: a function that chooses edges to delete.

    """
    ExtendXGraph(X,selectf,addf)
    ResetAgeOfEdges(X)
    ntry = 0
    while notdonef(X) and ntry<=limit:
        ntry = ntry + 1
        print ntry
        print X.graph().degree()
        removable_edges = filter(lambda e:e.age>0,X.edgelist)
        edgs = delf(X,removable_edges)
        for e in edgs:
            e.age = 0
            e.whendeleted = ntry
            print "Removing ",e.ends
        ExtendXGraph(X,selectf,addf)

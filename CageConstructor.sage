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
        self.g=3          # girth sougth
        self.k=3          # regularity sougth

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
    """Returns a new XGraph where we increment the age of all XEdges with
    positive age in X.
    
    Arguments:
    - `X`:
    - `increment`:
    """
    for e in X.edgelist:
        if e.age>0:
            e.increment_age(increment)

def TreeForCage(n,g,k):
    """Returns an XGraph with underlying graph a k-regular tree plus
    some isolated vertices. The tree will be the base for a
    (k,g)-cage.
    
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
                    # so that we can assign an age to the nonpermanent
                    # edges.
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

def XGraphWithEdgeAdded(X,method='cage:first'):
    """Given an XGraph, returns an XGraph with one more XEdge,
    according to some method.
    
    Arguments:
    - `X`: an X graph
    - `method`: method used to add the edge
    """
    edges_a_priori_elegible = filter(lambda e:e.age==0,X.edgelist)
    found = False
    if len(edges_a_priori_elegible) > 0:
        if method == 'cage:first':
            elegible_edges =\
                filter(lambda e:EdgeValidInCage(X.graph(),e.ends,X.g,X.k),\
                           edges_a_priori_elegible)
            if len(elegible_edges) > 0:
                new_edge = elegible_edges[0]
                found = True
        elif method == 'cage:maxdegreesum':
            elegible_edges =\
                filter(lambda e:EdgeValidInCage(X.graph(),e.ends,X.g,X.k),\
                           edges_a_priori_elegible)
            if len(elegible_edges) > 0:
                edgewithsums = \
                    sorted(elegible_edges,\
                               key=lambda e:X.graph().degree(e.ends[0])+\
                               X.graph().degree(e.ends[1]),reverse=True)
                new_edge = edgewithsums[0]
                found = True
        elif method == 'cage:maxdegreesum:notrecent':
            elegible_edges =\
                filter(lambda e:EdgeValidInCage(X.graph(),e.ends,X.g,X.k),\
                           edges_a_priori_elegible)
            if len(elegible_edges) > 0:
                edgewithsums = \
                    sorted(elegible_edges,\
                               key=lambda e:X.graph().degree(e.ends[0])+\
                               X.graph().degree(e.ends[1])-e.whendeleted,reverse=True)
                new_edge = edgewithsums[0]
                found = True
        if found:
            print "Adding ",new_edge.ends
            X.edgelist.remove(new_edge)
            IncrementAgeOfEdges(X)
            new_edge.age = 1
            X.edgelist.append(new_edge)

def ExtendXGraph(X,method='cage:first'):
    """Extend an XGraph according to some method.

    For example, 'cage:first' searches for a cage always choosing the
    first available edge. It works for g=3, k=4,5,6,7,8.
    
    Arguments:
    - `X`: the XGraph
    - `method`: the method used
    """
    m = X.graph().size()
    XGraphWithEdgeAdded(X,method)
    while X.graph().size() > m:
        m = m+1
        XGraphWithEdgeAdded(X,method)

def SearchForGraph(X,limit=200,method='cage:maxdegreesum:notrecent',delmethod='oldest:random'):
    """Continuously look for a graph with certain properties, deleting
    edges if necessary.
    
    Arguments:
    - `X`: an XGraph
    - `limit`: maximum number of tries
    - `method`: method for adding edges
    - `delmethod`: method for deleting edges
    """
    ExtendXGraph(X,method)
    ntry = 0
    if method == 'cage:first' or method == 'cage:maxdegreesum'\
            or method == 'cage:maxdegreesum:notrecent':
        while set(X.graph().degree())<>set([X.k]) and ntry<=limit:
            ntry = ntry + 1
            print ntry
            print X.graph().degree()
            removable_edges = filter(lambda e:e.age>0,X.edgelist)
            if delmethod == 'random':
                i = random.choice([1,2,3])
                if i < len(removable_edges):
                    edgs = random.sample(removable_edges,i)
                else:
                    edgs = removable_edges
            if delmethod == 'oldest:random':
                i = random.choice([1,2,3])
                if i < len(removable_edges):
                    # sorted by age, starting with the oldest
                    oedgs = sorted(removable_edges,\
                                       key=lambda e: e.age,reverse=True)
                    edgs = oedgs[:i]
            for e in edgs:
                e.age = 0
                e.whendeleted = ntry
                print "Removing ",e.ends
            ExtendXGraph(X,method)

            

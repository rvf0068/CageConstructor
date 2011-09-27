import os
import random

class XEdge:
    """The end points of an edge, as an ordered pair, together with
    information on the stage at which the edge last appeared and last
    disappeared.

    whenadded is >0 if and only if the edge is present in the graph
    but the edge could be removed. In this case it indicates the stage
    at which it was added. It is =0 if the edge is not present, but
    could be added. It is -1 if the edge is 'permanent'.

    whendeleted records the stage when an edge was deleted, so as to
    not include it again so soon.

    By convention 'permanent' edges will have whenadded=-1,
    whendeleted=-1.
    """
    def __init__(self,ends=(0,1),whenadded=-1,whendeleted=-1):
        """Initialize XEdge (xtra edge)
        """
        self.ends = ends
        self.whenadded = whenadded
        self.whendeleted = whendeleted

    def pretty_print(self):
        return [self.ends,self.whenadded,self.whendeleted]

class XGraph:
    """Xtra Graph. It is determined by the list of its XEdges, and the
    number of vertices.
    """
    def __init__(self,verts,g,k,edgelist=[],edgeperm=[],pos={}):
        """Initialize Xtra Graph
        """
        self.verts=verts        # an integer. Vertices are labeled
                                # 0,1,..,self.verts
        self.g=g                # girth sought
        self.k=k                # regularity sought
        self.edgelist=edgelist  # list of XEdges that could be added
        self.edgeperm=edgeperm  # list of XEdges that have to be added
        self.pos=pos

    def graph(self):
        """Returns the underlying graph of the XGraph
        """
        g = Graph(self.verts)
        for e in self.edgeperm:
            g.add_edge(e.ends)
        for e in self.edgelist:
            if e.whenadded>0:
                g.add_edge(e.ends)
        g.set_pos(self.pos)
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

def EdgeValidInCage(X,e,g,k):
    """Checks if the XEdge e could be added to X and still have a
    (k,g)-graph.
    
    Arguments:
    - `X`: an XGraph
    - `e`: an XEdge
    - `g`: girth
    - `k`: regularity
    """
    G = X.graph()
    e0 = e.ends[0]
    e1 = e.ends[1]
    return G.distance(e0,e1) >= g-1 and max([G.degree(e0),G.degree(e1)]) < k

def TreeForCage(n,g,k):
    """Returns an XGraph with n vertices and underlying graph a
    k-regular tree plus some isolated vertices. The tree will be the
    base for a (k,g)-cage.

    Arguments:
    - `n`: total number of vertices
    - `g`: girth sought, has to be >=4.
    - `k`: regularity sought
    """

    def vls(s):
        """Amount of vertices up to level s.
        """
        if s == -1 and is_odd(g):
            return 0
        else:
            upto = [2*(((k-1)^(s+1)-1)/(k-2)),(k*(k-1)^s-2)/(k-2)]
            return upto[g % 2]

    def thepos(verts,thelen,theheigh):
        """Return the pos dictionary for the vertices vert. thelen is
        the length between vertices.
        """
        nverts = len(verts)
        u = []
        i = 0
        while i < nverts:
            x = -(nverts-1)/2*thelen+i*thelen
            u.append(x)
            i = i+1
        d = dict(zip(verts,zip(u,[theheigh]*nverts)))
        return d

    T = XGraph(n,g,k,edgelist=[],edgeperm=[],pos={})
    l = [(g-2)/2,(g-1)/2][g % 2]

    for i in range(l):
        verts = range(vls(i-1),vls(i)) # the vertices in level i
        T.pos.update(thepos(verts,2^(l-i),2*i))
        for j in range(len(verts)):
            if is_odd(g) and i == 0:
                u = k
            else:
                u = k-1
            if is_even(g):
                T.edgeperm.append(XEdge((0,1)))
            for t in range(u):
                T.edgeperm.append(XEdge((vls(i-1)+j,vls(i)+j*(k-1)+t)))
    T.pos.update(thepos(range(vls(l-1),vls(l)),1,2*l))
    T.pos.update(thepos(range(vls(l),n),1,2*(l+2)))

    # we determine the edges that could possible be added to T
    nonedges = T.graph().complement().edges(labels=False)
    for e in nonedges:
        edge = XEdge(e)
        if EdgeValidInCage(T,edge,g,k):
            edge.whenadded=0
            T.edgelist.append(edge)

    return T

def degreeSum(X,e):
    """Return the sum of the degrees of the extremums of e.
    """
    g = X.graph()
    e0 = e.ends[0]
    e1 = e.ends[1]
    return g.degree(e0)+g.degree(e1)

def First(X,edgelist,ntry):
    """Returns the first edge available (the first in edgelist).
    """
    return edgelist[0]

def DegreeSumMax(X,edgelist,ntry):
    """Returns the edge with maximum degree sum of its extremes.
    """
    edgewithsums = sorted(edgelist,\
                              key=lambda e:degreeSum(X,e),reverse=True)
    return edgewithsums[0]

def DegreeSumMaxNotRecent(X,edgelist,ntry):
    """Returns the edge with maximum degree sum of its extremes,
    taking into account the last time it was deleted.
    """
    edgewithsums = sorted(edgelist,key=lambda e:e.whendeleted)
    edgewithsums = sorted(edgewithsums,\
                              key=lambda e:degreeSum(X,e),reverse=True)
    return edgewithsums[0]

def DegreeSumMaxNotJustRemoved(X,edgelist,ntry):
    """Returns the edge with maximum degree sum of its extremes,
    taking into account the last time it was deleted.
    """
    edgewithsums = sorted(edgelist,key=lambda e:e.whendeleted)
    edgewithsums = sorted(edgewithsums,\
                              key=lambda e:degreeSum(X,e),reverse=True)
    # print "Scores. DSums:",map(lambda e:(degreeSum(X,e),e.ends),edgewithsums),"whendeleted:",\
    #     map(lambda e:e.whendeleted,edgewithsums)
    i = 0
    while i<len(edgewithsums) and edgewithsums[i].whendeleted == ntry:
        #print "vertex ",i,"whendeleted:",edgewithsums[i].whendeleted,"ntry",ntry
        #print edgewithsums[i].ends," checked."
        i = i+1
    if i<len(edgewithsums):
        return edgewithsums[i]
    else:
        print "No other choice"
        return edgewithsums[i-1]

def Random(X,edgelist,ntry=1):
    """Choose edges randomly for deletion.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        edgs = random.sample(edgelist,i)
    else:
        edgs = edgelist
    return edgs

def OldAndRandom(X,edgelist,ntry=1):
    """Choose edges randomly for deletion, starting with the oldest.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        oedgs = sorted(edgelist,key=lambda e: ntry-e.whenadded,reverse=True)
        edgs = oedgs[:i]
    else:
        edgs = edgelist
    return edgs

def OldAndRandomNotRecentlyDeleted(X,edgelist,ntry):
    """Choose edges randomly for deletion, starting with the oldest,
    and considering whendeleted.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        oedgs = sorted(edgelist,key=lambda e:e.whendeleted)
        oedgs = sorted(edgelist,key=lambda e: ntry-e.whenadded,reverse=True)
        edgs = oedgs[:i]
    else:
        edgs = edgelist
    return edgs

def AlternateDeletionMode(X,edgelist,ntry):
    """Choose edge for deletion, alternating proportional to age with
    inversely proportional to age.
    """
    i = random.choice([1,2,3])
    if i < len(edgelist):
        if is_odd(ntry):
            print "Removing old"
            oedgs = sorted(edgelist,key=lambda e: ntry-e.whenadded,reverse=True)
            print "Ages: ",map(lambda e: ntry-e.whenadded,oedgs)
            edgs = oedgs[:i]
        else:
            print "Removing recent"
            oedgs = sorted(edgelist,key=lambda e: ntry-e.whenadded)
            print "Ages: ",map(lambda e: (e.ends,ntry-e.whenadded),oedgs)
            edgs = oedgs[:i]
        return edgs

def EdgesCageProblem(X,edgelist):
    """Returns the edges that can be added to X in the cage
    problem ('Elegible edges').
    """
    return filter(lambda e:EdgeValidInCage(X,e,X.g,X.k),edgelist)

def XGraphWithEdgeAdded(X,selectf,addf,ntry):
    """Given an XGraph, returns an XGraph with one more XEdge,
    according to some method.

    Arguments:
    - `X`: an X graph
    - `selectf`: the method used to select edges
    - `addf`: the method used to add edges
    - `ntry`: the stage
    """
    edges_a_priori_elegible = filter(lambda e:e.whenadded==0,X.edgelist)
    found = False
    if len(edges_a_priori_elegible) > 0:
        elegible_edges = selectf(X,edges_a_priori_elegible)
        if len(elegible_edges) > 0:
            new_edge = addf(X,elegible_edges,ntry)
            found = True
        if found:
            print "Adding ",new_edge.ends,"out of",len(elegible_edges)
            new_edge.whenadded = ntry

def ExtendXGraph(X,selectf,addf,ntry):
    """Extend an XGraph according to some method.

    Arguments:
    - `X`: the XGraph
    - `selectf`: the method used to select edges
    - `addf`: the method used to add edges
    - `ntry`: the stage
    """
    m = X.graph().size()
    XGraphWithEdgeAdded(X,selectf,addf,ntry)
    while X.graph().size() > m:
        m = m+1
        XGraphWithEdgeAdded(X,selectf,addf,ntry)

def IsNotCageYet(X):
    """Returns True whenever the graph X with girth X.g and max degree
    less than k is not a cage.
    """
    return set(X.graph().degree()) <> set([X.k])

OneEdgeF = (First,DegreeSumMax,DegreeSumMaxNotRecent,\
                DegreeSumMaxNotJustRemoved)

EdgesF = (Random,OldAndRandom,OldAndRandomNotRecentlyDeleted,\
              AlternateDeletionMode)

def SearchForGraph(X,limit=200,\
                       notdonef = IsNotCageYet,\
                       selectf = EdgesCageProblem,\
                       addf = OneEdgeF[3],\
                       delf = EdgesF[2]):
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
    ntry = 1
    ExtendXGraph(X,selectf,addf,ntry)
    while notdonef(X) and ntry<=limit:
        print ntry
        print X.graph().degree()
        removable_edges = filter(lambda e:e.whenadded>0,X.edgelist)
        edgs = delf(X,removable_edges,ntry)
        for e in edgs:
            e.whenadded = 0
            e.whendeleted = ntry
            print "Removing ",e.ends
        ExtendXGraph(X,selectf,addf,ntry)
        ntry = ntry + 1
    if notdonef(X):
        icon = '/usr/share/icons/gnome/48x48/emotes/face-crying.png'
        os.system("notify-send --icon "+icon+" 'Could not find a suitable graph'")
    else:
        icon = '/usr/share/icons/gnome/48x48/emotes/face-laugh.png'
        os.system("notify-send --icon "+icon+" 'Found a suitable graph!'")

# Local Variables:
# eval: (yas/minor-mode)
# End:

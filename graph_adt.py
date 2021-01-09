# modulo graph_adt

def newG(V):
    return [V,[]]

def add_edgesG(g,E):
    def edge(e):
        return len(e)==2 and e[0] in g[0] and e[1] in g[0]
    if all(map(edge,E)):
        return [g[0],E]
    else:
        return g
    
def get_nodesG(g):
    return g[0]

def get_edgesG(g):
    return g[1]

def emptyG(g):
    return g[1]==[]
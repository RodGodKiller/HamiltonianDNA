# modulo path_adt

def newP(v1,v2):
    return [v1,v2]

def firstP(p):
    return p[0]

def lastP(p):
    return p[-1]

def lenP(p):
    assert(len(p)>=2)
    return len(p)

def compP(p1,p2):
    return p1[-1]==p2[0]

def glueP(p1,p2):
    if compP(p1,p2):
        return p1+p2[1:]
    else:
        print("Warning (glueP). Caminhos incompatíveis.",p1,p2)
        return p1
    
def breakP(p,i): # breaks a path at a cycle
    assert 0<i<len(p)-1
    return p[:i+1],p[i:]
    
def crossesP(p,v):
    return v in p

def copyP(p):
    return p[:]

def showP(p):
    print(p)
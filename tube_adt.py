#modulo tube_adt

import path_adt as path
import ind_adt as ind

def newT():
    return {}

def existsT(t,i):
    return i in t.keys()

def addT(t,i):
    t[ind.idI(i)]=i
    return t

def delT(t,i):
    t.pop(i,None)
    return t

def get_indT(t,i):
    return t.get(i)

def get_idsT(t):
    return list(t.keys())

def neigT(t,ii):
    def daux(i):
        return ind.distI(i,i1)
    if ii in t.keys():
        i1=t[ii]
        cand=[ix for ix in t.values() if path.compP(ind.pathI(i1),ind.pathI(ix)) or path.compP(ind.pathI(ix),ind.pathI(i1))]
        cand.sort(key=daux)
        return cand[:5]
    else:
        return []
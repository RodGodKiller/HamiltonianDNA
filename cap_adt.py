# modulo cap_adt

import event_adt as event

def newS():
    return []

def nextS(s):
    assert len(s)>0
    return s[0]

def delS(s):
    assert len(s)>0
    return s[1:]

def addS(s,ev):
    i=0
    j=len(s)-1
    ev_t = event.timeE(ev)
    while (i<=j):
        m=(i+j)//2
        if ev_t < event.timeE(s[m]):
            j=m-1
        else:
            i=m+1
    return s[:j+1]+[ev]+s[j+1:]

def lenS(s):
    return len(s)

def showS(s):
    for e in s:
        print(event.timeE(e),event.kindE(e),event.idE(e))
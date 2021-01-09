# módulo ind_adt

from math import sqrt

def newI(p,x,y,z,i):
    return [p,x,y,z,i]

def pathI(i):
    return i[0]

def xposI(i):
    return i[1]

def yposI(i):
    return i[2]

def zposI(i):
    return i[3]

def idI(i):
    return i[4]

def distI(i1,i2):
    return sqrt((i1[1]-i2[1])**2+(i1[2]-i2[2])**2+(i1[3]-i2[3])**2)

def eqI(i1,i2):
    return i1[-1]==i2[-1]
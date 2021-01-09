import cap_adt as cp
import event_adt as ev
import graph_adt as gr
import ind_adt as ind
import path_adt as path
import tube_adt as tb
from random import random, choice, uniform, randint
from math import log
from itertools import islice


def random_chunk(li, min_chunk=1, max_chunk=3):
    it = iter(li)
    while True:
        nxt = list(islice(it, randint(min_chunk, max_chunk)))
        if nxt:
            yield nxt
        else:
            break


def exprandom(m):
    x = random()
    return -m*log(x)


def simu(G, o = 1, d = 4, k = 20, ht = 20, tbc = 1, tbd = 1, tbz = 1):
    nodeList=gr.get_nodesG(G)

    # 0 - PREPARAÇÃO

    # Cria-se uma solução e uma CAP novas (a CAP com um evento de cisão.)
    sol=tb.newT()
    cap=cp.newS()
    cap=cp.addS(cap, ev.newE(exprandom(tbz), "cis", 0))

    # Criam-se os indivíduos (k cópias de cada uma das arestas de g) e adicionam-se à solução
    id = 0
    for e in gr.get_edgesG(G):
        for _ in range(k):
            id += 1
            i = ind.newI(e, random(), random(), random(), id)
            sol = tb.addT(sol, i)

    # Adiciona-se um evento de concatenação e deslocamento por cada indivíduo à CAP.
    for x in tb.get_idsT(sol):
        cap = cp.addS(cap, ev.newE(exprandom(tbc), "conc", x))
        cap = cp.addS(cap, ev.newE(exprandom(tbd), "desl", x))

    # 1 - HIBRIDAÇÃO
    time = 0
    evt = cp.nextS(cap)
    cap = cp.delS(cap)
    
    while time < ht:
        id = ev.idE(evt)
        
        # Se o indivíduo correspondente ao evento atual já não estiver na solução, passar à frente.
        if tb.existsT(sol, id):

            # Concatenação
            if ev.kindE(evt) == "conc":
                viz = tb.neigT(sol, id)
                indi = tb.get_indT(sol, id)
                # Se houver indivíduos compatíveis, proceder com a concatenação
                if viz != []:
                    # Escolhe um de entre os cinco disponíveis.
                    chosen = choice(viz) 
                    if path.compP(ind.pathI(indi), ind.pathI(chosen)):
                        caminho = path.glueP(ind.pathI(indi), ind.pathI(chosen))
                    else:
                        caminho = path.glueP(ind.pathI(chosen), ind.pathI(indi))
                    # Cria-se o novo indivíduo com identificador +1 que o último da solução.
                    newInd = ind.newI(
                                        caminho, 
                                        (ind.xposI(indi)+ind.xposI(chosen))/2, 
                                        (ind.yposI(indi)+ind.yposI(chosen))/2, 
                                        (ind.zposI(indi)+ind.zposI(chosen))/2, 
                                        tb.get_idsT(sol)[-1] + 1
                                    )
                    sol = tb.addT(sol, newInd)
                    # Removem-se os dois indivíduos da solução
                    sol = tb.delT(sol, ind.idI(chosen))
                    sol = tb.delT(sol, ind.idI(indi))

                    # Adiciona-se um novo evento de concatenação e deslocamento para este novo indivíduo.
                    cap = cp.addS(cap, ev.newE(time + exprandom(tbc), "conc", ind.idI(newInd)))
                    cap = cp.addS(cap, ev.newE(time + exprandom(tbd), "desl", ind.idI(newInd)))

                # Se não houver, não se faz nada. 
                else:
                    pass

            # Deslocamento
            elif ev.kindE(evt) == "desl":
                indi = tb.get_indT(sol, id)
                c = path.lenP(ind.pathI(indi))
                dist = 1/c
                x = ind.xposI(indi)
                y = ind.yposI(indi)
                z = ind.zposI(indi)
                newInd = ind.newI(
                                    ind.pathI(indi), 
                                    uniform(x-dist,x+dist), 
                                    uniform(y-dist,y+dist), 
                                    uniform(z-dist,z+dist), 
                                    ind.idI(indi)
                                )
                # Adiciona-se um novo evento de deslocação para este indivíduo

                cap = cp.addS(cap, ev.newE(time + exprandom(tbd), "desl", ind.idI(newInd)))


            # Cisão
            elif ev.kindE(evt) == "cis":
                # Percorrem-se todos os indivíduos
                for i in tb.get_idsT(sol):
                    cami = ind.pathI(tb.get_indT(sol, i))
                    if path.lenP(cami) > len(nodeList):
                        # Cria-se uma lista das partições do caminho
                        paths = list(random_chunk(li, 1, len(nodeList)//2))
                        # Adicionam-se os indivíduos criados da partição dos caminhos.
                        for p in paths:
                            sol = tb.addT(sol, ind.newI(p, 
                                                         random(),
                                                         random(),
                                                         random(),
                                                         tb.get_idsT(sol)[-1]
                                                       )
                                         )
                        # Remove-se o indivíduo da solução
                        sol = tb.delT(sol, ind.idI(i))

                cap = cp.addS(cap, ev.newE(time + exprandom(tbz), "cis", 0))
        cap = cp.delS(cap)
        evt = cp.nextS(cap)
        time = ev.timeE(evt)

# 2 - CAMINHOS COM ORIGEM E DESTINO PRETENDIDOS
    for x in tb.get_idsT(sol):
        oT = ind.pathI(tb.get_indT(sol, x))[0]  
        dT = ind.pathI(tb.get_indT(sol, x))[-1]
        if oT != o or dT != d:
            sol = tb.delT(sol, x)

# 3 - CAMINHOS COM COMPRIMENTO PRETENDIDO
    for x in tb.get_idsT(sol):
        comp = len(ind.pathI(tb.get_indT(sol, x)))
        if comp != len(nodeList):
            sol = tb.delT(sol, x)

# 4 - CAMINHOS QUE PASSAM EM TODOS OS VÉRTICES
    for v in nodeList:
        for x in tb.get_idsT(sol):
            if not path.crossesP(ind.pathI(tb.get_indT(sol, x)), v):
                sol = tb.delT(sol, x)

# 5 - RESULTADO
    return list(set([ind.pathI(i) for i in [tb.get_indT(sol, x) for x in tb.get_idsT(sol)]]))


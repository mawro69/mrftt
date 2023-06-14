#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:02:40 2022

@author: abouabdallah
"""

import numpy as np
import tt 

import bases





#same init for all

def init_Bi_Ai(factors):
    """ initialise the shapes, modes and the cores of the tt matrices Ai
    arguments :
        factors : list, list of factors (length = m)
    returns :
        factorRankArray :   matrix, the ranks of the cores 
        factorModeSizeArray : vector, the modes of the cores 
        factorCoreArray : list of list, list of the cores 

    """
    # number of nores and factors
    n = len(factors[0])
    m = len(factors)

    factorRankArray = np.ones((n+1, m))
    # tableau des modes des tenseurs de facteurs
    factorModeSizeArray = np.zeros(n)  # stock the mode size
    # list to stock the cores of the futures Ais list of list of sape n
    factorCoreArray = bases.listoflist(n)
    for i in range(n):
        # for every Ais we need to initalize a list of n cores
        factorCoreArray[i] = bases.listoflist(m)

    for l in range(m):  # for each \ell in 1 : m
        for i in range(n):
            Gil = factors[l][i]  # extracte G_i^\ell
            [r1, n1, r2] = Gil.shape  # compute the size of  de G_i^\ell
            factorRankArray[i, l] = r1  # stock r_i
            factorModeSizeArray[i] = n1  # stock the modes of the tenseurs
            factorCoreArray[i][l] = Gil  # stock the G_i\ell for each \ell

    return (factorRankArray, factorModeSizeArray, factorCoreArray)



# ======================================================================




""" without cores assembly """
""" functions rhat not depend on the core fusion """

# ======================================================================


def stock_Ai(factors, i, di):
    """compute the TT-matrices  Ai = Kron (G_i^l)
    arguments : 
        factorCoreArray : list of list, each element factorCoreArray[i] give as the cores of tt matrice Ai[Zi]
        i : int, the index of Ai
        di= int, the mode of Ai ie Z_i
    return :
        Ai : tt-matrix : """
    factorRankArray, factorModeSizeArray, factorCoreArray = init_Bi_Ai(
        factors)  # initialize the differents mode, shapes and G_{i}^\ell
    m = len(factorCoreArray[0])

    rowRankArray = factorRankArray[i, :]  # stock the differents r_i of Ai[z_i]
    colRankArray = factorRankArray[i + 1, :]
    Ai = tt.ones(bases.multi(rowRankArray, colRankArray))
    Ai = tt.matrix(Ai, n=rowRankArray, m=colRankArray)
    coresAi = tt.matrix.to_list(Ai)  # create the TT-matrix Ai[z_i]
    for l in range(m):
        n1, n2 = factorCoreArray[i][l][:, di, :].shape
        coresAi[l] = np.reshape(factorCoreArray[i][l]
                                [:, di, :], (1, n1, n2, 1))  # reshape the cores to create Ai[z_i]
        # as list of G_i[z_i]

    Ai = Ai.from_list(coresAi)  # create the TT-matrix Ai

    return Ai

# ======================================================================

# ======================================================================


def stock_all_Ai(factors, di, n):
    """ compute all the TT-matrices  Ai[z_i] for i =1, ..., n
    arguments :
        factorCoreArray : list of list, each element factorCoreArray[i] give as the cores of tt matrice Ai[Zi]
        n : int, number of individuals 
        di= int, the mode of Ai ie Z_i
    return :
        listofAi : list of Ai[z_i] for i =1, ..., n

    """
    listofAi = []  # empty list
    for i in range(n):  # for i=0, ..., n-1
        listofAi.append(stock_Ai(factors, i, di))  # listofAi [i] = Ai[Z_i]
    return listofAi




"""then uniform compression """


# uniform compression of cores of Ai

# ======================================================================


# ======================================================================


def how_comp(n, k):
    """ get the number of cores after fusion :
        argument :
            m : int, number of links 
            k : number of compressed cores 
        return  :
            numcomp : vector, the number of corempressed cores for the new cores   
            comp :  int, number of new cores 

    """
    m=bases.get_m(n)
    lis = bases.getFactors(m)
    if (k in lis):
        numcomp = int(m/k)
        comp = k*np.ones(numcomp)
    else:
        return("error")
    return (numcomp, comp)

# ====================================================================== 


# coding maybe boring, commenting your code  is more devasting than the 12 season of The Big Bang Theory

def to_compress(comp, numcomp):
    """ give the index of the core that will be merged 
    """
    tobecomp = []#list of future core merged 
    tobecomp = bases.listoflist(numcomp) 
    a = 0
    for k in range(numcomp):
        tobecomp[k] = []
        for l in range(int(comp[k])):
            tobecomp[k].append(a)
            a += 1
    return tobecomp

# ======================================================================


def compressed_cores(factorCoreArray, di, i, futurcoresAi, comp, numcomp):
    tobecomp = to_compress(comp, numcomp)
    facs = factorCoreArray[i]
    for k in range(numcomp):
        n1, n2 = facs[tobecomp[k][0]][:, di, :].shape
        futurcoresAi[k] = np.reshape(
            facs[tobecomp[k][0]][:, di, :], (1, n1, n2, 1))
        for j in range(1, len(tobecomp[k])):
            n1, n2 = facs[tobecomp[k][j]][:, di, :].shape
            a = np.reshape(facs[tobecomp[k][j]][:, di, :], (1, n1, n2, 1))
            futurcoresAi[k] = np.kron(futurcoresAi[k], a)
    return futurcoresAi


# ======================================================================


def stock_Ai_compressed(factors, i, di, comp, numcomp):
    """calcul des matrices  Ai = Kron (G_i^l)
    arguments : 
        factorCoreArray list of list: Une liste de liste ou chaque factorCoreArray[i] nous donne 
        les éléments qui constituent la tt matrice Ai   
    return :
        Bi matrix : """
    factorRankArray, factorModeSizeArray, factorCoreArray = init_Bi_Ai(factors)
    #m = len(factorCoreArray[0])
    #length = int(m/numcomp)
    # stock les différents r_i jouant les numéros de lignes
    rowRankArray = factorRankArray[i, :]
    # %stock les différents r_i jouant les numéros de colonnes
    colRankArray = factorRankArray[i + 1, :]
    rowRankArraycomp = np.ones(numcomp)
    colRankArraycomp = np.ones(numcomp)
    res = to_compress(comp, numcomp)
    for k in range(numcomp):
        for j in range(len(res[k])):
            rowRankArraycomp[k] = rowRankArraycomp[k]*rowRankArray[res[k][j]]
            colRankArraycomp[k] = colRankArraycomp[k]*colRankArray[res[k][j]]
    Ai = tt.ones(bases.multi(rowRankArraycomp, colRankArraycomp))
    Ai = tt.matrix(Ai, n=rowRankArraycomp, m=colRankArraycomp)
    futurcoresAi = tt.matrix.to_list(Ai)
    futurcoresAi = compressed_cores(
        factorCoreArray, di, i, futurcoresAi, comp, numcomp)
    Ai = Ai.from_list(futurcoresAi)

    return Ai



# ======================================================================


# ======================================================================

def stock_all_Ai_compressed(factors, di, n, comp, numcomp):
    listofAi = []
    for i in range(n):
        listofAi.append(stock_Ai_compressed(factors, i, di, comp, numcomp))
    return listofAi







""" the second uniform compression method 
this second way legit the fact that m/k may not equal zeros
"""


def how_comp2(n, k):
    """ get the number of cores after fusion :
        argument :
            m : int, number of links 
            k : number of compressed cores 
        return  :
            numcomp : vector, the number of corempressed cores for the new cores   
            comp :  int, number of new cores 

    """
    m=bases.get_m(n)
    lis = bases.getFactors(m)
    numcomp, comp=0,0
    if (k in lis):
        numcomp = int(m/k)
        comp = k*np.ones(numcomp)
    else:
        numcomp = int(m/k) +1
        comp = k*np.ones(numcomp)
        add_core= m-k*(numcomp-1) 
        comp[numcomp-1] = add_core
    return (numcomp, comp)


# ====================================================================== 


# coding maybe boring, commenting your code  is more devasting than the 12 season of The Big Bang Theory

def to_compress2(comp, numcomp):
    """ give the index of the core that will be merged 
    """
    tobecomp = []#list of future core merged 
    tobecomp = bases.listoflist(numcomp) 
    a = 0
    for k in range(numcomp):
        tobecomp[k] = []
        for l in range(int(comp[k])):
            tobecomp[k].append(a)
            a += 1
    return tobecomp

# ======================================================================


def compressed_cores2(factorCoreArray, di, i, futurcoresAi, comp, numcomp):
    tobecomp = to_compress2(comp, numcomp)
    facs = factorCoreArray[i]
    for k in range(numcomp):
        n1, n2 = facs[tobecomp[k][0]][:, di, :].shape
        futurcoresAi[k] = np.reshape(
            facs[tobecomp[k][0]][:, di, :], (1, n1, n2, 1))
        for j in range(1, len(tobecomp[k])):
            n1, n2 = facs[tobecomp[k][j]][:, di, :].shape
            a = np.reshape(facs[tobecomp[k][j]][:, di, :], (1, n1, n2, 1))
            futurcoresAi[k] = np.kron(futurcoresAi[k], a)
    return futurcoresAi


# ======================================================================


def stock_Ai_compressed2(factors, i, di, comp, numcomp):
    """calcul des matrices  Ai = Kron (G_i^l)
    arguments : 
        factorCoreArray list of list: Une liste de liste ou chaque factorCoreArray[i] nous donne 
        les éléments qui constituent la tt matrice Ai   
    return :
        Bi matrix : """
    factorRankArray, factorModeSizeArray, factorCoreArray = init_Bi_Ai(factors)
    #m = len(factorCoreArray[0])
    #length = int(m/numcomp)
    # stock les différents r_i jouant les numéros de lignes
    rowRankArray = factorRankArray[i, :]
    # %stock les différents r_i jouant les numéros de colonnes
    colRankArray = factorRankArray[i + 1, :]
    rowRankArraycomp = np.ones(numcomp)
    colRankArraycomp = np.ones(numcomp)
    res = to_compress2(comp, numcomp)
    for k in range(numcomp):
        for j in range(len(res[k])):
            rowRankArraycomp[k] = rowRankArraycomp[k]*rowRankArray[res[k][j]]
            colRankArraycomp[k] = colRankArraycomp[k]*colRankArray[res[k][j]]
    Ai = tt.ones(bases.multi(rowRankArraycomp, colRankArraycomp))
    Ai = tt.matrix(Ai, n=rowRankArraycomp, m=colRankArraycomp)
    futurcoresAi = tt.matrix.to_list(Ai)
    futurcoresAi = compressed_cores(
        factorCoreArray, di, i, futurcoresAi, comp, numcomp)
    Ai = Ai.from_list(futurcoresAi)

    return Ai



# ======================================================================


# ======================================================================

def stock_all_Ai_compressed2(factors, di, n, comp, numcomp):
    listofAi = []
    for i in range(n):
        listofAi.append(stock_Ai_compressed2(factors, i, di, comp, numcomp))
    return listofAi




"""merge with uniformes volumes """


# ======================================================================

def build_symbolic_ttfactors(n, edges):
    """ calcul les gi et les dimensions non essentielles
    arguments :
        @modele dict : notre modèle stocké dans un dictionnaires
          @ A_\ell list of integers : ensembles des noeuds qui composent notre modèle
          @listedespsislist : ensemble des psi qui constituent notre modèle d'ising
        @wd  int: \in S vaut -1 ou +1 pour le modèle d'ising
    return :
        @ ttfacteurs : liste de liste  remplie avec les différents G_i"""

    m = edges.shape[0]
    facs = [[] for i in range(m)]
    for l in range(m):
        cfacs = ['empty' for i in range(n)]
        i = int(edges[l, 0])
        j = int(edges[l, 1])
        # print(l, list_of_edges[l], i, j)
        for k in range(i):
            cfacs[k] = '1'
        for k in range(j+1, n):
            cfacs[k] = '1'
        cfacs[i] = 'U'
        for k in range(i+1, j):
            cfacs[k] = 'I'
        cfacs[j] = 'V'
        facs[l] = cfacs
    return facs



# ======================================================================

def loadprod(load) :
    prod =0
    for k in range (len(load)) :
        prod= prod+np.log(load[k])
    return (prod)

# ======================================================================

def getFactors(n):
    # Create an empty list for factors
    factors=[];

    # Loop over all factors
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)

    # Return the list of factors
    return factors


# ======================================================================


def loadprod(load):
    """ Compute the log of volume of a matrix in tt format
        Keyword arguments:
            load -- list of volumes of cores 
        Return :
            prod : the log of product 
    """
    prod = 0
    for k in range(len(load)):
        prod = prod+np.log(load[k])
    return (prod)

# ======================================================================

def list_to_vec(liste):
    """ Transform a list to vector
        Keyword arguments:
            liste : liste
        Return :
            vec= the vector
            """
    f = len(liste)
    vec = np.zeros(f)
    for k in range(f):
        vec[k] = liste[k]
    return vec

# ======================================================================

def symbolic_factors(dic):
    """ Cpmpute symbolic factors from the dictionnary of graph
    Keyword arguments:
            dic -- dictionnary of edges 
        Return :
            sfactors : """
    list_of_edges = list(dic.keys())
    edges = np.zeros([len(dic.keys()), 2], dtype=int)
    m = len(list_of_edges)
    n = max(max(list_of_edges))+1

    #nodes_edges = [[] for l in range(n)]

    n = max(max(list_of_edges))+1
    for l in range(m):
        i, j = list_of_edges[l]
        edges[l, 0] = i
        edges[l, 1] = j

    sfactors = build_symbolic_ttfactors(n, edges)
    return sfactors

# ======================================================================

def nds_edges(dic):
    """ Fonction qui donne la nouvelle numérotation des liens \ell
        Keyword arguments:
            dic : dictionnaire du graphe des liens 
        Return :
            edges : liste des liens """
    list_of_edges = list(dic.keys())
    edges = np.zeros([len(dic.keys()), 2], dtype=int)
    m = len(list_of_edges)
    n = max(max(list_of_edges))+1

    nodes_edges = [[] for l in range(n)]

    for l in range(m):
        i, j = list_of_edges[l]
        edges[l, 0] = i
        edges[l, 1] = j
        nodes_edges[i].append(l)
        nodes_edges[j].append(l)
        # Compute the symbolic factors

    return n, m, nodes_edges, edges

# ======================================================================

def comp_load(sfactors, dic):

    n, m, nodes_edges, edges = nds_edges(dic)

    sA = [[] for i in range(n)]
    for i in range(n):
        cA = ['' for i in range(m)]
        for l in range(m):
            cA[l] = sfactors[l][i]
        sA[i] = cA
    comp = np.zeros([n, 4], dtype=int)
    load = np.zeros([n, n-1], dtype=int)
    for l, e in enumerate(sA):
        comp[l, 0] = e.count('1')
        comp[l, 1] = e.count('I')
        comp[l, 2] = e.count('U')
        comp[l, 3] = e.count('V')
        #print("A[", l, "] ", e)
        for k in range(n-1):
            e = nodes_edges[l][k]
            load[l, k] = edges[e, 1]-edges[e, 0]

    return load, comp, sA

# ======================================================================

def volume_cores_A_i(k, sA, m, q):
    load2 = np.zeros(m)
    for l, e in enumerate(sA[k]):
        if e == '1':
            load2[l] = np.log(1)
        elif e == 'I':
            load2[l] = np.log(q*q)
        elif e == 'U':
            load2[l] = np.log(q)
        elif e == 'V':
            load2[l] = np.log(q)
        else:
            #print("Wrong type ", e)
            exit(-1)

    return load2

# ======================================================================

def init_groupage_Ai(load2, nb_paquets):

    size = load2.sum()/nb_paquets
#    paquets = np.zeros(nb_paquets+2, dtype=int)
#    paquets_load = np.zeros(nb_paquets+1, dtype=int)
    paquets2 = []
    paquets_load2 = []
    count = 1
    loadp = 0
    for l in range(load2.size):
        loadp += load2[l]
       # aa=0
        if(loadp > size):

            paquets2.append(l-1)
            paquets_load2.append(loadp-load2[l])
            loadp = load2[l]
            count += 1
            # print(aa)
    paquets2.append(load2.size)
    paquets_load2.append(loadp)
    paquets_load2v = list_to_vec(paquets_load2)
    paquets2v = list_to_vec(paquets2)

    return paquets_load2v, paquets2v


# ======================================================================

def groupage_Ai_main(load2, nb_paquets):
    paquets_load2v, paquets2v = init_groupage_Ai(load2, nb_paquets)
    paquets2v1 = np.zeros(len(paquets2v))
    paquets2v1[0] = paquets2v[0]
    for k in range(1, len(paquets2v)):
        paquets2v1[k] = int(paquets2v[k]-paquets2v[k-1])
        paquets2v1[k] = int(paquets2v1[k])

    return paquets2v1

# ======================================================================

def volumes_groupes(load2, groupes):
    volumes = np.ones(len(groupes))
    a = 0
    for k in range(len(groupes)):
        volumes[k] = 0

        for l in range(int(groupes[k])):
            a = a+1
            volumes[k] = volumes[k]+load2[a]

    return volumes


# ======================================================================

def groupage_Ai_all(groupes, nmax):
    nvgp = []
    for k in range(len(groupes)):
        # print(k)
        if (groupes[k] < nmax+1):
            nvgp.append(groupes[k])
        else:
            a = int(groupes[k]//nmax)
            pn = int(groupes[k] % nmax)
            for k in range(a):
                nvgp.append(nmax)
            if (pn > 0):
                nvgp.append(pn)

    gpnv = list_to_vec(nvgp)

    return gpnv

# ======================================================================

def multi(vec1, vec2):
    n = len(vec1)
    vec3 = np.ones(n)
    for k in range(n):
        vec3[k] = vec1[k]*vec2[k]
    return vec3

# ======================================================================

def symbolic_product(Bis3):
    lis = []
    n = len(Bis3)
    for k in range(n-1, -1, -2):
        #print([k, k-1])
        Bk = Bis3[k].tt.r
        Bk1 = Bis3[k-1].tt.r
        aa = multi(Bk, Bk1)
        lis.append(aa)
    shapes = lis
    return shapes


# ======================================================================

def how_comp_bis(dic, nb_paquets, nmax , q) :#default nb_paquets = 300 and nmax=4
    sfactors=symbolic_factors(dic)
    n, m, nodes_edges, edges=nds_edges(dic)
    load, comp, sA =comp_load(sfactors, dic)
    max_ind = np.argmax(comp[:, 1])
    load2=volume_cores_A_i(max_ind,sA , m, q)
    paquets_load2v, paquets2v=init_groupage_Ai (load2, nb_paquets)
    groupes=groupage_Ai_main (load2, nb_paquets)
    bgroupes=groupage_Ai_all (groupes, nmax)
    len(bgroupes)
    numcorps=len(bgroupes)


    return bgroupes, numcorps 

"""
# ======================================================================


def loadprodadd1(load):
    prod = 0
    for k in range(len(load)):
        prod = prod+np.log(load[k]+1)
    return (prod)



# ======================================================================








def comp_load(sfactors, dic):

    n, m, nodes_edges, edges = nds_edges(dic)

    sA = [[] for i in range(n)]
    for i in range(n):
        cA = ['' for i in range(m)]
        for l in range(m):
            cA[l] = sfactors[l][i]
        sA[i] = cA
    comp = np.zeros([n, 4], dtype=int)
    load = np.zeros([n, n-1], dtype=int)
    for l, e in enumerate(sA):
        comp[l, 0] = e.count('1')
        comp[l, 1] = e.count('I')
        comp[l, 2] = e.count('U')
        comp[l, 3] = e.count('V')
        print("A[", l, "] ", e)
        for k in range(n-1):
            e = nodes_edges[l][k]
            load[l, k] = edges[e, 1]-edges[e, 0]

    return load, comp, sA


def volume_cores_A_i(k, sA, m, q):
    load2 = np.zeros(m)
    for l, e in enumerate(sA[k]):
        if e == '1':
            load2[l] = np.log(1)
        elif e == 'I':
            load2[l] = np.log(q*q)
        elif e == 'U':
            load2[l] = np.log(q)
        elif e == 'V':
            load2[l] = np.log(q)
        else:
            print("Wrong type ", e)
            exit(-1)

    return load2


def init_groupage_Ai(load2, nb_paquets):

    size = load2.sum()/nb_paquets
#    paquets = np.zeros(nb_paquets+2, dtype=int)
#    paquets_load = np.zeros(nb_paquets+1, dtype=int)
    paquets2 = []
    paquets_load2 = []
    count = 1
    loadp = 0
    for l in range(load2.size):
        loadp += load2[l]
       # aa=0
        if(loadp > size):

            paquets2.append(l-1)
            paquets_load2.append(loadp-load2[l])
            loadp = load2[l]
            count += 1
            # print(aa)
    paquets2.append(load2.size)
    paquets_load2.append(loadp)
    paquets_load2v = list_to_vec(paquets_load2)
    paquets2v = list_to_vec(paquets2)

    return paquets_load2v, paquets2v


def groupage_Ai_main(load2, nb_paquets):
    paquets_load2v, paquets2v = init_groupage_Ai(load2, nb_paquets)
    paquets2v1 = np.zeros(len(paquets2v))
    paquets2v1[0] = paquets2v[0]
    for k in range(1, len(paquets2v)):
        paquets2v1[k] = int(paquets2v[k]-paquets2v[k-1])
        paquets2v1[k] = int(paquets2v1[k])

    return paquets2v1


def volumes_groupes(load2, groupes):
    volumes = np.ones(len(groupes))
    a = 0
    for k in range(len(groupes)):
        volumes[k] = 0

        for l in range(int(groupes[k])):
            a = a+1
            volumes[k] = volumes[k]+load2[a]

    return volumes


def groupage_Ai_all(groupes, nmax):
    nvgp = []
    for k in range(len(groupes)):
        # print(k)
        if (groupes[k] < nmax+1):
            nvgp.append(groupes[k])
        else:
            a = int(groupes[k]//nmax)
            pn = int(groupes[k] % nmax)
            for k in range(a):
                nvgp.append(nmax)
            if (pn > 0):
                nvgp.append(pn)

    gpnv = list_to_vec(nvgp)

    return gpnv

"""
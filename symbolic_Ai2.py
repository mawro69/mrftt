#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 07:36:37 2022

@author: abouabdallah
"""


# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 14:05:28 2022

@author: abouabdallah
"""




import numpy as np
import tt


#import matplotlib.pyplot as plt
#import compressed_mrf_tt
#import symbolic_tt_mrf


def build_symbolic_ttfactors(n, edges):

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




def loadprod(load) :
    prod =0
    for k in range (len(load)) :
        prod= prod+np.log(load[k])
    return (prod)


def getFactors(n):
    # Create an empty list for factors
    factors=[];

    # Loop over all factors
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)

    # Return the list of factors
    return factors



def getFactors(n):
    # Create an empty list for factors
    factors = []

    # Loop over all factors
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)

    # Return the list of factors
    return factors


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


def multi(vec1, vec2):
    n = len(vec1)
    vec3 = np.ones(n)
    for k in range(n):
        vec3[k] = vec1[k]*vec2[k]
    return vec3


def symbolic_product(Bis3):
    lis = []
    n = len(Bis3)
    for k in range(n-1, -1, -2):
        print([k, k-1])
        Bk = Bis3[k].tt.r
        Bk1 = Bis3[k-1].tt.r
        aa = multi(Bk, Bk1)
        lis.append(aa)
    shapes = lis
    return shapes


def compute_Bis_bysteps(roundprec, Bis3):
    lis = []
    print('hi')
    n = len(Bis3)
    if (n % 2 == 1):
        for k in range(n-1, 0, -2):
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk1*Bk
            aa = aa.round(roundprec)
            lis.append(aa)
        lis.append(Bis3[0])
        step1 = lis
    else:
        lis = []
        n = len(Bis3)
        for k in range(n-1, -1, -2):
            print('hi')
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk1*Bk
            aa = aa.round(roundprec)
            lis.append(aa)
        step1 = lis

    return step1


def compute_Bis_bysteps_imp(roundprec, Bis3):
    lis = []
    print('hi')
    n = len(Bis3)
    if (n % 2 == 1):
        for k in range(n-1, 0, -2):
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk*Bk1
            aa = aa.round(roundprec)
            lis.append(aa)
        lis.append(Bis3[0])
        step1 = lis
    else:
        lis = []
        n = len(Bis3)
        for k in range(n-1, -1, -2):
            print('hi')
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk*Bk1
            aa = aa.round(roundprec)
            lis.append(aa)
        step1 = lis

    return step1


def compute_Bis_bysteps2(roundprec, Bis3):
    lis = []
    n = len(Bis3)
    if (n % 2 == 1):
        for k in range(n-1, 0, -2):
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk*Bk1
            aa = aa.round(roundprec)
            lis.append(aa)
        lis.append(Bis3[0])
        step1 = lis
    else:
        for k in range(n-1, -1, -2):
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk*Bk1
            aa = aa.round(roundprec)
            lis.append(aa)
        step1 = lis

    return step1


def compute_4Bis_bysteps(roundprec, Bis3):
    lis = []
    n = len(Bis3)
    if (n % 2 == 1):
        for k in range(n-1, 0, -4):
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)

            aa = Bk1*Bk
            aa = aa.round(roundprec)
            lis.append(aa)
        lis.append(Bis3[0])
        step1 = lis
    else:
        lis = []
        n = len(Bis3)
        for k in range(n-1, -4, -2):
            print([k, k-1])
            Bk = Bis3[k].round(roundprec)
            Bk1 = Bis3[k-1].round(roundprec)
            aa = Bk*Bk1
            aa = aa.round(roundprec)
            lis.append(aa)
            step1 = lis

    return step1


def compute_Bis_bysteps3(roundprec, Bis3):
    lis = []
    n = len(Bis3)
    for k in range(n-1, -1, -2):
        print([k, k-1])
        Bk = Bis3[k].round(roundprec)
        Bk1 = Bis3[k-1].round(roundprec)
        aa = Bk*Bk1
        aa = aa.round(roundprec)
        lis.append(aa)
    step1 = lis

    return step1


"""def compute_Bis_bysteps2 (q, n, roundprec,step1) :
    lis =[]
    n2=len(step1)
    for k in range (n2-1, -1, -2 ) :
        print(k)
        Bk=step1[k]
        Bk1=step1[k-1]
        aa=Bk*Bk1    
        aa=aa.round(roundprec)
        lis.append(aa)
        
    return lis

"""


def compute_z(factors, q, n, roundprec, listeoffullBi):
    """calcul de la constante de normalisation z 
    arguments : 
        ttfact  list of list: liste des facteurs gi du model 
    return :
        z real : constante de normalisation"""
    #listeoffullBi=compute_Bi(factors, q, n,roundprec)
    z = listeoffullBi[n-1]
    z = z.round(roundprec)
    lis = []
    lis.append(max(z.tt.r))

    for k in range(n-2, 0, -1):
        # print(listeoffullBi[k].shape)
        print(k)
        z = (listeoffullBi[k].round(roundprec))*z
        z = z.round(roundprec)
        print(z.erank)
        lis.append(max(z.tt.r))

    z = (listeoffullBi[0].round(roundprec))*z
    lis.append(max(z.tt.r))
    return z, lis

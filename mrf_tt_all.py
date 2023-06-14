#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 18:20:03 2022

@author: abouabdallah
"""

"""
code permettant de faire la procédure de novikov entière :
    Il y aurait les différentes méthodes pour le calcul des Ais et Bis 
    Calcul de la constante de normalisation
    Calcul des marginales 

"""
"""
Mrf
"""

""" 
Librairies 
"""
import numpy as np
#import math
import tt
import time

#my libs 
import bases
import model_gen
import texan_compute_marg

import  cores_melting
import compute_marginales
#import comparator_for_margins

#import symbolic_tt_mrf


# some bugs : function ais all mehods



""" some functions that may be used for Bis and Ais computation"""

# ======================================================================


def core_eye(r, n):
    """% Return 3-dimensional array filled with eye matrices,
    % which corresponds to non-essential dimension core."""
    dimensions = [r, n, r]  # the dimensions of the core
    # create core 3d matrix of zeros according to dimensions
    core = np.zeros(dimensions)
    for i in range(n):
        core[:, i, :] = np.eye(r)  # in every dimension compute eye matrix
    return (core)

# ======================================================================


def tt_zeros(vec):
    """
    initialise a tensor zeros tensor inspired bu tt.ones 
    """
    tenseur = tt.ones(vec)
    fur = tenseur.to_list(tenseur)
    for k in range(len(fur)):
        vic = fur[k].shape
        fur[k] = np.zeros(vic)
    nvtenseur = tt.ones(vec)
    nvtenseur = nvtenseur.from_list(fur)
    return(nvtenseur)


# fonctions pour le traitement tt de novikov


"""
begin cores :
    -svd 
    -add non essential dimensions 
    - compute model 

"""


def which_case(Lam) :
    cases=['assortative','hierarchical','ordered','dissortative','coreperiphery','hierarchical2']
    for k in range (5, -1, -1) : 
        Lam2=model_gen.thelam(cases[k])
        if (np.linalg.norm(Lam-Lam2)==0) :
            case=cases[k]

    return case

# ======================================================================


def tt_format2(psi):
    """ compute the tt format of a binarie factor 
    arguments :
        @psi matrix : matrix   
    return :
        @ listedescores list: liste of cores of the matrix"""

    psi11 = tt.tensor(psi)  # mise en tt format
    listedescores = tt.tensor.to_list(psi11)  # stocking the cores in a list

    return listedescores


# ======================================================================


def ajout_des_dimension_non_ess2(listedescores, nodesimp, n):
    """ initialize the non essentials dimensions : G_{k}^l \in \R^{n2\times r\times n2} (vectors matrix and vectors)

   arguments :
       @listedescores list : list of cores interfering in the link \ell  
       @nodesimp : index of nodes interfering in the link \ell  
       @ n  int , number of  individuals

   return :
       @ nvliste  liste : list of all cores of the matrix"""

    # \ell ={i,j}
    nvliste = []  # liste of all dimensions
    n0, r, n1 = listedescores[0].shape  # shape of G_i
    n1, r, n2 = listedescores[1].shape  # shape of G_j
    a = 0  # counter
    for k in range(n):  # for loop from 0 to n-1
        if k in nodesimp:  # if k =i or j
            # add the cores G_i or G_j in the list
            nvliste.append(listedescores[a])
            a += 1  # increase the counter by one so we start with G_i the G_j
        else:  # if k not equal i or j
            if k < nodesimp[0]:
                bla = core_eye(n0, r)  # create a vector (real mathematecly)
                # add it to list
                nvliste.append(bla)
            elif k > nodesimp[1]:
                bla = core_eye(n2, r)  # create a vector (real mathematecly)
                # add it to list
                nvliste.append(bla)
            else:  # if i<k<j
                bla = core_eye(n1, r)  # create à core_eye
                # add it to list
                nvliste.append(bla)

    return nvliste


# ======================================================================


def tt_factor_general(modele):
    """ compute the Gi and the non essentials dimensions
    arguments :
        @modele dict : mixed binaries unaries factors 
          @ keys list of integers : the links that compose our model ([i,j] for sbm)
          @values : psi that compose our model
        @wd  int: \in S vaut -1 ou +1 pour le modèle d'ising 
    return :
        @ ttfacteurs : liste de liste  remplie avec les différents G_i"""

    list_of_edges = list(modele.keys())  # [i,j] links of our model
    m = len(list_of_edges)  # number of links
    n = max(max(list_of_edges))+1  # number of individuals
    facs = bases.listoflist(m)  # empty list of list

    for l in range(m):
        psi = modele[list_of_edges[l]]  # extract psi
        listedescores = tt_format2(psi)  # compute svd
        facs[l] = ajout_des_dimension_non_ess2(
            listedescores, list_of_edges[l], n)  # add non essential dimensions

    return facs


""" generate the model from scratch and compute the factors """

"""generate the model from scratch """








def model_org(n, alpha, Lam, a):
    """generate the model from scratch and compute it's factors 
    arguments :
        @n : int, number of individuals 
        @Lam : matrix,  lambda 
        @alpha : list of parametre alpha 


    return :
        @ factors : liste of liste  of G_i"""

    q = len(alpha)  # number of classes
    # assign a classe to each individuals
    case= which_case(Lam) 
    Lam=model_gen.thelam(case)
    C,D= model_gen.mod (n,q, Lam, alpha)   
    SBM=model_gen.model(D, Lam, alpha)
    dica=SBM.multiply(a)
    dica=SBM.dica

    factors = tt_factor_general(dica)  # create factors

    return factors


# ======================================================================


def stock_all_Ai_q(factors, q, n):
    """ compute all the TT-matrices  Ai[z_i] for i =1, ..., n, z_i=0,q-1 
    arguments :
        factorCoreArray : list of list, each element factorCoreArray[i] give as the cores of tt matrice Ai[Zi]
        n : int, number of individuals 
        q= int, number of classes 
    return :
        listofAiq : list of Ai[z_i] for i =1, ..., n, Z_i=0, ..., q-1

    """

    listofAiq = []  # empty list
    for di in range(q):
        listofAiq.append(cores_melting.stock_all_Ai(factors, di, n))

    return listofAiq






def stock_all_Ai_q_compressed(factors, q, n, comp, numcomp):
    """ compute all the TT-matrices in melted form Ai[z_i] for i =1, ..., n, z_i=0,q-1 
    arguments :
        factorCoreArray : list of list, each element factorCoreArray[i] give as the cores of tt matrice Ai[Zi]
        n : int, number of individuals 
        q= int, number of classes 
        numcomp : number of new cores 
        comp : vectors, vectors with the number of old cores that constitue each new core
        
    return :
        listofAiq : list of Ai[z_i] in melted form for i =1, ..., n, Z_i=0, ..., q-1

    """

    listofAiq = []
    for di in range(q):
        listofAiq.append(cores_melting.stock_all_Ai_compressed(
            factors, di, n, comp, numcomp))

    return listofAiq

def stock_all_Ai_q_compressed_k(factors, q, n, k):
    """ compute all the TT-matrices in melted form Ai[z_i] for i =1, ..., n, z_i=0,q-1 
    arguments :
        factorCoreArray : list of list, each element factorCoreArray[i] give as the cores of tt matrice Ai[Zi]
        n : int, number of individuals 
        q= int, number of classes 
        k : nb of cors to melt
    return :
        listofAiq : list of Ai[z_i] in melted form for i =1, ..., n, Z_i=0, ..., q-1

    """

    numcomp, comp = cores_melting.how_comp(n, k)
    Ais =stock_all_Ai_q_compressed(factors, q, n, comp, numcomp)


    return Ais



def stock_all_Ai_q_compressed_k2(factors, q, n, k):
    """ compute all the TT-matrices in melted form Ai[z_i] for i =1, ..., n, z_i=0,q-1 
    arguments :
        factorCoreArray : list of list, each element factorCoreArray[i] give as the cores of tt matrice Ai[Zi]
        n : int, number of individuals 
        q= int, number of classes 
        k : nb of cors to melt
    return :
        listofAiq : list of Ai[z_i] in melted form for i =1, ..., n, Z_i=0, ..., q-1

    """

    numcomp, comp = cores_melting.how_comp2(n, k)
    Ais =stock_all_Ai_q_compressed(factors, q, n, comp, numcomp)


    return Ais



""" non uniform compression """
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


"""def build_symbolic_ttfactors(n, edges):

    m = edges.shape[0]
    facs = [[] for i in range(m)]
    for l in range(m):
        cfacs = ['empty' for i in range(n)]
        i = edges[l, 0]
        j = edges[l, 1]
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

"""


# ======================================================================


def symbolic_product(Bis3):
    lis = []
    n = len(Bis3)
    for k in range(n-1, -1, -2):
        print([k, k-1])
        Bk = Bis3[k].tt.r
        Bk1 = Bis3[k-1].tt.r
        aa = bases.multi(Bk, Bk1)
        lis.append(aa)
    shapes = lis
    return shapes



def compute_Ai_per_method(dic, method, k, nb_paquets, q, n) : 
    """compute Ai per method :
        It provide 3 methods :
            Ais without compression, with uniform compression and uniform volume comression.
        arguments :
            dic : mixed factors 
            q : number of classes
            n : number of individuals 
        return : 
            Ai : List of lists of tt matrix Ai[Z_i]
    """
    factors = tt_factor_general(dic)
    if (method == 0) :
        Ai = stock_all_Ai_q(factors, q, n)
    if (method == 1) :
        Ai = stock_all_Ai_q_compressed_k2(factors, q, n, k)
    if (method == 2) :
        nmax=5
        bgroupes, numcorps =cores_melting.how_comp_bis(dic, nb_paquets, nmax , q)
        Ai = stock_all_Ai_q_compressed(factors, q, n, bgroupes, numcorps )


        
    return Ai


def compute_Ai_per_method2(dic, q, n) : 
    """compute Ai per method :
        It provide 3 methods :
            Ais without compression, with uniform compression and uniform volume comression.
        arguments :
            dic : mixed factors 
            q : number of classes
            n : number of individuals 
        return : 
            Ai : List of lists of tt matrix Ai[Z_i]
    """
    method=int(input("0 : sans compression, 1 : uniforme, 2 : volumes uniformes"))
    methodes=["Without fusion","uniforme nb of cores fusion","uniform volumes fusion"]
    factors = tt_factor_general(dic)
    print(methodes[method])
    if (method == 0) :
        Ai = stock_all_Ai_q(factors, q, n)
    if (method == 1) :
        k=int(input("Nombre de cores compresses"))
        Ai = stock_all_Ai_q_compressed_k(factors, q, n, k)
    if (method == 2) :
        nmax=int(input("Nombre max de cores"))
        nb_paquets=int(input("Bombre de paquets"))
        bgroupes, numcorps =cores_melting.how_comp_bis(dic, nb_paquets, nmax , q)
        Ai = stock_all_Ai_q_compressed(factors, q, n, bgroupes, numcorps )


        
    return Ai



def compute_Ais_methods(dic, method, k, nb_paquets, q, n):
    n, m, nodes_edges, edges = nds_edges(dic)
    factors = tt_factor_general(dic)

    if (method == 1):
        Ais = stock_all_Ai_q(factors, q, n)
    elif (method == 2):
        m = bases.get_m(n)
        numcomp, comp = cores_melting.how_comp(n, k)
        Ais = stock_all_Ai_q_compressed(factors, q, n, comp, numcomp)
    else:
        sfactors = symbolic_factors(dic)
        n, m, nodes_edges, edges = nds_edges(dic)
        load, comp, sA = cores_melting.comp_load(sfactors, dic)
        max_ind = np.argmax(comp[:, 1])
        load2 = cores_melting.volume_cores_A_i(max_ind, sA, m, q)
        paquets_load2v, paquets2v = cores_melting.init_groupage_Ai(load2, nb_paquets)
        groupes = cores_melting.groupage_Ai_main(load2, nb_paquets)
        nmax = 7
        bgroupes = cores_melting.groupage_Ai_all(groupes, nmax)
        numcorps = len(bgroupes)
        Ais = stock_all_Ai_q_compressed(factors, q, n, bgroupes, numcorps)
    return Ais


# calcul des bis


def compute_Bi_f(factors, q, n, roundprec):  # normale
    Bi = []
    listofAiq = stock_all_Ai_q(factors, q, n)
    i = 0
    for i in range(n):
        rowRankArray = listofAiq[0][i].n
        colRankArray = listofAiq[0][i].m
        a = tt_zeros(bases.multi(rowRankArray, colRankArray))
        a = listofAiq[0][i]
#        a=a.round(roundprec)
        for di in range(1, q):
            # print(di)
            a = a+listofAiq[di][i]
            # .round(roundprec)
#            a=a.round(roundprec)
        Bi.append(a)
    return Bi


def compute_Bi_f_compressed(factors, q, n, roundprec, comp, numcomp):  # avec compression
    Bi = []
    listofAiq = stock_all_Ai_q_compressed(factors, q, n, comp, numcomp)
    i = 0
    for i in range(n):
        print(i)
        rowRankArray = listofAiq[0][i].n
        colRankArray = listofAiq[0][i].m
        a = tt_zeros(bases.multi(rowRankArray, colRankArray))
        a = listofAiq[0][i]
#        a=a.round(roundprec)
        for di in range(1, q):
            # print(di)
            a = a+listofAiq[di][i]
            # .round(roundprec)
#            a=a.round(roundprec)
        Bi.append(a)
    return Bi


def compute_Bi_all(listofAiq, q, n, roundprec):  # normale
    Bi = []

    for i in range(n):
        rowRankArray = listofAiq[0][i].n
        colRankArray = listofAiq[0][i].m
        a = tt_zeros(bases.multi(rowRankArray, colRankArray))
        a = listofAiq[0][i]
#        a=a.round(roundprec)
        for di in range(1, q):
            # print(di)
            a = a+listofAiq[di][i]
            # .round(roundprec)
#            a=a.round(roundprec)
        Bi.append(a)
    return Bi


def loads(Bis3):
    load = []
    for k in range(len(Bis3)):
        Bi = Bis3[k]
        ll = Bi.to_list(Bi)
        a = 0
        for l in range(len(ll)):
            a += np.log(ll[l].size)
        load.append(a)

    return load


def loadscores(Bi):
    load = []
    ll = Bi.to_list(Bi)
    for l in range(len(ll)):
        a = np.log(ll[l].size)
        load.append(a)

    return load






# compute z and to float


def compute_z(listeoffullBi, n, roundprec):
    """calcul de la constante de normalisation z 
    arguments : 
        ttfact  list of list: liste des facteurs gi du model 
    return :
        z real : constante de normalisation"""

    z = listeoffullBi[n-1]
    z = z.round(roundprec)
    for k in range(n-2, 0, -1):
        # print(listeoffullBi[k].shape)
        # print(k)
        z = listeoffullBi[k]*z
        z = z.round(roundprec)

    z = listeoffullBi[0]*z
    return z


def compute_acc_smalln(roundprec, Bis2, ff):

    listeoffullBi = Bis2
    n = len(listeoffullBi)
    z = listeoffullBi[n-1]
    z = z.round(roundprec)
    for k in range(n-2, 11, -1):
        # print(k)
        z = (listeoffullBi[k].round(roundprec))*z
        z = z.round(roundprec)

    z2 = (listeoffullBi[11].round(roundprec))
    z2 = z2.round(roundprec)
    for k in range(10, -1, -1):
        # print(k)
        z2 = (listeoffullBi[k].round(roundprec))*z2
        z2 = z2.round(roundprec)

    z2 = ff*z2
    z2 = z2.round(roundprec)
    aaaaa = z2*z

    return (aaaaa)


def to_float(Z):
    lis = Z.to_list(Z)
    lis2 = 1
    for k in range(len(lis)):
        lis[k] = np.reshape(lis[k], [lis[k].shape[0], lis[k].shape[3]])
        lis2 = np.dot(lis2, lis[k])
    lis2 = np.float(lis2)
    return (lis2)


################################################################################
""" calcul des marginales 
"""


def compute_acc_smalln_onemar_onestate(roundprec, Bis2, ff, Ais2, di, dj, i, j):

    listeoffullBi = Bis2
#    n=len(listeoffullBi)
    listeoffullBi[i] = Ais2[di][i]
    listeoffullBi[j] = Ais2[dj][j]
    lal = compute_acc_smalln(roundprec, listeoffullBi, ff)
    mar = to_float(lal)

    return (mar)


def compute_acc_smalln_onemar_allstate(roundprec, Bis2, q, ff, Ais2, i, j):

    aaa = np.ones((q, q))
    for di in range(q):
        for dj in range(q):
            aaaaa = compute_acc_smalln_onemar_onestate(
                roundprec, Bis2, ff, Ais2, di, dj, i, j)
            aaa[di, dj] = np.abs(aaaaa)

    return (aaa)


def comp_acc_smalln_allstate(roundprec, Bis2, ff, Ais2, q, n, w):
    thedic = {}
    lis = []

    print("it begins")
    for i in range(n):
        for j in range(i+1, n):
            print([i, j])
            aaa = np.ones((q, q))
            for di in range(q):
                for dj in range(q):
                    aaaaa = compute_acc_smalln_onemar_onestate(
                        roundprec, Bis2, ff, Ais2, di, dj, i, j)
                    aaa[di, dj] = aaaaa
            thedic[i, j] = aaa
            lis.append(aaa)
            print(aaa)
#            directory=str(i)+str(j)+'.txt'
#            np.savetxt(directory,aaa,fmt='%.2f')

    return (thedic, lis)


def comp_acc_smalln_allstate_unaries(roundprec, Bis2, ff, Ais2, q, n, w):
    thedic = {}

    print("it begins")
    for i in range(n):
        print(i)
        aaa = np.ones(q)
        for di in range(q):
            aaaaa = compute_acc_smalln_onemar_onestate(
                roundprec, Bis2, ff, Ais2, di, di, i, i)
            aaa[di] = aaaaa
            thedic[i, i] = aaa
            print(aaa)

    return thedic


# fixed rmx


def compute_z_timed(listeoffullBi, roundprec, n):
    start_time = time.time()
    W = 1
    rmax = []
    for k in range(n-1, -1, -1):

        #    print(k)
        Bi = listeoffullBi[k].round(roundprec)
        W = Bi*W
        W = W.round(roundprec)
        rmax.append(max(W.tt.r))
    end_time = time.time()
    rm = max(rmax)
    ti = end_time-start_time
    Z = to_float(W)
    return [ti, Z], rm


def compute_W_simple(Bis, roundprec, rmax, ff) :
    n=len(Bis)
    
    W=1
    for k in range (n-1, -1, -1) :
        Bi=Bis[k]*ff
        Bi=Bis[k].round(roundprec, 3)
        W=Bi*W
        W=W.round(roundprec, rmax)
    Z=to_float(W)

    return Z



def compute_z_timed_r(listeoffullBi, roundprec, n, rmax):
    start_time = time.time()
    W = 1
    for k in range(n-1, -1, -1):

        #    print(k)
        Bi = listeoffullBi[k].round(roundprec, rmax)
        W = Bi*W
        W = W.round(roundprec, rmax)
    end_time = time.time()
    ti = end_time-start_time
    Z = to_float(W)
    return ti, Z


def Bis_to_multiply(n, kmax):
    ls = []
    liscormult = []
    a = 0
    for k in range(n-1, -1, -1):
        if (a < kmax - 1):
            a += 1
            ls.append(k)
        else:
            ls.append(k)
            a = 0
            liscormult.append(ls)
            ls = []

    return liscormult


def compute_bis_by_step(kmax, listeoffullBi, rmax, ff, roundprec):
    n = len(listeoffullBi)
    liscormult = Bis_to_multiply(n, kmax)


#    listepfwi=[]
    W = 1
#    a=0
    mm = len(liscormult)
    liste = []
    for k in range(mm):
        W = 1
        for i in range(len(liscormult[k])):
            # print(liscormult[k][i])
            Bi = listeoffullBi[liscormult[k][i]] * ff
            W = Bi*W
            W = W.round(roundprec, rmax)
        liste.append(W)
    W = 1
    for k in range(len(liste)):
        W = liste[k]*W
        W = W.round(roundprec, 2*rmax)  # or rmax^2
    Zx = to_float(W)

    return Zx


def compute_bis_milieu(listeoffullBi, rmax, roundprec, ff):
    Bis3 = listeoffullBi
    n = len(listeoffullBi)

    wtest = Bis3[int(n/2)]
    wtest = wtest*ff
    for i in range(1, int(n/2)):
        print([int(n/2) - i, int(n/2), int(n/2) + i])
        wtest = wtest*Bis3[int(n/2) + i]
        wtest = wtest*ff
        wtest = wtest.round(roundprec, rmax)
        wtest = Bis3[int(n/2) - i]*wtest
        wtest = wtest*ff
        wtest = wtest.round(roundprec, rmax)
    wtest = Bis3[0]*wtest
    wtest = wtest*ff
    Zx2 = to_float(wtest)

    return Zx2


def compute_normes(listeoffullBi):
    Bis3 = listeoffullBi
    n = len(listeoffullBi)
    listeofnormes = []
    i = 0
    Bi = Bis3[i]

    for i in range(n):
        Bi = Bis3[i]
        listeofnormes.append(Bi.norm())

    return listeofnormes


def normalize_Bis(listeoffullBi):
    listeofnormes = compute_normes(listeoffullBi)
    n = len(listeofnormes)
    normBis = []
    for k in range(n):
        normBis.append((listeoffullBi[k])*(1./listeofnormes[k]))

    return normBis

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 13:59:31 2022

@author: abouabdallah
"""

import numpy as np
import math 

""" generate the model """

# ----------------------------------------------------------------------


def assig_multinom_classes(n, alpha):
    """ finding the classes and building Z

    arguments :
        n : int nombre d'individuals
        alpha : liste of probabilities of bellonging to classes
        outputs :
            C : vector \in R^{n} classes of individuals
            Z : matrix \in R^{n,b} classes of individuals


    """
    Z = np.random.multinomial(
        1, alpha, n)  # generattng the classes of individuals according to a multinomial law
    B = len(alpha)
    C = [0]*n  # vector of shape n
    for i in range(n):
        for k in range(B):
            if (Z[i, k] == 1):
                C[i] = k
    C = [int(k) for k in C]

    return (C, Z)


# ----------------------------------------------------------------------

def thelam(case):
    """ Creating the lambdas matrices

      inputs :
          case: char, the case of lambdas
          outputs :
              Lam : matrix \in \R^{QxQ} the Lambdas matrix
    """

    Lam = 2*np.ones((3, 3))
    if (case == 'assortative'):
        Lam[0, 0], Lam[1, 1], Lam[2, 2] = 2, 3, 4
        Lam[1, 0], Lam[2, 0], Lam[2, 1], Lam[0, 1], Lam[0,
                                                        2], Lam[1, 2] = 10, 10, 10, 10, 10, 10
    if (case == 'hierarchical'):
        Lam[0, 0], Lam[1, 1], Lam[2, 2] = 12, 6, 3
        Lam[1, 0], Lam[2, 0], Lam[2, 1], Lam[0,
                                             1], Lam[0, 2], Lam[1, 2] = 3, 1, 3, 3, 1, 3
    if (case == 'ordered'):
        Lam[0, 0], Lam[1, 1], Lam[2, 2] = 10, 10, 10
        Lam[1, 0], Lam[2, 0], Lam[2, 1], Lam[0,
                                             1], Lam[0, 2], Lam[1, 2] = 6, 2, 6, 6, 2, 6
    if (case == 'dissortative'):
        Lam[0, 0], Lam[1, 1], Lam[2, 2] = 8, 10, 9
        Lam[1, 0], Lam[2, 0], Lam[2, 1], Lam[0,
                                             1], Lam[0, 2], Lam[1, 2] = 3, 3, 3, 3, 3, 3
    if (case == 'coreperiphery'):
        Lam[0, 0], Lam[1, 1], Lam[2, 2] = 2, 4, 10
        Lam[1, 0], Lam[2, 0], Lam[2, 1], Lam[0,
                                             1], Lam[0, 2], Lam[1, 2] = 5, 8, 6, 5, 8, 6
    if (case == 'hierarchical2'):
        Lam[0, 0], Lam[1, 1], Lam[2, 2] = 3, 6, 10
        Lam[1, 0], Lam[2, 0], Lam[2, 1], Lam[0, 1], Lam[0,
                                                        2], Lam[1, 2] = 10, 15, 15, 10, 15, 15

    return Lam

def is_in_case(lam, cases) :
    for case in cases : 
        if (np.linalg.norm(lam-thelam(case))==0) :
            return True
    return False

# ----------------------------------------------------------------------

def mat_lam_aleatoire(B, dmax, dmin):
    """ generates random matrix Lambda
    inputs :
        B : int number of classes
        dmax : int distance max
        dmin : int distance min
    return :
        lam : distances betewen classes (randomly)
    """
    Lam = np.zeros((B, B))
    for i in range(B-1):
        for j in range(i, B):
            a = np.random.randint(dmin, dmax)
            Lam[i, j] = a
            if j != i:
                Lam[j, i] = Lam[i, j]

    return (Lam)

# ----------------------------------------------------------------------


def mat_distance(n, Lam, C):
    """ generate distance mat
    input :
        n : int number of individuals
        lam : distances between classes
        C  : classes des individuals
    output :
        D : matrix, distances matrix
    """
    D = np.zeros((n, n))
    for i in range(n-1):
        k = C[i]
        # k=int(k)
        # D[i,i]=np.random.poisson(lam=lam[k,k])
        for j in range(i+1, n):
            l = int(C[j])
            D[i, j] = np.random.poisson(lam=Lam[k, l])
            D[j, i] = D[i, j]
    #
    return (D)

# ======================================================================


"""initialisations used in ttem"""


def init_alpha_random(q):
    """ initialisation of \alpha randomly (used for TT-em algorithme)
    input :
        q : number of classes 
    output :
        alpha ; list of probabilities of appending to classes (initialisation)
    """
    alpha1 = np.random.random(q)
    m = np.sum(alpha1)
    div = 1/m
    alpha = div*alpha1
    return alpha

# ======================================================================


def init_lam(q, D):
    """ initialisation of \lambda randomly (used for TT-em algorithme), computed as the mean of Distance matrix
    input :
        q : number of classes 
        D : distance matrix 
    output :
        Lam ; distances betewen classes (random initialisation)
    """
    a = np.mean(D)
    Lam = np.zeros((q, q))
    for k in range(q):
        for l in range(q):
            Lam[k, l] = a

    return (Lam)



# ======================================================================


def which_case(Lam) :
    cases=['assortative','hierarchical','ordered','dissortative','coreperiphery','hierarchical2']
    for k in range (5, -1, -1) : 
        Lam2=thelam(cases[k])
        if (np.linalg.norm(Lam-Lam2)==0) :
            case=cases[k]

    return case

# ======================================================================


def mod (n,B, Lam, alpha) :
    C,Z=assig_multinom_classes(n,alpha)
    D=mat_distance(n,Lam, C) 
    return C,D
# ======================================================================


"""the computation of the factors : 
    -unaries factors  \phi_{i}
    -binaries factors psi_{ij}
    -mixed unaries binaries 
    
    """


# ======================================================================
def fact_value(Lam, d, dfac, b, bp):
    """ Compute  psi_{ij}[b,b']
    input : 
        Lam : matrix of paramètre lambda 
        d: int, distance 
        dfac: int, factoriel of a distance
        b: int, index
        bp: int, index
    output :
        y : real, psi_{ij}[b,b']
    """
    lam = Lam[b, bp]
    y = lam**d
    y = y/dfac
    y = y*math.exp(-lam)
    return y


# ======================================================================

def binary_factors(n, B, D, Lam):
    """ Compute the binaries factors
    input : 
        Lam : matrix of parametre lambda 
        D: int matrix, matrix of distance 
        n: int, number of individuals
        B: int, number of classes
    output :
        dicb : dict, dicionnary of binaries factors
             -Keys [i,j], i \in 1, ...n-1 and j =i+1, ...n
             -Values : binaries factors 
    """
    dicb = {}  # empty dictionnary of binaries
    for i in range(n-1):
        for j in range(i+1, n):
            d = int(D[i, j])
            dfac = math.factorial(d)  # compute D_{i,j}!
            psi = np.zeros((B, B))  # initializ a matrix of binaries
            for b in range(B):
                for bp in range(B):

                    psi[b, bp] = fact_value(Lam, int(d), dfac, int(
                        b), int(bp))  # Copute  psi_{ij}[b,b']
            dicb[(i, j)] = psi

    return dicb


# ======================================================================


def unary_factors(n, B, D, alpha):
    """ Compute unary factors 
    input : 
        alpha : list of parametre alpha 
        D: int matrix, matrix of distance 
        n: int, number of individuals
        B: int, number of classes
    output :
        dicu : dict, dicionnary of unaries factors
             -Keys [i], i \in 1, ...n
             -Values : unary factors 
    """
    dicu = {}  # empty dictionnary of unaries
    for i in range(n-1):
        phi = np.zeros(B)  # initialisaion of unary factor as vector
        for b in range(B):
            phi[b] = alpha[b]

        dicu[(i)] = phi

    return dicu


# ======================================================================


"""def updated_binary_factors(dicu, dicb):
    compute mixed binaries unaries 
    input : 
        dicu : dict, dicionnary of unaries 
        dicb : dict,dicionnary of  binaries
    output :
        dic : dict, dicionnary of mixed binaries unaries 
    n = max(max(dicb.keys()))
    B = dicb[0, 1].shape[0]
    dic = dicb
    for i in range(n-2):
        psi = dicb[i, i+1]
        phi = dicu[i]
        for b in range(B):
            for bp in range(B):
                # \psi'_{i,i+1}(Z_i,Z_{i+1}) = \psi_{i,i+1}(Z_i,Z_{i+1}) \phi_i(Z_i) if i=j
                psi[b, bp] = psi[b, bp]*phi[b]
                # \psi'_{i,j}(Z_i,Z_{j}) = \psi_{i,j}(Z_i,Z_{j})  else

        dic[i, i+1] = psi

    return dic
"""





def updated_binary_factors2(dicu, dicb):
    """ compute mixed binaries unaries 
    input : 
        dicu : dict, dicionnary of unaries 
        dicb : dict,dicionnary of  binaries
    output :
        dic : dict, dicionnary of mixed binaries unaries 
    """
    n = max(max(dicb.keys()))
    B = dicb[0, 1].shape[0]
    dic = dicb
    dic2=dict()
    for i in range (n-1) :
        psi = dicb[i, i+1]
        phi = dicu[i]
        psip=np.zeros((B,B))
        #print(dicb[i, i+1])
        #print(dicb[i, i+1] *10)
        #print(dic[i, i+1])
        for b in range(B):
            for bp in range(B):
                dic[i, i+1][b, bp] = psi[b, bp]*phi[b] # \psi'_{i,i+1}(Z_i,Z_{i+1}) = \psi_{i,i+1}(Z_i,Z_{i+1}) \phi_i(Z_i) if i=j
                psip[b, bp]=psi[b, bp]*phi[b]
        dic2[i,i+1]= psip
        #print(dicb[i, i+1])

        #print((i, i+1,dic[i, i+1]-dicb[i, i+1]))
    for i in range (n) :
        for j in range (i+2,n+1) :
            #print(i,j)
            dic2[i,j]=dicb[i,j]
            


    psi=dicb[n-1, n]    
    phi = dicu[n-2]
    phi2 = dicu[n-1]
    psip=np.zeros((B,B))
    for b in range(B):
        for bp in range(B):
            #print(psi[b, bp]*phi[b]*phi2[b])
            psip[b, bp] = psi[b, bp]*phi[b]*phi2[b]
    dic2[n-1, n]=psip


    return dic2



def fac_a(dic, a):
    """ Compute mixed factors * a
    usefull for hight n 
    input : 
        dic : dict, dicionnary of mixed binaries unaries 
    output :
        dic : dict, dicionnary of mixed binaries unaries * a
    """
    for k in range(len(list(dic.keys()))):
        dic[list(dic.keys())[k]] = dic[list(dic.keys())[k]]*a
    return dic







class model:
    """Define SBM model
    input : 
        D : array, distance matrix
        Lam : array, distances matrix 
        alpha : list, proba bellonging to classes
    output :
        the model :
            D : array, distance matrix
            Lam : array, distances matrix 
            alpha : list, proba bellonging to classes
            n : int, number of individuals
            q : int, number of classes 
            dicu : dict, unarie factors
            dicb : dict, binarie factors 
            
    """
    
    def __init__(self, *args ):# n, case, alpha
        D = args[0] 
        self.D=D
        n=D.shape[0]
        self.n = n
        self.m = n*(n-1)/2
        Lam = args[1]
        self.Lam= Lam
        cases=['assortative', 'hierarchical','ordered', 'dissortative','coreperiphery','hierarchical2']

        if (is_in_case(Lam, cases)) :
            self.case = which_case(Lam)
        self.case = "No case"
        alpha=args[2]
        self.alpha=alpha
        q=len(alpha)
        self.q=q
        dicu=unary_factors(n,q,D, alpha)
        dicb=binary_factors (n,q, D, Lam) 
        self.dicu=dicu
        self.dicb=dicb
        dic=updated_binary_factors2(dicu, dicb)        
        self.dic=dic
    def multiply(self, a):
        dic2=self.dic
        for k in range(len(list(dic2.keys()))):
            dic2[list(dic2.keys())[k]] = dic2[list(dic2.keys())[k]]*a
        self.dica=dic2
        

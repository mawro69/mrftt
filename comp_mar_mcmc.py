#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 12:19:46 2022

@author: abouabdallah
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 25 10:25:43 2022

@author: abouabdallah
"""

import numpy as np
import enamantem as ena
import math

# it begins


def unary_factors(n, Q, D, alpha):
    """ Compute unaries marginales 
    input : 
        alpha :  SBM model parameters alpha
        D: int matrix n*n, distance matrix
        n: int, number of individual
        Q: int, number of classes
    output :
        dicu : dict, dicionnairy of unaries factors 
    """
    dicu = {}
    for i in range(n):
        phi = np.zeros(Q)
        for b in range(Q):
            phi[b] = alpha[b]

        dicu[(i)] = phi

    return dicu


def dicu_tomat(dicu, n, Q):
    """ Transform unaries factors dictionnarie to matrix 
    input : 
        dicu: dict, dicionnairy of unaries factors 
        n: int, number of individual
        Q: int, number of classes

    output :
        unaries : np array, unaries factors in array form
    """
    unaries = np.ones((n, Q))
    for i in range(n):
        for b in range(Q):
            unaries[i, b] = dicu[i][b]

    return unaries


def all_fac(dicb):
    """ Compute binaries factors in another format
    input : 
        dicb: dict, dicionnairy of binaries factors 

    output :
        binaries : dict, binaries factors 
    """

    binaries = dict()
    n = max(max(list(dicb.keys())))+1
    for i in range(n-1):
        for j in range(i+1, n):
            binaries[(i, j)] = dicb[(i, j)]
            binaries[(j, i)] = np.transpose(dicb[(i, j)])

    return binaries

# we computed the binaries and unaries factors


# gibbs sampler initialisation
def assig_multinom_classes(n, alpha):
    """ finding the classes and building Z

    inputs :
            n : int number of individuals
            alpha : list, list of probas of belonging to classes
        outputs :
            C : vector \in R^{n}  classes of individuals
            Z : matrix \in R^{n,b}  classes of individuals


    """
    Z = np.random.multinomial(
        1, alpha, n)  # multinomial labels for individuals
    B = len(alpha)  # lenghth of alpha
    C = [0]*n  # vector full of zeros
    for i in range(n):  # loop on infividuals
        for k in range(B):  # loop on classess
            if (Z[i, k] == 1):
                C[i] = k
    C = [int(k) for k in C]  # liste of individuals classes

    return (C, Z)  # we return the a list of classes and the matrix of classes





#compute the conditionnal probab

def compute_proba(i, b, dicu, dicb, C):
    """ finding the classes and building Z
    inputs :
           i : integer, the number of the individual
           b : integer, the state of the individual 
           dicu :  array, unary factors
           dicb :  dic, binary factors 

    outputs :
        proba : float,  P_{\theta}(Z_i= z_i|\{Z_j= z_j\}_{j \neq i}, D)


    Noteaf

    1) écrire : proba : float,  P_{\theta}(Z_i= b|\{Z_j= C[j]\}_{j \neq i}, D)
    car z_i=b et z_j = C[j]

    2) la façon de retrouver le nbre d'individus est un peu compliquée ... je le mettrais en argument ... 
    """

    n = max(max(list(dicb.keys())))+1  # number of individuals
    proba = dicu[i, b]  # unary factor
    for j in range(n):  # loop on the individuals
        if (j != i):  # if i is not j  from product formula \prod_{j \neq i}
            proba = proba*dicb[i, j][b, C[j]]  # product formula
    return proba





#gibbs step

def echantilloneur_gibbs_step_cor(n, q, dicu, dicb, Cold):
    """ compute one stepps of Gibbs sampler

      inputs :
              n : int, number of individuals
              q : integer, the number of classes
              dicu :  array, unary factors
              dicb :  dic, binary factors 
              Cinit : classes computed in the last step


          outputs :
              C : vector \in R^{n}  classes of individuals
              Z : matrix \in R^{n,b}  classes of individuals

      Noteaf : 

      1) Z est une matrice R^{n,q}
      2) je mettrais un Cold et un Cnew, plutôt qu'un Cinit pour les deux ...
      3) le calcul de init me semble bien compliqué ...
      4) est-ce que Z = np.random.multinomial(1, pcond, 1)  n'est pas déjà un vecteur disjonctif complet, partout des 0 sauf un 1 en une position ?
      5) pas compris le role de vectolist() ... 

      """
    Cnew=Cold
    dicu1=np.zeros((n,q))
    for i in range(n):  # loop on the individuals
        pcond = np.zeros(q)
        for b in range(q):  # loop on the classes
            # we compute the proba=proba*dicb[i,j][b,C[j]]
            pcond[b] = compute_proba(i, b, dicu, dicb, Cnew)
        pcond = pcond/np.sum(pcond)  # normalization
        
        dicu1[i,]=pcond #ajoutée
        #print(pcond)

        # multinomial assignation of individual in classes
        Z = np.random.multinomial(1, pcond, 1)
        #print(Z)

        Z_pl = Z.tolist()[0]
        #print([i, Z_pl.index(1)])
        Cnew[i] = Z_pl.index(1)  # the max will be the assigned class
#    print(Cnew)

    return Cnew, dicu1




#gibbs one step

def echantilloneur_gibbs_allsteps(n, q, dicu, dicb, Niter):
    """ compute all stepps of Gibbs sampler

      inputs :
              n : int, number of individuals
              q : integer, the number of classes
              dicu :  array, unary factors
              dicb :  dic, binary factors 
              Niter : number of iterations

          outputs :
              C : vector \in R^{n}  classes of individuals
              Z : matrix \in R^{n,b}  classes of individuals

      noteaf
      1) pourquoi une taille à (q) [un tuple] et non ps simplement q [un entier] ? 
      """
    #np.random.dirichlet(np.ones(3),size=1)

    #alpha = np.random.randint(1, 10, (q))
    #alpha = alpha/alpha.sum()  # random init of alpha
    alpha = np.random.dirichlet(np.ones(q),size=1)
    alpha=alpha[0].tolist()
    
    
    Cinit, Zinit = ena.assig_multinom_classes(n, alpha)  # initialization of C
    #print(Cinit) #we print it

    dictot = q*np.ones((Niter, n))
    dicu1 = dicu
    for l in range(Niter):
        # we compute an iteratuin of gibbs sampler
        Cinit, dicu1 = echantilloneur_gibbs_step_cor(n, q, dicu1, dicb, Cinit)
        #print(Cinit)

        dictot[l] = Cinit  # we  store the result for the tests
    #print(Cinit)
    return Cinit, dicu



def echantilloneur_gibbs_allsteps_stockage(n, q, dicu, dicb, Niter):
    """ compute all stepps of Gibbs sampler

      inputs :
              n : int, number of individuals
              q : integer, the number of classes
              dicu :  array, unary factors
              dicb :  dic, binary factors 
              Niter : number of iterations

          outputs :
              C : vector \in R^{n}  classes of individuals
              Z : matrix \in R^{n,b}  classes of individuals

      noteaf
      1) pourquoi une taille à (q) [un tuple] et non ps simplement q [un entier] ? 
      """
    #np.random.dirichlet(np.ones(3),size=1)

    #alpha = np.random.randint(1, 10, (q))
    #alpha = alpha/alpha.sum()  # random init of alpha
    alpha = np.random.dirichlet(np.ones(q),size=1)
    alpha=alpha[0].tolist()
    
    
    Cinit, Zinit = ena.assig_multinom_classes(n, alpha)  # initialization of C
    #print(Cinit) #we print it

    dictot = q*np.ones((Niter, n))
    dicu1 = dicu
    for l in range(Niter):
        # we compute an iteratuin of gibbs sampler
        Cinit, dicu1 = echantilloneur_gibbs_step_cor(n, q, dicu1, dicb, Cinit)
        #print(Cinit)

        dictot[l] = Cinit  # we  store the result for the tests
    #print(Cinit)
    return dictot

#Cinit, dicu1=echantilloneur_gibbs_allsteps(n, q, dicu, dicb, Niter)
#C


def echantillons_token(n, q, dicu, dicb, Niter, firstone, latence) :
    dicmcmc=echantilloneur_gibbs_allsteps_stockage(n, q, dicu, dicb, Niter)
    size = int((Niter -firstone)/latence)
    echantillons = np.zeros((size,dicmcmc.shape[1]))
    for k in range (firstone , Niter, latence ) :
        echantillons[int((k-firstone)/latence), :] =dicmcmc[k, :]
    return echantillons


def Lechant(n,q,dicu, dicb, Niter, L) :
    """ compute L Gibbs sampler

      inputs :
              n : int, number of individuals
              q : integer, the number of classes
              dicu :  array, unary factors
              dicb :  dic, binary factors 
              Niter : int, number of iterations
              L : int, number of sampler

          outputs :
              dictot : array,  classes of individuals computed by L sampler.
      """

    dictot=q*np.ones((L, n))
    for l in range (L) :
        print(["initial classes for the sampler", l])
        Cinit, dicu= echantilloneur_gibbs_allsteps(n,q,dicu, dicb, Niter)
        #print(["classes obtained for the sampler", l])
        #print(Cinit)
        dictot[l]= Cinit
    
    return dictot


def margunaires(dictot) :
    q= int(np.max(dictot))+1
    L,n=dictot.shape
    unaries=np.zeros((n,q))
    for i in range (n) :
        for l in range (L) :
            for b in range (q) :
                if (dictot[l,i] == b) :
                    unaries[i,b]+=1
    unaries=unaries/L
    
    return unaries


def margbinaires(dictot) :
    q= int(np.max(dictot))+1
    L,n=dictot.shape
    binaries=dict()
    for i in range (n-1) :
        for j in range (i+1,n) :
            binaries[i,j] =np.zeros((q,q))
            for l in range (L) :

                for b in range (q) :
                    for bp in range (q) :
                        if (dictot[l,i] == b and dictot[l,j] == bp) :
                            
                            binaries[i,j][b,bp]+= 1
            binaries[i,j]=binaries[i,j]/L
                            
                    
    
    return binaries 



def mat_distance(n,Lam, C) :
     """ generate distance mat
     input :
         n : int nombre d'individus
         lam : distances entre les classes
         C  : classes des individus
     output :
         D : matrix, matrice des distances
     """
     D = np.zeros((n,n))
     for i in range (n-1):
         k = C[i]
         #k=int(k)
         #D[i,i]=np.random.poisson(lam=lam[k,k])
         for j in range (i+1,n) :
             l = int(C[j])
             D[i,j] = np.random.poisson(lam=Lam[k,l])
             D[j,i] = D[i,j]
     #
     return  (D)

# ======================================================================




"""factors """


def app_classes(matclases, n,q ) :
    classes=[]
    for k in range (n) :
        classes.append(list(matclases[k,]).index(max(matclases[k,])))
    return classes


def init_alpha_random(q) :
    alpha1 =np.random.random(q)
    m=np.sum(alpha1)
    div=1/m
    alpha=div*alpha1
    return alpha 

def init_lam(q, D) :
    a= np.mean(D)
    Lam=np.zeros((q,q))
    for k in range (q) :
        for l in range (q) :
            Lam[k,l]=a
            
            
    return (Lam)



# ======================================================================
def fact_value (Lam, d,dfac, b, bp) :
    """ calcul psi_{ij}[b,b']
    input : 
        Lam : matrice paramètre lambda d'un modèle SBM
        d: int, une distance 
        dfac: int, la factorielle d'une distance
        b: int, indice
        bp: int, indice
    output :
        y : real psi_{ij}[b,b']
    """
    lam=Lam[b,bp]
    y=lam**d
    y=y/dfac
    y=y*math.exp(-lam)
    return y


# ======================================================================

def binary_factors (n,B, D, Lam) :
    """ calcul des marginales bianires 
    input : 
        Lam : matrice paramètre lambda d'un modèle SBM
        D: int matrix, matrice des distance 
        n: int, nombre d'ndividus
        B: int, nombre de classes
    output :
        dicb : dict, dicionnaire de tous les facteurs binaires 
    """
    dicb={}
    for i in range (n-1) :
        for j in range (i+1,n) :
            d=D[i,j]
            dfac=math.factorial(d)
            psi =np.zeros((B,B))
            for b in range (B) :
                for bp in range (B) :
                    
                    psi[b,bp]=fact_value (Lam, int(d),dfac, int(b), int(bp))
            dicb[(i,j)]=psi
            
    return dicb

# ======================================================================


def unary_factors (n,B, D, alpha) :
    """ calcul des marginales bianires 
    input : 
        Lam : matrice paramètre lambda d'un modèle SBM
        D: int matrix, matrice des distance 
        n: int, nombre d'ndividus
        B: int, nombre de classes
    output :
        dicu : dict, dicionnaire de tous les facteurs unaires 
    """
    dicu={}
    for i in range (n) :
        phi=np.zeros(B)
        for b in range (B) :
            phi[b]=alpha[b]
        
        dicu[(i)]=phi
    
    return dicu


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


def thelam(case) :
    Lam=2*np.ones((3,3))
    if (case == 'assortative') :
        Lam[0,0], Lam[1,1],Lam[2,2]=2,3,4
        Lam[1,0],Lam[2,0],Lam[2,1], Lam[0,1],Lam[0,2],Lam[1,2]=10,10,10,10,10,10
    if (case == 'hierarchical'):
        Lam[0,0], Lam[1,1],Lam[2,2]=12,6,3
        Lam[1,0],Lam[2,0],Lam[2,1], Lam[0,1],Lam[0,2],Lam[1,2]=3,1,3,3,1,3
    if (case == 'ordered'):
        Lam[0,0], Lam[1,1],Lam[2,2]=10,10,10
        Lam[1,0],Lam[2,0],Lam[2,1], Lam[0,1],Lam[0,2],Lam[1,2]=6,2,6,6,2,6
    if (case == 'dissortative') :
        Lam[0,0 ], Lam[1,1 ],Lam[2,2 ]= 8 ,10,9
        Lam[1,0 ],Lam[2,0 ],Lam[2,1 ], Lam[0,1 ],Lam[0,2 ],Lam[1,2]=3,3,3,3,3,3
    if (case == 'coreperiphery') :
        Lam[0,0 ], Lam[1,1 ],Lam[2,2 ]= 2 ,4,10
        Lam[1,0 ],Lam[2,0 ],Lam[2,1 ], Lam[0,1 ],Lam[0,2 ],Lam[1,2]=5, 8,6,5, 8,6
    if (case == 'hierarchical2'):
        Lam[0,0 ], Lam[1,1 ],Lam[2,2 ]= 3,6, 10
        Lam[1,0 ],Lam[2,0 ],Lam[2,1 ], Lam[0,1 ],Lam[0,2 ],Lam[1,2]=10,15,15,10,15,15
    
    
    
    
    return Lam

def mod (n,B, Lam, alpha) :
    C,Z=assig_multinom_classes(n,alpha)
    D=mat_distance(n,Lam, C) 
    return C,D



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
        #self.case = which_case(Lam)
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
        
        




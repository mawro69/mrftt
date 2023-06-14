#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 10:07:24 2021

@author: abouabdallah
version 21.08.03
"""

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math

# ----------------------------------------------------------------------

def assig_multinom_classes(n,alpha) :
      """ finding the classes and building Z

      inputs :
              n : int nombre d'individus
              alpha : liste des probas d'appartenance aux classes
          outputs :
              C : vector \in R^{n} classes des individus
              Z : matrix \in R^{n,b} classes des individus


      """
      Z = np.random.multinomial(1, alpha, n)
      B = len(alpha)
      C = [0]*n
      for i in range (n) :
          for k in range (B):
              if (Z[i,k] ==1) :
                  C[i] = k
      C = [int(k) for k in C]

      return (C, Z)

# ----------------------------------------------------------------------

def mat_lam_aleatoire(B,dmax, dmin) :
     """ generates random matrix Lambda
     inputs :
         B : int nombre de classes
         dmax : int distance max
         dmin : int distance min
     return :
         lam : distances entre les classes (aleatoire)
     """
     Lam = np.zeros((B,B))
     for i in range (B-1):
         for j in range (i,B):
             a = np.random.randint(dmin,dmax)
             Lam[i,j] = a
             if j != i:
                 Lam[j,i] = Lam[i,j]

     return (Lam)

# ----------------------------------------------------------------------

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


# ======================================================================


def updated_binary_factors(dicu, dicb) :
    """ calcul des marginales bianires avec les marginales unaires 
    input : 
        dicu : dict, dicionnaire de tous les facteurs unaires 
        dicb : dict, dicionnaire de tous les facteurs binaires
    output :
        dic : dict, dicionnaire de tous les facteurs binaires avec les marginales unaires 
    """
    n=max(max(dicb.keys()))
    B=dicb[0,1].shape[0]
    dic=dicb
    for i in range (n-1):
        psi=dicb[i,i+1]
        phi=dicu[i]
        for b in range (B) :
            for bp in range (B) :
                psi[b,bp]= psi[b,bp]*phi[b]
        dic[i,i+1] = psi
    """psi=dicb[n-1,n]
    phi=dicu[n]
    for b in range (B) :
        for bp in range (B) :
            psi[b,bp]= psi[b,bp]*phi[b]
    dic[n-1,n] = psi
    """
    return dic 
    





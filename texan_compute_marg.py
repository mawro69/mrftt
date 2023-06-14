#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 26 13:49:38 2022

@author: abouabdallah
"""

import numpy as np




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
    dicu=np.zeros((n,B))
    for i in range (n) :
        phi=np.zeros(B)
        for b in range (B) :
            phi[b]=alpha[b]
        
        dicu[i,]=phi
    
    return dicu



def compute_tensor(q, dicb, dicu,n) :
    if n != 9 :
        return ("erreur developpé que pour n=9")
    Psi = np.zeros(q**(n),dtype=float)
    Psi.shape = (q,q,q,q,q,q,q,q,q)
    z=range(q)
    for a in range(len(z)) :
        z1=z[a]
        for b in range (len(z)) :
            z2=z[b]
            for c in range(len(z)) :
                z3=z[c]
                for d in range(len(z)) :
                    z4=z[d]
                    for e in range(len(z)) :
                        z5=z[e]
                        for f in range(len(z)) :
                            z6=z[f]
                            for g in range(len(z)) :
                                z7=z[g]
                                for h in range(len(z)) :
                                    z8=z[h]
                                    for i in range(len(z)) :
                                        z9=z[i]
                                        lis=[z1,z2,z3,z4,z5,z6,z7,z8,z9]
                                        psi1 =1
                                        for k in range (n -1) :
                                            for l in range (k+1,n):
                                                psi1=psi1*dicb[k,l][lis[k], lis[l]]
                                        for k in range (n) :
                                            psi1=psi1*dicu[k][lis[k]]
                                        Psi[a,b,c,d,e,f,g,h,i] = psi1
    
    
    return Psi






def compute_unaries(Psi) :
    
    unairesfull=np.zeros((len(Psi.shape),Psi.shape[0]))
    cst=np.sum(Psi)
    
    unairesfull[0,]= np.sum(Psi, axis=(1,2, 3, 4, 5, 6, 7, 8))/cst
    unairesfull[1,]= np.sum(Psi, axis=(0,2, 3, 4, 5, 6, 7, 8))/cst
    unairesfull[2,]= np.sum(Psi, axis=(1,0, 3, 4, 5, 6, 7, 8))/cst
    unairesfull[3,]= np.sum(Psi, axis=(1,2, 0, 4, 5, 6, 7, 8))/cst
    unairesfull[4,]= np.sum(Psi, axis=(1,2, 3, 0, 5, 6, 7, 8))/cst
    unairesfull[5,]= np.sum(Psi, axis=(1,2, 3, 4, 0, 6, 7, 8))/cst
    unairesfull[6,]= np.sum(Psi, axis=(1,2, 3, 4, 5, 0, 7, 8))/cst
    unairesfull[7,]= np.sum(Psi, axis=(1,2, 3, 4, 5, 6, 0, 8))/cst
    unairesfull[8,]= np.sum(Psi, axis=(1,2, 3, 4, 5, 6, 7, 0))/cst
    
    return unairesfull


def binaries_margins(Psi,n) :
    ll= [i for i in range(n)]
    ll2=ll
    binairesfull=dict()
    w=np.sum(Psi)
    for i in range (n-1) :
        for j in range (i+1,n) :
            ll2.remove(i)
            ll2.remove(j)
            if ([i, j]==[6,7]) : 
                print(ll2)
            binairesfull[(i,j)]=np.sum(Psi, axis=tuple(ll2))/w
            ll2=[i for i in range(n)]

    return binairesfull


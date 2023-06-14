#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 28 14:25:23 2022

@author: abouabdallah
"""



import numpy as np
import math

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



def objectif (alpha, D, lam, n, Q, dicu):
    
    
    value=np.zeros((n,Q))
    for i in range (n) :
        for q in range (Q) :
            
            
            value[i,q] = alpha[q]
            for j in range (n) :
                if (j!=i) :
                    for l in range (Q) :
                        a=((lam[q,l]**D[i,j])*math.exp(-lam[q,l])/math.factorial(D[i,j]))
                        a=a**dicu[j][l]
                        value[i,q]=value[i,q]*a
    for i in range (n) :
        f=np.sum(value[i,:])
        for q in range (Q) :
            
            value[i,q]=value[i,q]/f
    
    return value 






def tolerance (vec, vec2, epsilone) :
    a= np.linalg.norm(vec-vec2)
    b =np.linalg.norm(vec2)
    res = a/b<epsilone
    return res



def tolerance_calculee (vec, vec2) :
    print(vec-vec2)
    a= np.linalg.norm(vec-vec2)
    print(a)
    b =np.linalg.norm(vec2)
    print(b)
    res = a/b
    return res





def pt_fixe (alpha, D, lam, n, Q, epsilone, dicu ):
    dic=dicu
    a=0
    #print(a)
    dicu1=objectif (alpha, D, lam, n, Q, dicu)
    while (tolerance(dicu1, dic, epsilone)==False and a<2000) :
        a=a+1
       # print(a)
        dic=dicu1
        dicu1=objectif (alpha, D, lam, n, Q, dicu1)
        #print("iteration")
 #       print(dicu1)
        #print(tolerance_calculee (dic, dicu1))
        if (a==1999) :
            print("no conv")
            dicu1=np.zeros((n,Q))
            break
    
    
    
    return dicu1



def pt_fixe_allit (alpha, D, lam, n, Q, epsilone, dicu ):
    dic=dicu
    a=0
    #print(a)
    lis=[]
    dicu1=objectif (alpha, D, lam, n, Q, dicu)
    lis.append(dicu1)

    while (tolerance(dicu1, dic, epsilone)==False and a<2000) :
        a=a+1
       # print(a)
        dic=dicu1
        dicu1=objectif (alpha, D, lam, n, Q, dicu1)
        lis.append(dicu1)
        #print("iteration")
 #       print(dicu1)
        #print(tolerance_calculee (dic, dicu1))
        if (a==1999) :
            print("no conv")
            dicu1=np.zeros((n,Q))
            break
    
    
    
    return lis[1900:1999]



def normalisation(res) :
    n,q=res.shape
    for k in range (n) :
        for l in range (q):
            aaa=np.sum(np.sum(res[k,:]))
            if (aaa > 0) :
                res[k,l]=res[k,l]/aaa
            if (res[k,l]<1e-10):
                res[k,l]=0
    return res





def normalisation2(res) :
    n,q=res.shape
    for k in range (n) :
        div=np.sum(res[k,:])
        for l in range (q):
            res[k,l]=res[k,l]/div
            if (res[k,l]<1e-10):
                res[k,l]=0
    return res



def binaries (res) :
    n,q=res.shape
    dicb={}
    for i in range (n-1):
        for j in range (i+1, n) :
            mat=np.zeros((q,q))
            for k in range (q):
                for l in range(q) :
                    mat[k,l]=res[i,k]*res[j,l]
            dicb[i,j]=mat
    return dicb



def unaires_binaires_ptfixe(alpha, D, Lam, n, B) :
    dicu=unary_factors (n,B, D, alpha)            

#    obj=ptfixe.objectif (alpha, D, Lam, n, B, dicu)
    ob2=pt_fixe (alpha, D, Lam, n, B, 1e-2, dicu )
    res=normalisation(ob2)
    res2=binaries (res)
    
    
    return res, res2
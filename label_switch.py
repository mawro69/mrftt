#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 17:00:25 2023

@author: abouabdallah
"""




import numpy as np


def app_classes(matclases, n,q ) :
    classes=[]
    for k in range (n) :
        classes.append(list(matclases[k,]).index(max(matclases[k,])))
    return classes


def contingeance(partition1, partition2) :
    conting=np.zeros((int(max(partition1) +1),int(max(partition2) +1)))
    for k in range (len(partition1)) :
        i=int(partition1[k])
        j=int(partition2[k])
        conting[i,j]=conting[i,j]+1
    
    return conting


def rotation(vec1, vec2) :
    mat=contingeance(vec1, vec2) 
    q=mat.shape[0]
    rota=[]
    for k in range (q) :
        for l in range (q) :
            if (mat[k,l] !=  0) :
                rota.append(l)
    
    return rota



def rotated_res(res,C, n, B) :
    vec1=C
    vec2=app_classes(res, n,B) 
    rota = rotation(vec1, vec2)
    res=res[:,rota]
    
    return res




"""def label_switch_mcmc(dicmcmc,C, n, q) :
    dicmcmc2=np.zeros((dicmcmc.shape))
    dicmcmc3= dicmcmc.astype(int)
    rota=rotation(dicmcmc3[0,:], C)
    for k in range (dicmcmc.shape[0]) :
        for l in range (dicmcmc.shape[1]) :
            

    for k in range (dicmcmc.shape[0]) :
        clas=dicmcmc3[k,:]
        type(clas)
        rotation(clas, C)
        dicmcmc2[k,:]=dicmcmc3[k,:][rotation(clas, C)]
    return dicmcmc2
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 10:48:15 2022

@author: abouabdallah
"""

import pandas as pd

import seaborn as sns
import numpy as np

def data_to_df(uns1,uns2, title1, title2) :#the format vem, mcmc, TT
    n, q = uns1.shape
    df = np.zeros((n*q, 2))
    a=0
    for k in range (n) : 
        for l in range (q) :
            
            #print(a)
            df[a,:]= [uns1[k,l],uns2[k,l]]
            a+=1

#    df = pd.DataFrame(df)
    my_frame = pd.DataFrame(data={title1: df[:,0],title2: df[:,1]})
    return my_frame


def marginales_df(uns1, bns1, bns2, title1, title2) :
    n, q = uns1.shape
    df = np.zeros((int(n*(n-1)*q*q/2), 2))
    a=0
    for i in range (n-1) : 
        for j in range (i+1,n) :
            for k in range (q) :
                for l in range (q) :
                    df[a,:]= [bns1[i,j][k,l],bns2[i,j][k,l]]
                    a+=1
    my_frame = pd.DataFrame(data={title1: df[:,0],title2: df[:,1]})

    return my_frame



def swarm(uns1,uns2, title1, title2, whichones) :
    if (whichones==1) :
        df=data_to_df(uns1,uns2, title1, title2)
        swarm_plot = sns.pairplot(df)
    if (whichones==2) :
        df=marginales_df(uns1,uns2, title1, title2)
        swarm_plot = sns.pairplot(df)
    else :
        return ("Impossible")
    return swarm_plot



def data_to_df2(uns1,uns2, uns3, title1, title2, title3) :#the format vem, mcmc, TT
    n, q = uns1.shape
    df = np.zeros((n*q, 3))
    a=0
    for k in range (n) : 
        for l in range (q) :
            
            #print(a)
            df[a,:]= [uns1[k,l],uns2[k,l], uns3[k,l]]
            a+=1

#    df = pd.DataFrame(df)
    my_frame = pd.DataFrame(data={title1: df[:,0],title2: df[:,1],title3: df[:,2]})
    return my_frame


def marginales_df2(uns1, bns1, bns2, bns3, title1, title2, title3) :
    n, q = uns1.shape
    df = np.zeros((int(n*(n-1)*q*q/2), 3))
    a=0
    for i in range (n-1) : 
        for j in range (i+1,n) :
            for k in range (q) :
                for l in range (q) :
                    df[a,:]= [bns1[i,j][k,l],bns2[i,j][k,l],bns3[i,j][k,l]]
                    a+=1
    my_frame = pd.DataFrame(data={title1: df[:,0],title2: df[:,1], title3: df[:,2]})

    return my_frame


def swarm2(uns1,uns2,uns3, title1, title2, title3, whichones) :
    if (whichones==1) :
        df=data_to_df2(uns1,uns2, uns3, title1, title2, title3)
        swarm_plot = sns.pairplot(df)
    #if (whichones==2) :
    #    df=marginales_df2(uns1, bns1, bns2, bns3, title1, title2, title3)
    #    swarm_plot = sns.pairplot(df)
    else :
        return ("Impossible")
    return swarm_plot

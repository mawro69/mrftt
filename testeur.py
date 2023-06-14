#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 14:05:15 2023

@author: abouabdallah
"""

import pandas  as pd 
import numpy as np
import comp_mar_mcmc

df= pd.read_csv('brollon2/Dfnn60assortative112.csv')


df.shape



C=np.loadtxt('C50nnassortative112.txt')

C.shape

marguna=np.zeros((50, 3))

for k in range (n) :
    for l in range (3) :
        marguna[k,l]=df['TT'][3*k+l]

df['TT']


comp_mar_mcmc.app_classes(marguna, 50, 3)==C






df= pd.read_csv('brollon2/Dfnn60hierarchical2112.csv')


df.shape



C=np.loadtxt('C50nnhierarchical2112.txt')

C.shape

marguna=np.zeros((50, 3))

for k in range (n) :
    for l in range (3) :
        marguna[k,l]=df['TT'][3*k+l]

df['TT']


comp_mar_mcmc.app_classes(marguna, 50, 3)==C


df= pd.read_csv('brollon2/Dfnn60hierarchical112.csv')


df.shape



C=np.loadtxt('C50nnhierarchical112.txt')

C.shape

marguna=np.zeros((50, 3))

for k in range (n) :
    for l in range (3) :
        marguna[k,l]=df['TT'][3*k+l]

df['TT']


comp_mar_mcmc.app_classes(marguna, 50, 3)==C






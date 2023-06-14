#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  4 10:30:51 2023

@author: abouabdallah
"""


import mrf_tt_all

import numpy as np



#import tt
import time

#import bases
#import model_gen

#import cores_melting
#import compressed_mrf_tt
#import symbolic_tt_mrf



#####generer un model
alpha=[1/3,1/3,1/3]
case="assortative"
Lam=mrf_tt_all.model_gen.thelam(case)

n=9
q=3
Q=q
C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   


SBM=mrf_tt_all.model_gen.model(D, Lam,alpha)

D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q


dicu=SBM.dicu
dicb=SBM.dicb

dic=SBM.dic

dic=mrf_tt_all.model_gen.updated_binary_factors2(dicu, dicb)
print(D)
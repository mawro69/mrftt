#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 19:06:30 2023

@author: abouabdallah
"""

import sys
import comp_mar_mcmc
import matplotlib.pyplot as plt
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time

import visualisation

import label_switch


case ="assortative"

ff=1


rs1=0
alpha=[1/3,1/3,1/3]
#case="coreperiphery"
n=50
q=3
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
dd=0


C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   
SBM=mrf_tt_all.model_gen.model(D, Lam, alpha)
dic10=SBM.multiply(a)
D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q
dic10=SBM.dica
dicu=comp_mar_mcmc.unary_factors(n,q,D, alpha)
dicb=comp_mar_mcmc.binary_factors (n,q, D, Lam) 
epsilone=1e-1
dicu=pt_fixe.unary_factors (n,q, D, alpha)            
rs=pt_fixe.pt_fixe (alpha, D, Lam, n, q, epsilone, dicu )
rs1=np.linalg.norm(rs)
print(rs1)
if (rs1 ==0) :
    rs2=pt_fixe.pt_fixe_allit (alpha, D, Lam, n, q, epsilone, dicu )



vva=[]
for k in range (1,len(rs2)) :
    vva.append(np.linalg.norm(rs2[k]))


plt.hist(vva)
plt.savefig("hisosci.png")


plt.plot(vva)
plt.savefig("plotpsci.png")


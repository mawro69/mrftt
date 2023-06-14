#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 16:12:11 2023

@author: abouabdallah
"""






import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time

import visualisation

import label_switch

rs1=0

alpha=[1/3,1/3,1/3]
#alpha=[0.4,0.35,0.25]

case="dissortative"
n=50
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
dd=0

while (rs1 ==  0) :
    C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   
    print(np.array(C))
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
    dd+=1
    print(dd)



res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)


print(np.array(comp_mar_mcmc.app_classes(res, n ,q)))

print(np.array(C))



label_switch.contingeance(np.array(comp_mar_mcmc.app_classes(res, n ,q)),np.array(C))
res_switched=label_switch.rotated_res(res,C, n, q)

dicb=comp_mar_mcmc.all_fac (dicb)
dicu= comp_mar_mcmc.dicu_tomat(dicu,n,q)


dicmcmc=comp_mar_mcmc.echantilloneur_gibbs_allsteps(n, q, dicu, dicb, Niter=500)


dicmcmc=comp_mar_mcmc.echantilloneur_gibbs_allsteps_stockage(n, q, dicu, dicb, Niter=9000)

Niter=35000
firstone= 9000
latence = 500



dicmcmc=comp_mar_mcmc.echantillons_token(n, q, dicu, dicb, Niter, firstone, latence)


unsmcmc=comp_mar_mcmc.margunaires(dicmcmc)


np.array(comp_mar_mcmc.app_classes(unsmcmc, n ,q))

np.array(comp_mar_mcmc.app_classes(res, n ,q))
np.array(C)


label_switch.contingeance(np.array(comp_mar_mcmc.app_classes(unsmcmc, n ,q)),np.array(C))
res_switched_mcmc=label_switch.rotated_res(unsmcmc,C, n, q)


np.array(C)


ao=time.time()
aaa=comp_mar_mcmc.Lechant(n,q,dicu, dicb, Niter=500, L=1)
ao2=time.time()

binaries= comp_mar_mcmc.margbinaires(aaa)
unaries1 = comp_mar_mcmc.margunaires(aaa)

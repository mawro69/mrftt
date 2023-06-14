#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 09:57:20 2023

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


import texan_compute_marg

rs1=0

num=0

alpha=[1/3,1/3,1/3]
#alpha=[0.4,0.35,0.25]

case="coreperiphery"
n=9
#print(n)
q=3
#print(q)


alltimes=[]

for num in range (100) :
    print(num)
    a=10
    Lam=mrf_tt_all.model_gen.thelam(case)
    dd=0
    rs1 = 0
    while (rs1 ==  0) :
        C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   #print(np.array(C))
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
    np.savetxt('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/D9_'+ str(num) + '_'+ case +'.txt',D,fmt='%.2f')
    np.savetxt('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/C9_'+ str(num) + '_'+ case +'.txt',C,fmt='%.2f')
    Ctest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/C9_'+ str(num) + '_'+ case +'.txt')
    Dtest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/D9_'+ str(num) + '_'+ case +'.txt')
    start_time = time.time()
    res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)


    end_time= time.time()
    ti=end_time-start_time
    tipf=ti


    ll=list(dic10.keys())
    factors =mrf_tt_all.model_org(n, alpha, Lam, a)
    roundprec=5e-2
    rmax=q**3
    Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
    Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)
    start_time = time.time()
    Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 27, 1)    
    end_time= time.time()
    ti=end_time-start_time
    Zeq2=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    
    Zeq3=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 3, 1)    
    start_time = time.time()

    dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq)
    binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)
    dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
    unaires=dicub['unary']
    binaires =dicub['binary']
    end_time= time.time()
    ti=end_time-start_time
    titt=ti
    lol=mrf_tt_all.compute_marginales.other_format(dicub)
    uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)
    df=visualisation.data_to_df(uns,res, "TT", "Pf")
    df.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")
    df2.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)
    start_time = time.time()
    tens=texan_compute_marg.compute_tensor(q, dicb, dicu, n)
    unsexact=texan_compute_marg.compute_unaries(tens)
    binsexact=texan_compute_marg.binaries_margins(tens, n)
    end_time= time.time()
    ti=end_time-start_time
    tiexact=ti

    dfall=visualisation.data_to_df2(uns,res,unsexact, "TT", "Pf", "Exact")
    dfall.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2all=visualisation.marginales_df2(uns,lol,res2, binsexact, "TT", "Pf", "Gibbs")
    df2all.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    alltimes.append([tiexact, titt, tipf])


np.savetxt('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/times'+ '_'+ case +'.txt',alltimes,fmt='%.2f')



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan  8 19:04:32 2023

@author: abouabdallah
"""



import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time

import visualisation

rs1=0

alpha=[1/3,1/3,1/3]
#alpha=[0.4,0.35,0.25]

case="ordered"
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



np.savetxt('/home/abouabdallah/Bureau/docurespres/article/resultats/D50nn'+ case +'.txt',D,fmt='%.2f')



np.savetxt('/home/abouabdallah/Bureau/docurespres/article/resultats/C50nn'+ case +'.txt',C,fmt='%.2f')



Ctest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/article/resultats/C50nn'+ case +'.txt')
Dtest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/article/resultats/D50nn'+ case +'.txt')



res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)


ll=list(dic10.keys())
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**2

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
#stock_all_Ai_q_compressed_k2(factors, q, n, k)


print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
start_time = time.time()
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 27, 1)    
end_time= time.time()
ti=end_time-start_time

Zeq2=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    

Zeq3=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 3, 1)    

print(Zeq3, Zeq2, Zeq)
print(["computation time", ti])

print("End Z")

q
start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,q,Zeq)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time

print(["computation time", ti])


print("end margins")


#unemar=mrf_tt_all.compute_marginales.compute_margin_smart(Ais1eq, Bis1eq, 15, 15, n, q, roundprec, rmax) 

uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)



print("Comparaisons preliminaires")

print("CM vs classes")
print(C==np.array(comp_mar_mcmc.app_classes(res, n,q )))

print("CM vs TT")

print(np.array(comp_mar_mcmc.app_classes(res, n,q ))==np.array(comp_mar_mcmc.app_classes(uns, n,q )))

print("Classes vs TT")
print(C==np.array(comp_mar_mcmc.app_classes(uns, n,q )))



df=visualisation.data_to_df(uns,res, "TT", "Pf")


df.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/Dfnn60'+case+'.csv')


df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")

df2.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/Dfnn60'+case+'bin.csv')

swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)



dicb=comp_mar_mcmc.all_fac (dicb)
dicu= comp_mar_mcmc.dicu_tomat(dicu,n,q)

ao=time.time()
aaa=comp_mar_mcmc.Lechant(n,q,dicu, dicb, Niter=500, L=1000)
ao2=time.time()

binaries= comp_mar_mcmc.margbinaires(aaa)
unaries1 = comp_mar_mcmc.margunaires(aaa)
unaries1
np.array(comp_mar_mcmc.app_classes(unaries1, n,q ))==C

df3=visualisation.data_to_df2(uns,res,unaries1, "TT", "Pf", "Gibbs")

swarm_plot=visualisation.swarm(uns,unaries1, "TT", "Pf", 1)
swarm_plot=visualisation.swarm2(uns,res,unaries1, "TT", "Pf","Gibbs", 1)



dfall=visualisation.data_to_df2(uns,res, unaries1, "TT", "Pf", "Gibbs")


dfall.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/Dfnn60'+case+'.csv')


df2all=visualisation.marginales_df2(uns,bins,res2, binaries, "TT", "Pf", "Gibbs")

df2all.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/Dfnn60'+case+'bin.csv')

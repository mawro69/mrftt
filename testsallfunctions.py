#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 22 16:16:08 2022

@author: abouabdallah
"""

##################################### Model generator 


import comp_mar_mcmc
import pt_fixe


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




""" texans methode """
ducut=mrf_tt_all.texan_compute_marg.unary_factors(n,Q, D, alpha)

testtensor=mrf_tt_all.texan_compute_marg.compute_tensor(Q, dicb, ducut,n)

unaries=mrf_tt_all.texan_compute_marg.compute_unaries(testtensor)

binaries=mrf_tt_all.texan_compute_marg.binaries_margins(testtensor,n)
binaries[7,8].sum()



""" Novikov"""


psi= dic[0,1]
gij=mrf_tt_all.tt_format2(psi)
nonessent=mrf_tt_all.ajout_des_dimension_non_ess2(gij, [0,1], n)

nonessentall=mrf_tt_all.tt_factor_general(dicb)

Ais=mrf_tt_all.stock_all_Ai_q(nonessentall, q, n)
Ais=mrf_tt_all.compute_Ai_per_method(dic, 0, None, None, q, n) 

rmax=9999
roundprec= 5e-2
Bis=mrf_tt_all.compute_Bi_all(Ais, q, n, roundprec)

Zeq=mrf_tt_all.compute_W_simple(Bis, roundprec, 9999, 1)    

start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis,roundprec,rmax,Zeq)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais, q,roundprec, rmax, Zeq)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time


##### n=40


alpha=[1/3,1/3,1/3]
case="assortative"
n=40
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   

SBM=mrf_tt_all.model_gen.model(D, Lam, alpha)
dic10=SBM.multiply(a)
D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q
dic10=SBM.dica

factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**3

print("without cores melting")



print("begin Ai")
Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 0, None, None, q, n) 
print("end Ai")

print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")

Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    

print("End Z")


Lzeq=mrf_tt_all.compute_marginales.compute_log_Z(Zeq, n)


print("begin margins")

start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time

print(["computation time", ti])
print("end margins")







print("with cores melting")



print("begin Ai")
Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n) 
print("end Ai")

print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")

Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    

print("End Z")


Lzeq=mrf_tt_all.compute_marginales.compute_log_Z(Zeq, n)


print("begin margins")

start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time

print(["computation time", ti])
print("end margins")





#test uniforme volume 2

case="assortative"
n=45
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   

SBM=mrf_tt_all.model_gen.model(D, Lam, alpha)
dic10=SBM.multiply(a)
D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q
dic10=SBM.dica

factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**3


k=3

Ais1eq=mrf_tt_all.stock_all_Ai_q_compressed_k2(factors, q, n, k)



print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
Zeq1=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
print("End Z")




##########uniforme volumes
print("begin Ai")

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 5, None, q, n)  

print("end Ai")

Ais1eq= mrf_tt_all.stock_all_Ai_q_compressed_k(factors, q, n, k)


print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
print("End Z")




##########uniforme volumes
print("begin Ai")

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 2, 0, 300, q, n)  

print("end Ai")


print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
Zeq2=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
print("End Z")


print([Zeq, Zeq1, Zeq2])






######################################################################################

""" calcul des marginales """
rmax2=Q**2

print("begin margins")

alpha=[1/3,1/3,1/3]

case="dissortative"
n=62
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
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

res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)

C
ll=list(dic10.keys())
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**2

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 3, None, q, n)
#stock_all_Ai_q_compressed_k2(factors, q, n, k)



print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
start_time = time.time()
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
end_time= time.time()
ti=end_time-start_time

m=SBM.m
print((ti*m*q**2)/3600)
print(["computation time", ti])

print("End Z")


start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time
n
print(["computation time", ti])


print("end margins")

unaires[0]

uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)


C==np.array(comp_mar_mcmc.app_classes(res, n,q ))

np.array(C)
np.array(comp_mar_mcmc.app_classes(res, n,q ))
np.array(comp_mar_mcmc.app_classes(uns, n,q ))

import visualisation


df=visualisation.data_to_df(uns,res, "TT", "Pf")

swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)

C==comp_mar_mcmc.app_classes(uns, n,q )

SBM.dicu
SBM.dicb


dicb=comp_mar_mcmc.all_fac (dicb)
dicu= comp_mar_mcmc.dicu_tomat(dicu,n,q)

ao=time.time()
aaa=comp_mar_mcmc.Lechant(n,q,dicu, dicb, Niter=1000, L=1000)
ao2=time.time()


binaries= comp_mar_mcmc.margbinaires(aaa)
unaries1 = comp_mar_mcmc.margunaires(aaa)

C==comp_mar_mcmc.app_classes(unaries1, n,q )
np.array(C)
np.array(comp_mar_mcmc.app_classes(unaries1, n,q ))



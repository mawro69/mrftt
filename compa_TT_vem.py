#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 29 12:20:19 2022

@author: abouabdallah
"""

import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time



alpha=[1/3,1/3,1/3]

case="coreperiphery"
n=50
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   

np.savetxt('/home/abouabdallah/Bureau/docurespres/D50coreperiphery.txt',D,fmt='%.2f')



np.savetxt('/home/abouabdallah/Bureau/docurespres/C50coreperiphery.txt',C,fmt='%.2f')



Ctest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/C50coreperiphery.txt')
Dtest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/D50coreperiphery.txt')



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
np.array(comp_mar_mcmc.app_classes(res, n,q ))==np.array(comp_mar_mcmc.app_classes(uns, n,q ))

import visualisation


df=visualisation.data_to_df(uns,res, "TT", "Pf")


df.to_csv('/home/abouabdallah/Bureau/docurespres/Df'+case+'.csv')


df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")

df2.to_csv('/home/abouabdallah/Bureau/docurespres/Df'+case+'bin.csv')

swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)



##########################################################################################
""" assortative """

import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time


rs1=0


alpha=[1/3,1/3,1/3]
alpha=[0.3,0.2,1/2]

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

C


np.savetxt('/home/abouabdallah/Bureau/docurespres/D50nn'+ case +'.txt',D,fmt='%.2f')



np.savetxt('/home/abouabdallah/Bureau/docurespres/C50nn'+ case +'.txt',C,fmt='%.2f')



Ctest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/C50nn'+ case +'.txt')
Dtest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/D50nn'+ case +'.txt')






res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)

C
ll=list(dic10.keys())
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**3

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
#stock_all_Ai_q_compressed_k2(factors, q, n, k)



print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
start_time = time.time()
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
end_time= time.time()
ti=end_time-start_time

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

unaires[0][16]
unaires[1][16]
unaires[2][16]


#unemar=mrf_tt_all.compute_marginales.compute_margin_smart(Ais1eq, Bis1eq, 15, 15, n, q, roundprec, rmax) 

uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)


C==np.array(comp_mar_mcmc.app_classes(res, n,q ))

np.array(C)
np.array(comp_mar_mcmc.app_classes(res, n,q ))==np.array(comp_mar_mcmc.app_classes(uns, n,q ))

import visualisation


df=visualisation.data_to_df(uns,res, "TT", "Pf")


df.to_csv('/home/abouabdallah/Bureau/docurespres/Dfnn'+case+'.csv')


df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")

df2.to_csv('/home/abouabdallah/Bureau/docurespres/Dfnn'+case+'bin.csv')

swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)





#####################
import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time


rs1=0


alpha=[1/3,1/3,1/3]
alpha=[0.4,0.35,0.25]

case="hierarchical"
n=50
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
dd=0
while (rs1 ==  0) :
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
    dd+=1
    print(dd)

C


np.savetxt('/home/abouabdallah/Bureau/docurespres/D50nn'+ case +'.txt',D,fmt='%.2f')



np.savetxt('/home/abouabdallah/Bureau/docurespres/C50nn'+ case +'.txt',C,fmt='%.2f')



Ctest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/C50nn'+ case +'.txt')
Dtest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/D50nn'+ case +'.txt')






res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)

C
ll=list(dic10.keys())
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**3

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
#stock_all_Ai_q_compressed_k2(factors, q, n, k)



print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
start_time = time.time()
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
end_time= time.time()
ti=end_time-start_time

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

unaires[0][16]
unaires[1][16]
unaires[2][16]


#unemar=mrf_tt_all.compute_marginales.compute_margin_smart(Ais1eq, Bis1eq, 15, 15, n, q, roundprec, rmax) 

uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)


C==np.array(comp_mar_mcmc.app_classes(res, n,q ))

np.array(C)
np.array(comp_mar_mcmc.app_classes(res, n,q ))==np.array(comp_mar_mcmc.app_classes(uns, n,q ))

import visualisation


df=visualisation.data_to_df(uns,res, "TT", "Pf")


df.to_csv('/home/abouabdallah/Bureau/docurespres/Dfnn'+case+'.csv')


df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")

df2.to_csv('/home/abouabdallah/Bureau/docurespres/Dfnn'+case+'bin.csv')

swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)





######################n=60
import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time



rs1=0


alpha=[1/3,1/3,1/3]
alpha=[0.4,0.35,0.25]

case="hierarchical"
n=55
print(n)
q=3
print(q)
a=10
Lam=mrf_tt_all.model_gen.thelam(case)
dd=0
while (rs1 ==  0) :
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
    dd+=1
    print(dd)

C


np.savetxt('/home/abouabdallah/Bureau/docurespres/D60nn'+ case +'.txt',D,fmt='%.2f')



np.savetxt('/home/abouabdallah/Bureau/docurespres/C60nn'+ case +'.txt',C,fmt='%.2f')



Ctest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/C60nn'+ case +'.txt')
Dtest=np.loadtxt('/home/abouabdallah/Bureau/docurespres/D60nn'+ case +'.txt')


Ctest.shape



res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)

C
ll=list(dic10.keys())
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=1e-1
rmax=q**3

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
#stock_all_Ai_q_compressed_k2(factors, q, n, k)


rmax=9
print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")
start_time = time.time()
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
end_time= time.time()
ti=end_time-start_time

Zeq2=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 27, 1)    

Zeq3=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 3, 1)    

print(["computation time", ti])

print("End Z")


start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,3,Zeq)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time
n
print(["computation time", ti])


print("end margins")

unaires[0][16]
unaires[1][16]
unaires[2][16]


#unemar=mrf_tt_all.compute_marginales.compute_margin_smart(Ais1eq, Bis1eq, 15, 15, n, q, roundprec, rmax) 

uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)


C==np.array(comp_mar_mcmc.app_classes(res, n,q ))

np.array(C)
np.array(comp_mar_mcmc.app_classes(res, n,q ))==np.array(comp_mar_mcmc.app_classes(uns, n,q ))

import visualisation


df=visualisation.data_to_df(uns,res, "TT", "Pf")


df.to_csv('/home/abouabdallah/Bureau/docurespres/Dfnn60'+case+'.csv')


df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")

df2.to_csv('/home/abouabdallah/Bureau/docurespres/Dfnn60'+case+'bin.csv')

swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)
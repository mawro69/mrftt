#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 22:26:20 2022

@author: abouabdallah
"""

import mrf_tt_all
import numpy as np
import time

import bases
import model_gen

import cores_melting
import compressed_mrf_tt

import compute_marginales






#########generer notre model 
alpha=[1/3,1/3,1/3]
case="assortative"
Lam=model_gen.thelam(case)
n=42
q=3
C,D= model_gen.mod (n,q, Lam, alpha)   


SBM=model_gen.model(D, Lam, alpha)



dic10=SBM.multiply(10)

D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q


dic = SBM.dic
############################# compute factors 
a=10
factors =mrf_tt_all.model_org(n, alpha, Lam, a)




###################### compute z and margins without melting

Ais=mrf_tt_all.stock_all_Ai_q(factors, q, n)


roundprec=5e-2
Bis=mrf_tt_all.compute_Bi_all(Ais, q, n, roundprec)


rmax=9
ff=1
n
start_time = time.time()

Z=mrf_tt_all.compute_W_simple(Bis, roundprec, rmax, ff)    
end_time= time.time()

ti=end_time-start_time


Ais1eq=mrf_tt_all.compute_Ai_per_method(dic, 0, 3, 0, q, n) 

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    

start_time = time.time()
dicbij=compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq)

binas=compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq)

dicub=compute_marginales.bina_una(binas) 


unaires=dicub['unary']
binaires =dicub['binary']

end_time= time.time()

ti=end_time-start_time
ll=list(dic.keys())

unairesf,binairesf =compute_marginales.formatptfixe(unaires, binaires, ll , q)



####################################### uniform number of cores melting


Ais2eq=mrf_tt_all.compute_Ai_per_method(dic, 1, 3, 0, q, n) 
Bis2eq=mrf_tt_all.compute_Bi_all(Ais2eq, q, n, roundprec)
Z2eq=mrf_tt_all.compute_W_simple(Bis2eq, roundprec, 9, 1)    



rmax=9
ff=1
n



start_time = time.time()
print("begin graph")
dicbij=compute_marginales.graphdesBi(Bis2eq,roundprec,rmax,Zeq)

print("end graph")

print("begin margins")
binas=compute_marginales.all_binaries_all_states (dicbij,Ais2eq, q,roundprec, rmax, Zeq)
print("end margins")

dicub=compute_marginales.bina_una(binas) 


unaires=dicub['unary']
binaires =dicub['binary']

end_time= time.time()

ti=end_time-start_time
ll=list(dic.keys())

unairesf,binairesf =compute_marginales.formatptfixe(unaires, binaires, ll , q)







############################################## uniform volumes melting 
Ais3=mrf_tt_all.compute_Ai_per_method2(dic, q, n) # choix 2, 5 et 300
Bis3=mrf_tt_all.compute_Bi_all(Ais3, q, n, roundprec)
Z3=mrf_tt_all.compute_W_simple(Bis3, roundprec, 9, 1)    




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 11:18:28 2022

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

# test sur les fonctions de bases
d = 5
listoflis = bases.listoflist(d)


vec1 = np.random.rand(6)
vec2 = np.random.rand(6)

vecres = bases.multi(vec1, vec2)



#########generer notre model 
alpha=[1/3,1/3,1/3]
case="assortative"
Lam=model_gen.thelam(case)
n=42
q=3
C,D= model_gen.mod (n,q, Lam, alpha)   


SBM=model_gen.model(D, Lam, alpha)

type(SBM.case)
SBM.alpha

type(SBM.Lam)
type(SBM.D)
dic=SBM.dic
dic[7,8]

dic10=SBM.multiply(1)
dic[7,8]

dic10=SBM.multiply(10)
dic[7,8]
SBM.dica[7,8]
SBM.dic[7,8]

D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q



dir(SBM)
############################# compute factors 
a=10
factors =mrf_tt_all.model_org(n, alpha, Lam, a)

###################### compute z without melting

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


Ais1=mrf_tt_all.compute_Ai_per_method2(dic, q, n)#choix 0
Bis1=mrf_tt_all.compute_Bi_all(Ais1, q, n, roundprec)
Z=mrf_tt_all.compute_W_simple(Bis1, roundprec, 9, 1)    

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic, 0, 3, 0, q, n) 
Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    


Ais2=mrf_tt_all.compute_Ai_per_method2(dic, q, n)#choix 1, 3
Bis2=mrf_tt_all.compute_Bi_all(Ais2, q, n, roundprec)
Z2=mrf_tt_all.compute_W_simple(Bis2, roundprec, 9, 1)    

Ais2eq=mrf_tt_all.compute_Ai_per_method(dic, 1, 3, 0, q, n) 
Bis2eq=mrf_tt_all.compute_Bi_all(Ais2eq, q, n, roundprec)
Z2eq=mrf_tt_all.compute_W_simple(Bis2eq, roundprec, 9, 1)    



Ais3=mrf_tt_all.compute_Ai_per_method2(dic, q, n) # choix 2, 5 et 300
Bis3=mrf_tt_all.compute_Bi_all(Ais3, q, n, roundprec)
Z3=mrf_tt_all.compute_W_simple(Bis3, roundprec, 9, 1)    

Ais3eq=mrf_tt_all.compute_Ai_per_method(dic, 2, 0, 300, q, n) 
Bis3eq=mrf_tt_all.compute_Bi_all(Ais3eq, q, n, roundprec)
Z3eq=mrf_tt_all.compute_W_simple(Bis3eq, roundprec, 9, 1)    



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


n
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




Dtest=np.loadtxt('/home/abouabdallah/Bureau/programmettem/D18dissort.txt')
Ctest=np.loadtxt('/home/abouabdallah/Bureau/programmettem/C18dissort.txt')

alpha=[1/3,1/3,1/3]
case="dissortative"
Lam=model_gen.thelam(case)
n=18
q=3


SBM=model_gen.model(Dtest, Lam, alpha)

type(SBM.case)
SBM.alpha

type(SBM.Lam)
type(SBM.D)
dic=SBM.dic
dic[7,8]

dic10=SBM.multiply(1)
dic[7,8]


D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q

dir(SBM)
############################# compute factors 
a=1
factors =mrf_tt_all.model_org(n, alpha, Lam, a)


Ais1eq=mrf_tt_all.compute_Ai_per_method(dic, 0, 3, 0, q, n) 
Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    


def app_classes(matclases, n,q ) :
    classes=[]
    for k in range (n) :
        classes.append(list(matclases[k,]).index(max(matclases[k,])))
    return classes





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
Ctest
app_classes(unairesf, n, q)




Ctest=np.loadtxt('/home/abouabdallah/Bureau/programmettem/C18assort.txt')
Dtest=np.loadtxt('/home/abouabdallah/Bureau/programmettem/D18assort.txt')


alpha=[1/3,1/3,1/3]
case="assortative"
Lam=model_gen.thelam(case)
n=18
q=3


SBM=model_gen.model(Dtest, Lam, alpha)

type(SBM.case)
SBM.alpha

type(SBM.Lam)
type(SBM.D)
dic=SBM.dic
dic[7,8]

dic10=SBM.multiply(1)
dic[7,8]


D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q

dir(SBM)
############################# compute factors 
a=1
factors =mrf_tt_all.model_org(n, alpha, Lam, a)

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
Ctest
app_classes(unairesf, n, q)


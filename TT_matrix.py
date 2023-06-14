#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 09:28:10 2022

@author: abouabdallah
"""

import tt
import mrf_tt_all
import numpy as np
import time

import bases
import model_gen

import cores_melting
import compressed_mrf_tt


import symbolic_tt_mrf

alpha=[1/3,1/3,1/3]
case="assortative"
Lam=model_gen.thelam(case)
n=10
q=3
C,D= model_gen.mod (n,q, Lam, alpha)   


SBM=model_gen.model(D, Lam,alpha)

D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q


D1= np.random.randint(10,size=(10, 10))
D_rand= D1+ np.transpose(D1)

for k in range (n) :
    D_rand[k,k]=0
print(D_rand)
    

############################# compute factors 
a=1


case=mrf_tt_all.which_case(Lam)
factors =mrf_tt_all.model_org(n, alpha, Lam, a)

###################### compute z without melting

Ais=mrf_tt_all.stock_all_Ai_q(factors, q, n)


#Ais[0][1].full()

#Ais[0][3].full().shape



Ai=Ais[0][3]

Ai.tt
Ai.m
Ai.n
Ai.tt.r


import pandas as pd 

data = [['tom', 10], ['nick', 15], ['juli', 14]]
m=int(SBM.m)

data[0][0]
type(data)

def how(n,m,i):
    data=bases.listoflist(m)
    for k in range (len (Ais[0][i].to_list(Ais[0][i]))) :
        #print([list(SBM.dic.keys())[k], Ais[0][i].to_list(Ais[0][i])[k].shape ])
        data[k]=[list(SBM.dic.keys())[k], Ais[0][i].to_list(Ais[0][i])[k].size ]
    return data


i=6
def diff_num(n,m,i) :
    data =how(n,m,i)
    df = pd.DataFrame(data, columns = ['Link', 'size'])
    lol=list(df['size'] )
    numbers=np.zeros(3)
    for k in range (len(lol)) :
        if (lol[k] == 1) :
            numbers[0]+= 1
        else :
            if (lol[k] == 3) :
                numbers[1]+= 1
            else :
                numbers[2]+= 1
    return numbers

n
for k in range (n) :
    print(diff_num(n,m,k))


        





for k in range (len (Ais[0][3].to_list(Ais[0][3]))) :
    print(Ais[0][4].to_list(Ais[0][4])[k].shape)



roundprec=1e-3
Bis=mrf_tt_all.compute_Bi_all(Ais, q, n, roundprec)

Bis[3].full().shape

Bis[3].tt.r
for k in range (len (Bis[3].to_list(Bis[3]))) :
    print(Bis[3].to_list(Bis[3])[k].shape)


Bis[3].to_list(Bis[3])[2]

Ais[0][3].to_list(Ais[0][3])[2]
Ais[1][3].to_list(Ais[1][3])[2]
Ais[2][3].to_list(Ais[2][3])[2]

B=Ais[0][3]+Ais[1][3]

B.to_list(B)[2].shape


B=Bis[2]*Bis[3]

Bis[2].tt.r
Bis[3].tt.r
B.tt.r


Bis[2].n
Bis[3].n
B.n


Bis[2].m
Bis[3].m
B.m



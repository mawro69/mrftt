#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 15:43:51 2022

@author: abouabdallah
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:04:04 2022

@author: abouabdallah
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 00:57:34 2022

@author: abouabdallah
"""



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 11:40:09 2022

@author: abouabdallah
"""



import numpy as np
import time

import bases
import model_gen

import cores_melting
import compressed_mrf_tt

import compute_marginales
import os 

import mrf_tt_all


import pandas as pd

import seaborn as sns

def data_to_df(uns1,uns2) :#the format vem, mcmc, TT
    n, q = uns1.shape
    df = np.zeros((n*q, 2))
    a=0
    for k in range (n) : 
        for l in range (q) :
            
            #print(a)
            df[a,:]= [uns1[k,l],uns2[k,l]]
            a+=1

#    df = pd.DataFrame(df)
    my_frame = pd.DataFrame(data={'eps': df[:,0],'No eps': df[:,1]})
    return my_frame


def marginales_df(uns1, bns1, bns2) :
    n, q = uns1.shape
    df = np.zeros((int(n*(n-1)*q*q/2), 2))
    a=0
    for i in range (n-1) : 
        for j in range (i+1,n) :
            for k in range (q) :
                for l in range (q) :
                    df[a,:]= [bns1[i,j][k,l],bns2[i,j][k,l]]
                    a+=1
    my_frame = pd.DataFrame(data={'eps': df[:,0],'No eps': df[:,1]})

    return my_frame



# test sur les fonctions de bases


alpha=[1/3,1/3,1/3]
case="coreperiphery"
n=42
print(n)
q=3
print(q)

Lam=model_gen.thelam(case)
C,D= model_gen.mod (n,q, Lam, alpha)   
SBM=model_gen.model(D, Lam, alpha)
dic10=SBM.multiply(10)
D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q
dic10=SBM.dica
a=10
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**2

print("without cores melting")

q=3
Lam=model_gen.thelam(case)
C,D= model_gen.mod (n,q, Lam, alpha)   
SBM=model_gen.model(D, Lam, alpha)
dic10=SBM.multiply(10)
D=SBM.D
alpha=SBM.alpha
Lam=SBM.Lam
n=SBM.n
q=SBM.q
dic10=SBM.dica
a=10
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
rmax=q**2

print("begin Ai")

len(factors[0][1])



Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 0, 3, 0, q, n) 

print("end Ai")

print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")

Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 27, 1)    

print("End Z")

Z=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 1000, 1)    

print("begin margins")

start_time = time.time()
dicbij=compute_marginales.graphdesBi(Bis1eq,roundprec,27,Zeq)
binas=compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, 27, Zeq)
dicub=compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti=end_time-start_time

print(["computation time", ti])
print("end margins")

ll=list(dic10.keys())
unsnoeps, binnoeps=compute_marginales.formatptfixe(unaires, binaires, ll , q)


print("begin margins")

start_time = time.time()
dicbij=compute_marginales.graphdesBi(Bis1eq,roundprec,1000,Zeq)
binas=compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, 1000, Zeq)
dicub=compute_marginales.bina_una(binas) 
unaires2=dicub['unary']
binaires2 =dicub['binary']

end_time= time.time()
ti=end_time-start_time

print(["computation time", ti])
print("end margins")


unseps, bineps=compute_marginales.formatptfixe(unaires2, binaires2, ll , q)



df=data_to_df(unseps,unsnoeps)

df.to_csv('/home/abouabdallah/Bureau/graphiques_manuscrit/unaries.csv')

df2=marginales_df(unseps,bineps, binnoeps)

df2.to_csv('/home/abouabdallah/Bureau/graphiques_manuscrit/binaries.csv')

swarm_plot = sns.pairplot(df)

swarm_plot = sns.pairplot(df2)



print("with volume cores melting")

print("begin Ai")

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 2, 0, 300, q, n)  

print("end Ai")

Ais1eq[0]
print("begin Bi")

Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("end Bi")

print("begin Z")

Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    
print("End Z")


print("begin margins")

start_time = time.time()
dicbij=compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq)
binas=compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, 6, Zeq)
dicub=compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
end_time= time.time()
ti2=end_time-start_time
print(["computation time", ti2])
print("end margins")


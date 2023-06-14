#!/usr/bin/env python3


# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 12:21:53 2023

@author: abouabdallah
"""


import sys
import comp_mar_mcmc
import pt_fixe
import mrf_tt_all
import numpy as np
#import tt
import time
import pandas as pd
import visualisation

import label_switch



def other_format(dicub) :
    dicused=dicub['binary']
    ll=list(dicused[0,0].keys())
    dicti_bi_format={}
    q=max(max(dicused.keys()))+1
    for k in range (len(ll)) :
        dicti_bi_format[ll[k]]=np.zeros((q,q))
        for zi in range (q) :
            for zj in range (q) :
                dicti_bi_format[ll[k]][zi, zj]=abs(dicused[zj,zi][ll[k]])
                if (dicti_bi_format[ll[k]][zi, zj] < 1e-20) :
                    dicti_bi_format[ll[k]][zi, zj] =0
        
        av=np.sum(dicti_bi_format[ll[k]]) 
        if (av != 0 ) :
            dicti_bi_format[ll[k]]=  dicti_bi_format[ll[k]]*(1./ np.sum(dicti_bi_format[ll[k]]))  
    return dicti_bi_format






def ij_inclasse (ll, binaires1,binaries, res2,q) :
    m=len(ll)
    matinclasse =np.zeros((m, 3))
    for k in range (m) :
        lop=[np.sum(np.diag(binaires1[ll[k]])), np.sum(np.diag(binaries[ll[k]])),np.sum(np.diag(res2[ll[k]]))]
        matinclasse[k, :]=lop
    data = {"i,j" : ll,
            "TT": matinclasse[:, 0],
            "Gibbs": matinclasse[:, 1],
            "Pf" : matinclasse[:, 2]}
    df3 = pd.DataFrame(data)

    return df3


def ij_inclassettgibs (ll, binaires1,binaries,q) :
    m=len(ll)
    matinclasse =np.zeros((m, 2))
    for k in range (m) :
        lop=[np.sum(np.diag(binaires1[ll[k]])), np.sum(np.diag(binaries[ll[k]]))]
        matinclasse[k, :]=lop
    data = {"i,j" : ll,
            "TT": matinclasse[:, 0],
            "Gibbs": matinclasse[:, 1]}
    df3 = pd.DataFrame(data)

    return df3



ff   = sys.argv[1]
case = sys.argv[2]
n = sys.argv[3]

n=int(n)



print("begin compute")

print(case)
print(ff)
print(n)
rs1=0
alpha=[1/3,1/3,1/3]
#case="coreperiphery"



print(["experiernce no", ff, case, n])
q=3
a=10
if (case =='hierarchical') :
    a=10**(4/5)

Lam=mrf_tt_all.model_gen.thelam(case)
dd=0

C,D= mrf_tt_all.model_gen.mod (n,q, Lam, alpha)   
SBM=mrf_tt_all.model_gen.model(D, Lam, alpha)

np.savetxt('D50nn' +str(n)+ case + str(ff) +'.txt',D,fmt='%.2f')
np.savetxt('C50nn'+ str(n)+case + str(ff)+'.txt',C,fmt='%.2f')

print("distances matrix and classes saved")

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



tempsdecalcul=[]
print("start pt fixe")

dicu=pt_fixe.unary_factors (n,q, D, alpha)            

rs=pt_fixe.pt_fixe (alpha, D, Lam, n, q, epsilone, dicu )
rs1=np.linalg.norm(rs)
print(rs1)

if (rs1 ==0) :
    rs2=pt_fixe.pt_fixe_allit (alpha, D, Lam, n, q, epsilone, dicu )
    print("citeration saved and no conergence")
    #np.savetxt('respf'+ case + str(ff) +'.txt',rs2,fmt='%.2f')
    np.save('respf'+ case + str(ff)+'.npy', np.array(rs2, dtype=object), allow_pickle=True)
else :
    start_time = time.time()

    res, res2=pt_fixe.unaires_binaires_ptfixe(alpha, D, Lam, n, q)
    print("result obtained and conergence")
    end_time = time.time()
    tempsdecalcul.append(end_time- start_time)
    print(tempsdecalcul)



#b = np.load('respf'+ case + str(ff)+'.npy', allow_pickle=True)


print("end pt fixe")
print("start TT")





ll=list(dic10.keys())
factors =mrf_tt_all.model_org(n, alpha, Lam, a)
roundprec=5e-2
#rmax=q**3
rmax=2*(q**3)


print("Compute Ai and Bi")

Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)

print("rounding parameters ")

print(["rounding prec",roundprec,"rmax", rmax])



print(["comute z","rmax = 3 9 27 " + str(rmax)])


start_time = time.time()
Zeq4=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    
end_time= time.time()


start_time = time.time()
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 27, 1)    
end_time= time.time()
ti=end_time-start_time
Zeq2=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 9, 1)    
Zeq3=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, 3, 1)    



print(Zeq3, Zeq2, Zeq, Zeq4)




print("start comp marginales")


print("comute the graph")



Ais1eq=mrf_tt_all.compute_Ai_per_method(dic10, 1, 2, None, q, n)
Bis1eq=mrf_tt_all.compute_Bi_all(Ais1eq, q, n, roundprec)
Zeq=mrf_tt_all.compute_W_simple(Bis1eq, roundprec, rmax, 1)    

start_time = time.time()
dicbij=mrf_tt_all.compute_marginales.graphdesBi(Bis1eq,roundprec,rmax,Zeq4)
binas=mrf_tt_all.compute_marginales.all_binaries_all_states (dicbij,Ais1eq, q,roundprec, rmax, Zeq4)
dicub=mrf_tt_all.compute_marginales.bina_una(binas) 
unaires=dicub['unary']
binaires =dicub['binary']
ti=end_time-start_time
uns, bins=mrf_tt_all.compute_marginales.formatptfixe(unaires, binaires, ll , q)
print("end TT")
end_time = time.time()
tempsdecalcul.append(end_time- start_time)
print(tempsdecalcul)



#print(uns)

#df=visualisation.data_to_df(uns,res, "TT", "Pf")
#    df.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/coreperiphery/Dfnn60'+case + str(ff)+'.csv')
#df2=visualisation.marginales_df(uns,bins,res2, "TT", "Pf")
#    df2.to_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/coreperiphery/Dfnn60'+case + str(ff)+'bin.csv')
#swarm_plot=visualisation.swarm(uns,res, "TT", "Pf", 1)

print("start gibbs marginales")
start_time = time.time()
dicb=comp_mar_mcmc.all_fac (dicb)
dicu= comp_mar_mcmc.dicu_tomat(dicu,n,q)

Niter=210000

firstone= 10000

latence = 200
dicmcmc=comp_mar_mcmc.echantillons_token(n, q, dicu, dicb, Niter, firstone, latence)
np.array(C)
binaries= comp_mar_mcmc.margbinaires(dicmcmc)
unaries1 = comp_mar_mcmc.margunaires(dicmcmc)
end_time = time.time()
tempsdecalcul.append(end_time- start_time)
print(tempsdecalcul)

print("end gibbs marginales")

#df3=visualisation.data_to_df2(uns,res,unaries1, "TT", "Pf", "Gibbs")

if (np.linalg.norm(rs)==0) :
    df=visualisation.data_to_df(uns, unaries1, "TT", "Gibbs")
    df.to_csv('Dfnn60'+case + str(ff)+str(n)+'.csv')
    df2=visualisation.marginales_df(uns,bins,binaries, "TT", "Gibbs")
    df2.to_csv('Dfnn60'+case + str(ff)+str(n)+'bin.csv')
    binaires1 =other_format(dicub)
    df3 =ij_inclassettgibs (ll, binaires1,binaries,q)
    df3.to_csv('ij_inclasse'+case +'_' + str(ff)+'_'+str(n)+'bin.csv')

else :
    dfall=visualisation.data_to_df2(uns,res, unaries1, "TT", "Pf", "Gibbs")
    dfall.to_csv('Dfnn60'+case + str(ff)+str(n)+'.csv')
    df2all=visualisation.marginales_df2(uns,bins,res2, binaries, "TT", "Pf", "Gibbs")
    df2all.to_csv('Dfnn60'+case + str(ff)+str(n)+'bin.csv')
    binaires1 =other_format(dicub)
    df3 =ij_inclasse (ll, binaires1,binaries, res2,q)
    df3.to_csv('ij_inclasse'+case +'_' + str(ff)+'_'+str(n)+'bin.csv')



print(df3.head())


#binaries #gibbs 
#binaires #tt
#res2




np.savetxt('tempsdecalcul'+ case + str(ff) +'.txt',tempsdecalcul,fmt='%.2f')






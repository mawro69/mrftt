#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:16:15 2022

@author: abouabdallah
"""

import tt 
import numpy as np

import mrf_tt_all

def listoflist(d):
    """ initialise a liste of liste of shape d

   arguments :
       @ d  int : shape of litse

   return :
       @ listeflists : lempty list of list"""

    listoflists = []
    a_list = [1, 5]
    for k in range(d):
        listoflists.append(a_list)
    return(listoflists)




def get_m(n):
    return int(n*(n-1)/2)



def getFactors(n):
    # Create an empty list for factors
    factors=[];

    # Loop over all factors
    for i in range(1, n + 1):
        if n % i == 0:
            factors.append(i)

    # Return the list of factors
    return factors


def getFactors_clean(n) :
    lis=getFactors(n)
    f=len(lis)
    lis.remove(lis[f-1])  
    lis.remove(lis[0])
    return lis

def how_comp(m,k) :
    lis=getFactors(m)
    if (k in lis) :
        numcomp=int(m/k)
        comp=k*np.ones(numcomp)
    else : 
        return("error")
    return (numcomp,comp)


def to_compress(comp, numcomp):
    tobecomp=[]
    tobecomp=listoflist(numcomp)
    a=0
    for k in range (numcomp) :
        tobecomp[k]=[]
        for l in range (int(comp[k])) :
            tobecomp[k].append(a)
            a+=1
    return tobecomp

def compressed_cores(factorCoreArray, di, i, futurcoresAi, comp, numcomp) :
    tobecomp=to_compress(comp, numcomp)
    facs=factorCoreArray[i]
    for k in range (numcomp) :
        n1,n2=facs[tobecomp[k][0]][:,di, :].shape
        futurcoresAi[k]=np.reshape(facs[tobecomp[k][0]][:,di, :], (1,n1,n2,1))
        for j in range (1,len(tobecomp[k])) :
            n1,n2=facs[tobecomp[k][j]][:,di, :].shape
            a=np.reshape(facs[tobecomp[k][j]][:,di, :], (1,n1,n2,1))
            futurcoresAi[k]=np.kron(futurcoresAi[k],a)
    return futurcoresAi




def stock_Ai_compressed( factors, i, di,comp, numcomp): 
    """calcul des matrices  Ai = Kron (G_i^l)
    arguments : 
        factorCoreArray list of list: Une liste de liste ou chaque factorCoreArray[i] nous donne 
        les éléments qui constituent la tt matrice Ai   
    return :
        Bi matrix : """
    factorRankArray, factorModeSizeArray , factorCoreArray=mrf_tt_all.init_Bi_Ai(factors) 
    m=len(factorCoreArray[0])
    length=int(m/numcomp)
    rowRankArray = factorRankArray[i, :]#stock les différents r_i jouant les numéros de lignes 
    colRankArray = factorRankArray[i+ 1, :]# %stock les différents r_i jouant les numéros de colonnes 
    rowRankArraycomp = np.ones(numcomp)
    colRankArraycomp = np.ones(numcomp) 
    res=to_compress(comp, numcomp)
    for k in range (numcomp) :
        for j in range (len(res[k])) :
            rowRankArraycomp[k]=rowRankArraycomp[k]*rowRankArray[res[k][j]]
            colRankArraycomp[k]=colRankArraycomp[k]*colRankArray[res[k][j]]
    Ai=tt.ones(mrf_tt_all.multi(rowRankArraycomp,colRankArraycomp))
    Ai=tt.matrix(Ai,n=rowRankArraycomp, m=colRankArraycomp)
    futurcoresAi=tt.matrix.to_list(Ai)
    futurcoresAi=compressed_cores(factorCoreArray, di, i, futurcoresAi, comp, numcomp)
    Ai=Ai.from_list(futurcoresAi)

            
    return Ai



def stock_all_Ai_compressed(factors,di,n,comp, numcomp):
    listofAi=[]
    for i in range (n):
        listofAi.append(stock_Ai_compressed(factors, i, di,comp, numcomp))
    return listofAi


def stock_all_Ai_q_compressed(factors,q,n,comp, numcomp):
    listofAiq=[]
    for di in range (q ):
        listofAiq.append(stock_all_Ai_compressed(factors,di,n,comp, numcomp))
    
    return listofAiq


def compute_Bi_f_compressed(factors, q, n,roundprec,comp, numcomp) :
    Bi=[]
    listofAiq=stock_all_Ai_q_compressed(factors,q,n,comp, numcomp)
    i=0
    for i in range (n) :
        print(i)
        rowRankArray=listofAiq[0][i].n
        colRankArray=listofAiq[0][i].m
        a=mrf_tt_all.tt_zeros(mrf_tt_all.multi(rowRankArray,colRankArray))
        a=listofAiq[0][i]
#        a=a.round(roundprec)
        for di in range (1,q) :
            #print(di)
            a=a+listofAiq[di][i]
            #.round(roundprec)
#            a=a.round(roundprec)
        Bi.append(a)
    return Bi


def compute_z( q, n, roundprec,listeoffullBi) :
    """calcul de la constante de normalisation z 
    arguments : 
        ttfact  list of list: liste des facteurs gi du model 
    return :
        z real : constante de normalisation"""
    #listeoffullBi=compute_Bi(factors, q, n,roundprec)
    z=listeoffullBi[n-1]
    z=z.round(roundprec)
    for k in range (n-2,0, -1) :
        #print(listeoffullBi[k].shape)
        print(k)
        z=(listeoffullBi[k].round( roundprec))*z
        z=z.round( roundprec)
        print(z.erank)


    z=(listeoffullBi[0].round(roundprec))*z
    return z




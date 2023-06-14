#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:21:59 2022

@author: abouabdallah
"""

import numpy as np 
import tt 

import mrf_tt_all



def compute_log_Z(Z, n) :
    val = np.log10(Z) - n*(n-1)/2
    return val

def compute_margin(Ais, Bis, i, j, n, Q, roundprec) :
    rmax= Q**2
    mar=np.zeros((Q,Q))
    listocomp=Bis
    
    for di in range (Q) :
        for dj in range (Q) :
            listocomp[i]=Ais[di][i]
            listocomp[j]=Ais[dj][j]
            #print([di,dj])
            mar[di,dj]= mrf_tt_all.compute_W_simple(listocomp, roundprec, rmax, 1)   
    return mar


def message_passing(Ais, Bis, i, j, n, Q, roundprec) :
    rcl=[]
    rmax= Q**2
    listocomp=Bis
    
    
    W=1
    for k in range (i-1,-1, -1) :
        print(k)
        Bi=listocomp[k]
        W=W*Bi
        W=W.round(roundprec, rmax)
    
    rcl.append(W)
    print(i)
    W=1
    for k in range (j-1,i, -1) :
        print(k)
        Bi=listocomp[k]
        W=W*Bi
        W=W.round(roundprec, rmax)
    
    rcl.append(W)
    print(j)

    W=1
    for k in range (n-1, j, -1) :
        print(k)
        Bi=listocomp[k]
        W=Bi*W
        W=W.round(roundprec, rmax)
    
    rcl.append(W)
        

    return rcl



def leftside(Ais, Bis, i, j, n, Q, roundprec) :
    left=1
    rmax=Q**2
    for k in range (n-1, j-1, -1) :
        print(k)
        Bi=Bis[k]
        left=Bi*left
        left=left.round(roundprec, rmax)

    
    return left


def compute_on_margin_smart(rcl, Ais, Bis, i, j,zi,zj,  n, rmax, roundprec ) :
    rcl2=[]
    rcl2.append(rcl[0])
    rcl2.append(Ais[zi][i])
    rcl2.append(rcl[1])
    rcl2.append(Ais[zj][j])
    rcl2.append(rcl[2])
    themar=mrf_tt_all.compute_W_simple(rcl2, roundprec, rmax, 1)

    return themar


def compute_on_margin_smart2(rcl, Ais, Bis, i, j,zi,zj,  n, rmax, roundprec ) :
    rcl2=rcl[2]
    rcl2=rcl2*Ais[zj][j]
    rcl2=rcl2.round(roundprec, rmax)
    for k in range (j-1,i, -1) :
#        print(k)
        Bi=Bis[k]
        rcl2=rcl2*Bi
        rcl2=rcl2.round(roundprec, rmax)
    
    rcl2=rcl2*Ais[zi][i]
    rcl2=rcl2.round(roundprec, rmax)
    rcl2=rcl2*rcl[0]
    rcl2=rcl2.round(roundprec, rmax)

    
    themar=mrf_tt_all.to_float(rcl2)




    return themar




def compute_margin_smart(Ais, Bis, i, j, n, Q, roundprec, rmax) :
    rmax= Q**2
    mar=np.zeros((Q,Q))
    rcl=message_passing(Ais, Bis, i, j, n, Q, roundprec)

    
    for di in range (Q) :
        for dj in range (Q) :
            print([di,dj])
            mar[di,dj]= compute_on_margin_smart(rcl, Ais, Bis, i, j,di,dj,  n, rmax, roundprec ) 
    return mar



###############################################################message passing 


def graphdesBi(Bis,roundprec,rmax,w):
     #
     dic_Cij = {} # dictionnaire vide
     n = len(Bis) #n nombre d'individus
     dic_Cij[0,n-1]=w #format tt
     for k in range (n) :
         dic_Cij[(k,k)] = Bis[k] # on crée la base de notre
     for i in range(n):
         #
         print(i)
         for j in range (i+1,n):
             if ([i,j] != [0,n-1]) : 
                 
                 C_ij =  (dic_Cij[(i,j-1)].round(roundprec, rmax)) * (dic_Cij[(j,j)].round(roundprec, rmax))
                 C_ij =C_ij.round(roundprec,rmax)        
                 dic_Cij[(i,j)] = C_ij

     #
     return dic_Cij

def dictofleftprimes(dicti,Ais, di,roundprec, rmax ) :#modifié format alain
    n=len(Ais[di])
    i=0
    dictleft={}
    dictleft[i]=np.reshape(1,(1,1))
    dictleft[i+1]=dicti[i,i]
    for k in range (1,n) :
        dictleft[k+1]= dicti[i,k]
    dictleftprime={}
    dictleftprime[i]=Ais[di][i].round(roundprec, rmax)
    for k in range (1,n) :
        dictleftprime[k]=(dictleft[k].round(roundprec, rmax))*(Ais[di][k].round(roundprec, rmax))
        dictleftprime[k]=dictleftprime[k].round(roundprec, rmax)
    return(dictleftprime,dictleft)








def dictofrightprimes(dicti,Ais, di,roundprec, rmax ) :#modifié format alain
    n=len(Ais[di])
    i=n-1
    dictright={}
    for k in range (0,n-2) :
        dictright[k]= dicti[k+1,i]
    dictright[i-1]=dicti[i,i]
    dictright[i]=1
    dictrightprime={}
    for k in range (0,n) :
        
        dictrightprime[k]=(Ais[di][k].round(roundprec, rmax))* dictright[k]
        dictrightprime[k]=dictrightprime[k].round(roundprec, rmax)

        

    return (dictrightprime,dictright)






def compute_binaries (dic_Cij,Ais, di1,di2,roundprec,rmax,w):
    n=len(Ais[di1])
    dictofrightprime,dictright=dictofrightprimes(dic_Cij,Ais, di1,roundprec, rmax )
    dictofleftprime,dictleft=dictofleftprimes(dic_Cij,Ais, di2,roundprec, rmax )
    dicomarg={}
    z=w
    div=1./z
    for i in range (n-1) :
        dicomarg[i,i+1]=(dictofleftprime[i].round(roundprec, rmax))*(dictofrightprime[i+1].round(roundprec, rmax))
        aa=div*dicomarg[i,i+1].round(roundprec, rmax)
        dicomarg[i,i+1]=mrf_tt_all.to_float(aa)
        dicomarg[i+1,i]=dicomarg[i,i+1]
        #print("srep1")
    for i in range (n) :
        temp=dictofleftprime[i]*dictright[i]*div
        dicomarg[i,i]=mrf_tt_all.to_float(temp)
        #print("srep2")
    for i in range (n-2) :
        for j in range (i+2,n) :
            a=(dictofleftprime[i].round(roundprec, rmax))*(dic_Cij[i+1,j-1].round(roundprec, rmax))
            a=a.round(roundprec, rmax)
            b=a*(dictofrightprime[j].round(roundprec, rmax))
            b=b*div            
            dicomarg[i,j]=mrf_tt_all.to_float(b)
            dicomarg[j,i]=dicomarg[i,j]
            #print("srep3")
    return (dicomarg)


def all_binaries_all_states (dic_Cij,Ais, Q,roundprec, rmax, w) :
    thedic ={}
    print("it begins")
    for di1 in range (Q) :
        for di2 in range (Q) :
            print([di1, di2])
            thedic[di1, di2] = compute_binaries (dic_Cij,Ais, di1,di2,roundprec,rmax,w)
    
    return thedic


def compute_binaries2 (dic_Cij,Ais, di1,di2,roundprec,rmax):
    n=len(Ais[di1])
    dictofrightprime,dictright=dictofrightprimes(dic_Cij,Ais, di1,roundprec,rmax)
    dictofleftprime,dictleft=dictofleftprimes(dic_Cij,Ais, di2, roundprec,rmax)
    dicomarg={}
    z=mrf_tt_all.to_float(dic_Cij[0,n-1])
    div=1./z
    for i in range (n-1) :
        dicomarg[i,i+1]=dictofleftprime[i]*dictofrightprime[i+1]
        aa=div*dicomarg[i,i+1]
        dicomarg[i,i+1]=mrf_tt_all.to_float(aa)
        dicomarg[i+1,i]=dicomarg[i,i+1]
        
        
    for i in range (n) :
        temp=dictofleftprime[i]*dictright[i]*div
        dicomarg[i,i]=mrf_tt_all.to_float(temp)
    for i in range (n-2) :
        for j in range (i+2,n) :
            a=dictofleftprime[i]*dic_Cij[i+1,j-1]
            a=a.round(roundprec)
            b=a*dictofrightprime[j]
            b=b*div            
            dicomarg[i,j]=mrf_tt_all.to_float(b)
            #dicomarg[j,i]=dicomarg[i,j] under check
    return (dicomarg)# for check



def all_binaries_all_states2 (dic_Cij,Ais, q,roundprec, rmax) :
    thedic ={}
    for di1 in range (q) :
        for di2 in range (q) :
            thedic[di1, di2] = compute_binaries2 (dic_Cij,Ais, di1,di2,roundprec,rmax)

    
    return thedic





def bina_una(thedic) :
    n=max(max(thedic[0,0].keys()))+1
    dicm ={}
    dicu={}
    B=max(max(thedic.keys()))+1
    for di1 in range (B) :
        dicu[di1] = {}
        for i in range (n) :
            dicu[di1][i] = thedic[di1, di1][i,i]
        for di2 in range (B) :
            dicm[di1,di2] = {}
            for i in range (n-1) :
                for j in range (i,n) :
                    dicm[di1,di2][i,j]=thedic[di1, di2][i,j]
    
    
    dicub={}         
    dicub['unary']=dicu
    dicub['binary']=dicm
    
    return (dicub)
    


"""dicub=mrf_tt.bina_una(thedic)
    unaires=dicub['unary']
    binaires =dicub['binary']
"""



def bianries_update(binaires,n,q) :
    aa=list(binaires.keys())
    binaries_updates=binaires
    for di in range (q) :
        for dj in range (q):
            for i in range (n-1):
                for j in range (i+1,n) :
                    binaries_updates[dj,di][j,i]=binaires[di,dj][i,j]
        
    return (binaries_updates)


def bina_una2(thedic) :
    n=max(max(thedic[0,0].keys()))+1
    dicm ={}
    dicu={}
    B=max(max(thedic.keys()))+1
    for di1 in range (B) :
        dicu[di1] = {}
        for i in range (n) :
            dicu[di1][i] = thedic[di1, di1][i,i]
        for di2 in range (B) :
            dicm[di1,di2] = {}
            for i in range (n-1) :
                for j in range (i+1,n) :
                    dicm[di1,di2][i,j]=thedic[di1, di2][i,j]
    
    
    dicub={}         
    dicub['unary']=dicu
    dicub['binary']=dicm
    
    return (dicub)



########normalisations 



def formatptfixe(unaires, binaires, ll , q) : 
    n=max(list(unaires[0].keys()))+1
    restt=np.ones((n,q))   
    for k in range (n) :
        for l in range (3) :
            restt[k,l]=abs(unaires[l][k]) #beug lie à la reprsentation vigue flottante -0.0 pour une tres faible valeur
    restt=normalisation(restt)
    res2tt=dict()
    for k in range (len(ll)) :
        res2tt[ll[k]]=np.zeros((q,q))
        for i in range (q) :
            for j in range (q) :
                res2tt[ll[k]][i,j]=abs(binaires[(i,j)][ll[k][0],ll[k][1]]) #beug lie à la reprsentation vigue flottante -0.0 pour une tres faible valeur
        mat=res2tt[ll[k]]
        aaa=np.sum(mat)
        if (aaa > 0) :
            mat=mat/aaa
        res2tt[ll[k]]=mat




    return restt, res2tt



def formatptfixe2(unaires, binaires, ll , q) : 
    n=max(list(unaires[0].keys()))+1
    restt=np.ones((n,q))   
    for k in range (n) :
        for l in range (3) :
            restt[k,l]=abs(unaires[l][k]) #beug lie à la reprsentation vigue flottante -0.0 pour une tres faible valeur
    restt=normalisation(restt)
    res2tt=dict()
    for k in range (len(ll)) :
        res2tt[ll[k]]=np.zeros((q,q))
        for i in range (q) :
            for j in range (q) :
                res2tt[ll[k]][i,j]=abs(binaires[(i,j)][ll[k][0],ll[k][1]]) #beug lie à la reprsentation vigue flottante -0.0 pour une tres faible valeur
        mat=res2tt[ll[k]]
        mat=mat/np.sum(mat)
        res2tt[ll[k]]=mat




    return restt, res2tt

def normalisation(res) :
    n,q=res.shape
    for k in range (n) :
        for l in range (q):
            aaa=np.sum(np.sum(res[k,:]))
            if (aaa > 0) :
                res[k,l]=res[k,l]/aaa
            if (res[k,l]<1e-10):
                res[k,l]=0
    return res



def normalisation2(res) :
    n,q=res.shape
    for k in range (n) :
        for l in range (q):
            res[k,l]=res[k,l]/np.sum(res[k,:])
            if (res[k,l]<1e-10):
                res[k,l]=0
    return res


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
        dicti_bi_format[ll[k]]=  dicti_bi_format[ll[k]]*(1./ np.sum(dicti_bi_format[ll[k]]))  
    return dicti_bi_format



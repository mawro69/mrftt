#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 15:21:55 2023

@author: abouabdallah
"""


import numpy as np 
import pandas as pd 
from scipy.stats import pearsonr
import matplotlib.pyplot as plt


cases=[]
sdspf=[]
mnspf=[]

mnstt=[]
sdtt=[]

case="assortative"
n=9
#print(n)
q=3
#print(q)



llte=[]
llpe=[]
llte2=[]
llpe2=[]

for num in range (100) :
    #print(num) 
    df1=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    TT=list(df1['TT'])
    Pf= list(df1['Pf'])
    Exact= list(df1['Exact'])
    corr, _ = pearsonr(TT, Exact)
    llte.append(corr)
    corr, _ = pearsonr(Pf, Exact)
    llpe.append(corr)
    TT2=list(df2['TT'])
    Pf2= list(df2['Pf'])
    Exact2= list(df2['Gibbs'])
    corr, _ = pearsonr(TT2, Exact2)
    llte2.append(corr)
    corr, _ = pearsonr(Pf2, Exact2)
    llpe2.append(corr)


print(case)

print(np.mean(llte2))
print(np.std(llte2))    

print(np.mean(llpe2))
print(np.std(llpe2))    

cases.append(case)
sdspf.append(np.std(llpe2))
mnspf.append(np.mean(llpe2))

mnstt.append(np.mean(llte2))
sdtt.append(np.std(llte2))


plt.hist(llpe2)
plt.title('Histo cor entre pt fixe et exact  cas ' + case)
plt.hist(llte2)
plt.title('Histo cor entre TT et exact cas ' + case)


case="dissortative"
n=9
#print(n)
q=3
#print(q)

llte=[]
llpe=[]
llte2=[]
llpe2=[]

for num in range (100) :
    ##print(num) 
    df1=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    TT=list(df1['TT'])
    Pf= list(df1['Pf'])
    Exact= list(df1['Exact'])
    corr, _ = pearsonr(TT, Exact)
    llte.append(corr)
    corr, _ = pearsonr(Pf, Exact)
    llpe.append(corr)
    TT2=list(df2['TT'])
    Pf2= list(df2['Pf'])
    Exact2= list(df2['Gibbs'])
    corr, _ = pearsonr(TT2, Exact2)
    llte2.append(corr)
    corr, _ = pearsonr(Pf2, Exact2)
    llpe2.append(corr)



print(case)
print(np.mean(llte2))
print(np.std(llte2))    

print(np.mean(llpe2))
print(np.std(llpe2))    

cases.append(case)
sdspf.append(np.std(llpe2))
mnspf.append(np.mean(llpe2))

mnstt.append(np.mean(llte2))
sdtt.append(np.std(llte2))


plt.hist(llpe2)
plt.title('Histo cor entre pt fixe et exact  cas ' + case)
plt.hist(llte2)
plt.title('Histo cor entre TT et exact cas ' + case)

case="coreperiphery"
n=9
#print(n)
q=3
#print(q)

llte=[]
llpe=[]
llte2=[]
llpe2=[]

for num in range (100) :
    #print(num) 
    df1=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    TT=list(df1['TT'])
    Pf= list(df1['Pf'])
    Exact= list(df1['Exact'])
    corr, _ = pearsonr(TT, Exact)
    llte.append(corr)
    corr, _ = pearsonr(Pf, Exact)
    llpe.append(corr)
    TT2=list(df2['TT'])
    Pf2= list(df2['Pf'])
    Exact2= list(df2['Gibbs'])
    corr, _ = pearsonr(TT2, Exact2)
    llte2.append(corr)
    corr, _ = pearsonr(Pf2, Exact2)
    llpe2.append(corr)



print(case)

print(np.mean(llte2))
print(np.std(llte2))    

print(np.mean(llpe2))
print(np.std(llpe2))    

cases.append(case)
sdspf.append(np.std(llpe2))
mnspf.append(np.mean(llpe2))

mnstt.append(np.mean(llte2))
sdtt.append(np.std(llte2))

plt.hist(llpe2)
plt.title('Histo cor entre pt fixe et exact  cas ' + case)
plt.hist(llte2)
plt.title('Histo cor entre TT et exact cas ' + case)



case="hierarchical"
n=9
#print(n)
q=3
#print(q)

llte=[]
llpe=[]
llte2=[]
llpe2=[]

for num in range (100) :
    #print(num) 
    df1=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    TT=list(df1['TT'])
    Pf= list(df1['Pf'])
    Exact= list(df1['Exact'])
    corr, _ = pearsonr(TT, Exact)
    llte.append(corr)
    corr, _ = pearsonr(Pf, Exact)
    llpe.append(corr)
    TT2=list(df2['TT'])
    Pf2= list(df2['Pf'])
    Exact2= list(df2['Gibbs'])
    corr, _ = pearsonr(TT2, Exact2)
    llte2.append(corr)
    corr, _ = pearsonr(Pf2, Exact2)
    llpe2.append(corr)




print(case)

print(np.mean(llte2))
print(np.std(llte2))    

print(np.mean(llpe2))
print(np.std(llpe2))    

cases.append(case)
sdspf.append(np.std(llpe2))
mnspf.append(np.mean(llpe2))
mnstt.append(np.mean(llte2))
sdtt.append(np.std(llte2))

plt.hist(llpe2)
plt.title('Histo cor entre pt fixe et exact  cas ' + case)
plt.hist(llte2)
plt.title('Histo cor entre TT et exact cas ' + case)


case="hierarchical2"
n=9
#print(n)
q=3
#print(q)

llte=[]
llpe=[]
llte2=[]
llpe2=[]

for num in range (100) :
    #print(num) 
    df1=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    TT=list(df1['TT'])
    Pf= list(df1['Pf'])
    Exact= list(df1['Exact'])
    corr, _ = pearsonr(TT, Exact)
    llte.append(corr)
    corr, _ = pearsonr(Pf, Exact)
    llpe.append(corr)
    TT2=list(df2['TT'])
    Pf2= list(df2['Pf'])
    Exact2= list(df2['Gibbs'])
    corr, _ = pearsonr(TT2, Exact2)
    llte2.append(corr)
    corr, _ = pearsonr(Pf2, Exact2)
    llpe2.append(corr)


print(case)

print(np.mean(llte2))
print(np.std(llte2))    

print(np.mean(llpe2))
print(np.std(llpe2))    

cases.append(case)
sdspf.append(np.std(llpe2))
mnspf.append(np.mean(llpe2))

mnstt.append(np.mean(llte2))
sdtt.append(np.std(llte2))

plt.hist(llpe2)
plt.title('Histo cor entre pt fixe et exact  cas ' + case)
plt.hist(llte2)
plt.title('Histo cor entre TT et exact cas ' + case)

case="ordered"
n=9
#print(n)
q=3
#print(q)

llte=[]
llpe=[]
llte2=[]
llpe2=[]

for num in range (100) :
    #print(num) 
    df1=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'.csv')
    df2=pd.read_csv('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/Dfnn9_'+ str(num) + '_'+ case +'bin.csv')
    TT=list(df1['TT'])
    Pf= list(df1['Pf'])
    Exact= list(df1['Exact'])
    corr, _ = pearsonr(TT, Exact)
    llte.append(corr)
    corr, _ = pearsonr(Pf, Exact)
    llpe.append(corr)
    TT2=list(df2['TT'])
    Pf2= list(df2['Pf'])
    Exact2= list(df2['Gibbs'])
    corr, _ = pearsonr(TT2, Exact2)
    llte2.append(corr)
    corr, _ = pearsonr(Pf2, Exact2)
    llpe2.append(corr)


print(case)

print(np.mean(llte2))
print(np.std(llte2))    

print(np.mean(llpe2))
print(np.std(llpe2))


cases.append(case)
sdspf.append(np.std(llpe2))
mnspf.append(np.mean(llpe2))

mnstt.append(np.mean(llte2))
sdtt.append(np.std(llte2))

plt.hist(llpe2)
plt.title('Histo cor entre pt fixe et exact  cas ' + case)
plt.hist(llte2)
plt.title('Histo cor entre TT et exact cas ' + case)


cases


for k in range (6) :
    print(k)
    sdspf[k]= round(sdspf[k], 2)
    mnspf[k]= round(mnspf[k], 2)
    mnstt[k]= round(mnstt[k], 2)
    sdtt[k]= round(sdtt[k], 2)



cases
mnspf
sdspf
mnstt
sdtt


cases
timeexact=[]
timett=[]
timepf=[]
k=5
timescasespf=[]
timescasestt=[]
timescasesexact=[]
for k in range (6) :
    alltimes=np.loadtxt('/home/abouabdallah/Bureau/docurespres/article/resultats/n9/times'+ '_'+ cases[k] +'.txt')
    timeexact=[]
    timett=[]
    timepf=[]

    for num in range (100) :
        timeexact.append(alltimes[num][0])
        timett.append(alltimes[num][1])
        timepf.append(alltimes[num][2])
    timescasespf.append([cases[k],round(np.mean(timepf),2),round(np.std(timepf) ,2)])
    timescasestt.append([cases[k],round(np.mean(timett),2),round(np.std(timett),2)])
    timescasesexact.append([cases[k],round(np.mean(timeexact),2),round(np.std(timeexact),2)])


        




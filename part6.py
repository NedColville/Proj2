# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 11:31:36 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
import math
test1=(123.2375145522137, 127.15848347691396)
HP=[1.4,125.1,470/(math.sqrt(2*np.pi))]

def gaus(x,sig,mean,norm):
    norm=norm/sig
    return norm*np.exp(-(x-mean)**2/(2*sig**2))

def trap(x,f,sig,mean,norm,N):
    a=x[0]
    b=x[len(x)-1]
    h=(b-a)/(N)
    val=0
    for i in range(0,N):
        val+=0.5*h*(f(a+i*h,sig,mean,norm)+f(a+(i+1)*h,sig,mean,norm))
    return val

def rNb(ml,mu):
    return math.sqrt((1500*20)*(np.exp((125.1-ml)/20)-np.exp((125.1-mu)/20)))

def main(ml,mu,HP,nList=1000,nInteg=100):
    m=np.linspace(ml,mu,nList)
    rootNB=rNb(ml,mu)
    nH=trap(m,gaus,*HP,nInteg)
    return nH/rootNB
vals=[]
ms=[]
Nhs=[]
errs=[]
for i in range(-40,40):
    dm=i/40*0.2
    HP=[1.4,125.1+dm,470/(math.sqrt(2*np.pi))]
    ms.append(dm)
    vals.append(main(*test1,HP))
    Nhs.append(trap(np.linspace(test1[0],test1[1],1000),gaus,*HP,1000))
    
plt.plot(ms,vals)
errs.append(np.std(Nhs))
plt.xlabel(r"$Change in Mass, GeV/c^{2}$")
plt.ylabel("Significance")
plt.show()
plt.plot(ms,Nhs)
plt.xlabel(r"Change in Mass, $GeV/c^{2}$")
plt.ylabel("Number of Higgs Produced")
plt.show()
vals=[]
ms=[]
Nhs=[]
for i in range(0,50):
    perc=i*0.04/50
    avM=perc*124.5+(1-perc)*125.1
    avSig=perc*2.6+(1-perc)*1.4
    HP=[avSig,avM,470/(math.sqrt(2*np.pi))]
    ms.append(perc)
    vals.append(main(*test1,HP))
    Nhs.append(trap(np.linspace(test1[0],test1[1],1000),gaus,*HP,1000))
plt.xlabel("Percentage of Affected Photons")
plt.ylabel("Significance")
plt.plot(ms,vals)
errs.append(np.std(Nhs))
plt.show()
plt.plot(ms,Nhs)
plt.xlabel("Percentage of Affected Photons")
plt.ylabel("Number of Higgs Produced")
plt.show()
Nhs=[]
ms=[]
vals=[]
for i in range(-40,40):
    dN=i/40*0.03*470
    HP=[1.4,125.1,(470+dN)/(math.sqrt(2*np.pi))]
    ms.append(dN)
    vals.append(main(*test1,HP))
    Nhs.append(trap(np.linspace(test1[0],test1[1],1000),gaus,*HP,1000))


errs.append(np.std(Nhs))
plt.plot(ms,Nhs)
plt.xlabel("Change in Expected Number of HB Created")
plt.ylabel("Number of Higgs Produced")
plt.show()
plt.plot(ms,vals)
plt.xlabel("Change in Expected Number of HB Created")
plt.ylabel("Significance")
plt.show()
plt.plot(Nhs)

print(errs)
HP=[1.4,125.1,470/(math.sqrt(2*np.pi))]
Nh=trap(np.linspace(test1[0],test1[1],1000),gaus,*HP,1000)
err=0
for i in range(0,3):
    err+=(errs[i])**2
err=err
plt.show()

def getP(ml,mu):
    Nb=(rNb(ml,mu))**2
    S=main(ml,mu,HP)
    N=trap(np.linspace(ml,mu,10000),gaus,*HP,1000)+Nb
    print(N**0.5,Nb**0.5)
    sigS=((N+err)/Nb+((N-Nb)**2)/(4*Nb**2))**0.5
    x=np.linspace(0,10,10000)
    y=gaus(x,sigS,S,1/(2*np.pi)**0.5)
    plt.plot(x,y)
    plt.xlabel("Significance")
    plt.ylabel("")
    plt.axvline(S,linestyle='dashed', color='red')
    x1=np.linspace(0,5,10000)
    plt.fill_between(x1,0,gaus(x1,sigS,S,1/(np.sqrt(2*np.pi))), alpha=0.2, color='g')
    plt.legend(["Distribution", "Mean ="+str(round(S, 2)) , "S<5"], loc='upper right')
    area=trap(x1,gaus,sigS,S,1/(np.sqrt(2*np.pi)),10000)
    plt.title("Updated Significance PDF")
    return (1-(area)),sigS
prob=getP(*test1)
print(prob)




    

    


    

    
    
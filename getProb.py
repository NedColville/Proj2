# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 19:45:06 2022

@author: Ned

"""

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib
from scipy.optimize import curve_fit
matplotlib.rcParams['axes.formatter.useoffset'] = False
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
    return math.sqrt(1.56183e7*((np.exp(-ml/20))-np.exp(-mu/20)))

def main(ml,mu,nList=1000,nInteg=1000):
    m=np.linspace(ml,mu,nList)
    rootNB=rNb(ml,mu)
    nH=trap(m,gaus,*HP,nInteg)
    return nH/rootNB
def getVals(mu1,mu2,ml1,ml2,n):
    gridy=np.linspace(ml1,ml2,n)
    gridx=np.linspace(mu1,mu2,n)
    meshX,meshY=np.meshgrid(gridx,gridy)
    vals=np.zeros((len(gridx),len(gridy),2))
    for i in range(0,len(gridx)):
        for j in range(0,len(gridy)):
            vals[i][j]=[gridx[i],gridy[j]]
    bigval=np.zeros((len(gridx),len(gridy)))      
    for i in range(0,len(vals)):
        for j in range(0,len(vals[0])):
            bigval[i][j]=main(vals[i][j][0],vals[i][j][1])
    return gridx,gridy,bigval



def find_nearest(array, value):
    array = np.asarray(array)
    pos = (np.abs(array - value)).argmin()
    return pos

def biLinInterp(x1eval,x2eval,x1,x2,f):
    i=find_nearest(x1,x1eval)
    j=find_nearest(x2,x2eval)
    In11=((x1[i+1]-x1eval)*f[j][i]+((x1eval-x1[i])*f[j][i+1]))/(x1[i+1]-x1[i])
    In12=((x1[i+1]-x1eval)*f[j+1][i]+((x1eval-x1[i])*f[j+1][i+1]))/(x1[i+1]-x1[i])
    final=((x2[j+1]-x2eval)*In11+(x2eval-x2[j])*In12)/(x2[j+1]-x2[j])
    return final
   
def getmRange(crit,n):
    check1=check2=1
    ml1,ml2=122,125
    mu1,mu2=127,128
    counter=0
    while check1>crit and check2>crit:
        data=getVals(ml1,ml2,mu1,mu2,n)
        maxes=np.unravel_index(data[2].argmax(), data[2].shape)
        ml1=data[0][maxes[0]-1]
        ml2=data[0][maxes[0]+1]
        mu1=data[1][maxes[1]-1]
        mu2=data[1][maxes[1]+1]
        check1=mu2-mu1
        check2=ml2-ml1
        counter+=1
    print(counter)
        
    return [[ml1,ml2],[mu1,mu2]]
#test=getmRange(0.0001,30)
test=[[123.23748128395927, 123.23754914948938], [127.15847166826197, 127.15849429010534]]
"""
data1=getVals(*test[0],*test[1],100)
gridx=data1[0]
gridy=data1[1]
bigval=data1[2]
meshX,meshY=np.meshgrid(gridx,gridy)
plt.contour(gridx,gridy,bigval)
plt.xlabel("$m_l$", fontsize=20)
plt.ylabel("$m_u$", fontsize =20)
plt.xticks(rotation=90)
plt.title("Significance Contour")
plt.scatter(meshX,meshY,c=bigval,s=1 )
plt.colorbar()
"""

def D2CFD2(y21,y01,y12,y10,h):
    cfd1=(y21-y01)/h
    cfd2=(y12-y10)/h
    return [cfd1,cfd2]


def nextXY(ml,mu,mlList,muList,bigvals,h,alpha,crit):
    vals=[]
    check=0
    count=0
    while check==0:
        y21=biLinInterp(ml+h,mu,mlList,muList,bigvals)
        y01=biLinInterp(ml-h,mu,mlList,muList,bigvals)
        y12=biLinInterp(ml,mu+h,mlList,muList,bigvals)
        y10=biLinInterp(ml,mu-h,mlList,muList,bigvals)
        grad=D2CFD2(y21,y01,y12,y10,h)
        nextPos=np.array([ml,mu])+alpha*np.array(grad)
        SCur=biLinInterp(ml, mu, mlList,muList,bigvals)
        STemp=biLinInterp(nextPos[0],nextPos[1],mlList,muList,bigvals)    
        if abs(SCur-STemp)<crit and count>1000:
            check=1
        plt.quiver(ml,mu,grad[0],grad[1])
        ml,mu=nextPos[0],nextPos[1]
        vals.append([ml,mu])
        count+=1
    return ml,mu



#test1=nextXY(gridx[45],gridy[5],gridx,gridy,bigval,1e-8,0.01,1e-9)
test1=(123.2375145522137, 127.15848347691396)
plt.show()
print((rNb(*test1))**2)
def getP(ml,mu):
    Nb=(rNb(ml,mu))**2
    S=main(ml,mu)
    N=trap(np.linspace(ml,mu,1000),gaus,*HP,1000)+Nb
    sigS=(N/Nb+((N-Nb)**2)/(4*Nb**2))**0.5
    x=np.linspace(0,10,1000)
    y=gaus(x,sigS,S,1/(2*np.pi)**0.5)
    plt.plot(x,y)
    plt.xlabel("Significance")
    plt.ylabel("")
    plt.axvline(S,linestyle='dashed', color='red')
    x1=np.linspace(0,5,10000)
    plt.fill_between(x1,0,gaus(x1,sigS,S,1/(np.sqrt(2*np.pi))), alpha=0.2, color='g')
    plt.legend(["Distribution", "Mean ="+str(round(S, 2)) , "S<5"], loc='upper right')
    area=trap(x1,gaus,sigS,S,1/(np.sqrt(2*np.pi)),10000)
    plt.title("Significance PDF")
    return (1-(area)),sigS
prob=getP(*test1)
print(prob)



    
    



    




        
    

    

        
        













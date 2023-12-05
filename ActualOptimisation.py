# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 14:08:56 2022

@author: Ned
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:02:26 2022

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

def main(ml,mu,nList=1000,nInteg=100):
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
test=getmRange(0.00001,30)
print(test)
data1=getVals(*test[0],*test[1],50)
gridx=data1[0]
gridy=data1[1]
bigval=data1[2]
meshX,meshY=np.meshgrid(gridx,gridy)
plt.contour(gridx,gridy,bigval)
plt.xlabel(r"$m_l, GeV/c^{2}$", fontsize=14)
plt.xticks(fontsize=7)
plt.ylabel(r"$m_u, GeV/c^{2}$", fontsize=14)
plt.title("Significance Contour")
plt.scatter(meshX,meshY,c=bigval,s=1 )
plt.colorbar()
def D2CFD(x1,x2,y,h1,h2):        
    cfd1=(y(x1+h1,x2)-y(x1-h1,x2))/h1
    cfd2=(y(x1,x2+h2)-y(x1,x2-h2))/h2
    return [cfd1,cfd2]

def D2CFD2(y21,y01,y12,y10,h):
    cfd1=(y21-y01)/h
    cfd2=(y12-y10)/h
    return [cfd1,cfd2]


def nextXY(ml,mu,mlList,muList,bigvals,h,alpha,N):
    vals=[]
    for i in range(0,N):
        y21=main(ml+h,mu)
        y01=main(ml-h,mu)
        y12=main(ml,mu+h)
        y10=main(ml,mu-h)
        grad=D2CFD2(y21,y01,y12,y10,h)
        
        nextPos=np.array([ml,mu])+alpha*np.array(grad)
        plt.quiver(ml,mu,grad[0],grad[1],linewidth=0.05)
        ml,mu=nextPos[0],nextPos[1]
        vals.append([ml,mu])
    return ml,mu


print(gridx[30]-gridx[29])
test1=nextXY(gridx[49],gridy[0],gridx,gridy,bigval,0.00000005,0.05,200)
print(test1)
print((rNb(*test1))**2)
print(main(*test1))
#sigNh=np.sqrt((1.4+rNb(*test1))**2+(rNb(*test1))**2)
#x=np.linspace(*test1,1000)
#Nh=trap(x,gaus,*HP,1000)
maxes=np.unravel_index(bigval.argmax(), bigval.shape)
print(gridx[maxes[0]],gridy[maxes[1]])
print(gridx[maxes[0]-1],gridx[maxes[0]+1])
print(gridy[maxes[1]-1],gridy[maxes[1]+1])







        
    

    

        
        













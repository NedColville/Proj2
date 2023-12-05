# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 10:45:00 2022

@author: Ned
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib
from matplotlib import cm

HP=[1.4,125.1,470/(math.sqrt(2*np.pi))]
def rectangle(x,y,param,N):
    val=0
    a=x[0]
    b=x[len(x)-1]
    h=(b-a)/N
    x1=a+h/2
    for i in range(0,N):
        val+= h*y(x1,*param)
        x1+=h
    return val
def gaus(x,sig,mean,norm):
    norm=norm/sig
    return norm*np.exp(-(x-mean)**2/(2*sig**2))

def rNb(ml,mu):
    return math.sqrt((1500*20)*(np.exp((125.1-ml)/20)-np.exp((125.1-mu)/20)))

def main(ml,mu,nList=1000,nInteg=100):
    m=np.linspace(ml,mu,nList)
    rootNB=rNb(ml,mu)
    nH=rectangle(m,gaus,HP,nInteg)
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
            bigval[i][j]=main(vals[j][i][0],vals[j][i][1],1000,1000)
    return gridx,gridy,bigval



data1=getVals(115,125.05,125.15,136.15,50)

gridx=data1[0]
gridy=data1[1]
bigval=data1[2]
meshX,meshY=np.meshgrid(gridx,gridy)
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(meshX,meshY, bigval,linewidth=0, cmap='viridis',antialiased=False)
ax.set_xlabel(r"$m_l, GeV/c^{2}$")
ax.set_ylabel(r"$m_u, GeV/c^{2}$")
ax.set_zlabel("Significance")
plt.colorbar()
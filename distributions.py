# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 18:17:52 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
import math

def basic(m):
    return 1500*np.exp(-(m-125.1)/(20))
m=np.linspace(100,160,5000)

plt.plot(m,basic(m), linestyle='--', color='blue')
plt.ylim(0,3500)
plt.xlim(110,150)
def gaus(x,sig,mean,norm):
    norm=norm/sig
    return norm*np.exp(-(x-mean)**2/(2*sig**2))

def higgsDistrib(exp,gaus,m):
    return exp(m)+gaus(m,1.4,125.1,470*((2*np.pi)**-0.5))

def rNb(ml,mu):
    return math.sqrt((1500*20)*(np.exp((125.1-ml)/20)-np.exp((125.1-mu)/20)))
plt.plot(m,higgsDistrib(basic,gaus,m), color='blue')

def trap(x,f,sig,mean,norm,N):
    a=x[0]
    b=x[len(x)-1]
    h=(b-a)/(N)
    val=0
    for i in range(0,N):
        val+=0.5*h*(f(a+i*h,sig,mean,norm)+f(a+(i+1)*h,sig,mean,norm))
    return val
ml, mu=123, 127
plt.axvline(ml,color='r', linestyle='dashed', label='ml')
plt.axvline(mu,color='r', linestyle='dashed', label='mu')
plt.legend(["Background", "Higgs Curve + Background"])
plt.xlabel(r"Photon Pair Mass (GeV/$c^{2}$)")
plt.ylabel(r"Number per GeV/$c^{2}$")
rNb=rNb(ml,mu)
Nh=trap(m,gaus,1.4,125.1,470*(2*np.pi)**-0.5,100)
print(Nh/rNb)

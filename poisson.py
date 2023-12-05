# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:56:14 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
x=np.arange(0,21)

def poissons(x,mean):
    vals=[]
    for i in range(0,len(x)):
        vals.append(((mean**(x[i]))*np.exp(-1*mean))/(np.math.factorial(x[i])))
    return vals
plt.plot(x,poissons(x,5))
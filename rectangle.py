# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 20:13:35 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
x=np.linspace(0,5,1000)
def rectangle(x,y,N):
    val=0
    a=x[0]
    b=x[len(x)-1]
    h=(b-a)/N
    x1=a+h/2
    for i in range(0,N):
        val+= h*y(x1)
        x1+=h
    return val




    
    
    
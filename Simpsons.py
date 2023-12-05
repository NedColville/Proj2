# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 18:11:50 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt


def simp(x,f,N):
    a=x[0]
    b=x[len(x)-1]
    h=(b-a)/N
    val=0
    xs=0
    ms=0
    for i in range(1,N):
        xs= xs+f(a+i*h)  
    for i in range(0,N):
        ms= ms+f(a+(h/2)+i*h)
    val=(h/6)*(f(a)+4*ms+2*xs+f(b))
    return val




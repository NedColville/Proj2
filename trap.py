# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 11:16:19 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
def trap(x,f,N):
    a=x[0]
    b=x[len(x)-1]
    h=(b-a)/(N)
    val=0
    for i in range(0,N):
        val+=0.5*h*(f(a+i*h)+f(a+(i+1)*h))
    return val
   


   
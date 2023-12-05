# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 16:27:01 2022

@author: Ned
"""
import numpy as np
import matplotlib.pyplot as plt
import random

def MCInteg1d(x,f,N):
    tot=0
    V=x[len(x)-1]-x[0]
    for j in range(0,N):
        i=random.randint(0,len(x)-1)
        tot+=f(x[i])
    return tot*(V/N)               

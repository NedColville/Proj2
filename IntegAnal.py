# -*- coding: utf-8 -*-
"""
Created on Tue Dec 13 20:52:10 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
from rectangle import rectangle
from Simpsons import simp
from trap import trap
from monte import MCInteg1d as monte
import math
import statistics as stat
erf=math.erf
methods=[rectangle,trap,simp,monte]
methodsArg=["Rectangle","Trapezium","Simpsons","Monte Carlo"]
def gaus(x,sig=1,mean=0,norm=(2*np.pi)**-0.5):
    norm = norm/sig
    return norm*(np.exp((-(x-mean)**2)/(2*sig**2)))
def getx(a=5):
    return np.linspace(0,a,1000)
def getVals(x,a=5):
    values=[]
    for m in methods:
        values.append(abs((0.5*math.erf(a)-m(x,gaus,100)))/0.5*math.erf(a))
    return values
plt.plot(getx(),gaus(getx()),1000)
plt.ylim(0,1)
print(getVals(getx()))
values=[]
best=[]
bestarg=[]
aList=[]
for i in range(0,100):
    a=3+0.05*i
    aList.append(a)
    temp=getVals(getx(a),a)
    values.append(temp)
    best.append(np.min(temp))
    bestarg.append(np.argmin(temp))
freqs=[bestarg.count(0),bestarg.count(1),bestarg.count(2),bestarg.count(3)]
plt.show()
xpos=[i for i, _ in enumerate(methods)]
plt.bar(xpos,freqs)
plt.ylabel("Frequency of Highest Accuracy")
plt.xticks(xpos,methodsArg)
plt.show()
a=[]
val1=[]
val2=[]
val3=[]
val4=[]
aList=[]
big=[]
for i in range(0,8):
    a=5+0.5*i
    aList.append(a)
    temp=getVals(getx(a),a)
    big.append(temp)
    val1.append(temp[0])
    val2.append(temp[1])
    val3.append(temp[2])
    val4.append(temp[3])
plt.plot(aList,val1, 'x')
plt.plot(aList,val2,'x')
plt.plot(aList,val3,'x')
plt.legend(["Rectangle", "Trapezium", "Simpsons"])
plt.yscale('log')
plt.xlabel("Upper Integration Limit")
plt.ylabel("Fractional Error")

plt.show()
      




    


    
    




    

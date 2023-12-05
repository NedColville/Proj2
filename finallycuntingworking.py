# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 19:53:10 2022

@author: Ned
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')

X1=np.linspace(0,5,100)
X2=np.linspace(0,5,100)
x1,x2 = np.meshgrid(X1,X2)
def f(x1,x2):
    return x1**2+5*np.sin(x2)
Z=f(x1,x2)
plt.contour(x1,x2,Z, levels=20)



    



#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.plot_surface(x1,x2, f(x1,x2), cmap=cm.coolwarm,linewidth=0, antialiased=False)

#test=Data2CFD(X1,X2,*start,Z,1)
grad=np.gradient(Z)
#test=np.array(test)
#arrow=X1[start[0]],X2[start[1]],test[0],test[1]
#plt.quiver(*arrow)
plt.colorbar()



def Data2CFD(x1,x2,x1i,x2i,f,alpha):
    grad1=(f[x2i+1][x1i]-f[x2i-1][x1i])/(x2[x2i+1]-x2[x2i-1])
    grad2=(f[x2i][x1i+1]-f[x2i][x1i-1])/(x1[x1i+1]-x1[x1i-1])
    return [grad2,grad1]

start=[40,19]
test=Data2CFD(X1,X2,start[0],start[1],Z,1)
arrow=X1[start[0]],X2[start[1]],test[0],test[1]
plt.quiver(*arrow)





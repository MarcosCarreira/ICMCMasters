# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 19:29:07 2016

@author: marcoscscarreira
"""

import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
#import scipy as sc
#import seaborn as sns

def grid(N,T,Nj,S,dx,vol,r,q,K,output='V'):
    dt=T/N
    if dx<vol*np.sqrt(3*dt):
        print('Convergence error')
    mu=(r-q)-0.5*vol**2
    emdx=np.exp(-dx)
    ppu=0.5*dt*((vol/dx)**2+mu/dx)
    ppm=1.0-dt*(vol/dx)**2-r*dt
    ppd=0.5*dt*((vol/dx)**2-mu/dx)
    St=np.full(2*Nj+1,0.)
    St[0]=S*np.exp(Nj*dx)
    for j in np.arange(1,2*Nj+1):
        St[j]=St[j-1]*emdx
    Ct=np.full((2*Nj+1,N+1),np.nan)
    for j in np.arange(0,2*Nj+1):
        Ct[j,N]=np.max([St[j]-K,0])
    for i in np.arange(N-1,-1,-1):
        for k in np.arange(N-i,2*Nj+1-(N-i)):
            Ct[k,i]=ppu*Ct[k-1,i+1]+\
                ppm*Ct[k,i+1]+ppd*Ct[k+1,i+1]
        #Boundary
    #Boundary
    if output=='G':
        return Ct
    else:
        return Ct[Nj,0]

grid(3,1,3,100,0.2,0.20,0.06,0.03,100,'G')
grid(3,1,3,100,0.2,0.20,0.06,0.03,100)  

grid(6,1,6,100,0.2,0.20,0.06,0.03,100,'G')
grid(6,1,6,100,0.2,0.20,0.06,0.03,100)

grid(1000,1,1000,100,0.015,0.20,0.06,0.03,100)
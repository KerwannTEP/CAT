#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 21:11:30 2020

@author: kerwann
"""

import numpy as np
import matplotlib.pyplot as plt

def f(r,L):
    return 1/pow(1+r**2,1/2)-L**2/(2*r**2)

def fp(r,L):
    return -r/pow((1+r**2),3/2) + L**2/(r**3)

def fpp(r,L):
    return -(pow((1+r**2),3/2)-3*r**2*pow((1+r**2),1/2))/(1+r**2)**3 - 3*L**2/(r**4)

def newton(L,eps=pow(10,-6)):
    r = pow(L,2/3)
    while(fp(r+eps,L) > 0):
        r = r - fp(r,L)/fpp(r,L)
        #print(r)
        #print(fp(r,L))
    return r+eps/2

def Ec(L,eps=pow(10,-6)):
    return f(newton(L,eps),L)

def CircOrbits(Lmax,N=100):
    L = np.linspace(0,Lmax,N)
    E = [1]
    for i in range(1,N):
        E.append(Ec(L[i]))
        
    plt.figure()
    plt.plot(E,L)
    plt.xlabel("Binding energy per unit mass")
    plt.ylabel("Angular momentum per unit mass")
    plt.title("Circular orbits in (E,L) parameter space")
    plt.xlim(0,1)
    plt.ylim(0,Lmax)
    plt.legend()
    plt.show()
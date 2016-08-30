#!/usr/bin/python
import numpy as np
import math as m
import argparse
import matplotlib
from pylab import *
def rescale_Y(data,nu,beta,L):

    new_data = np.zeros(len(data))
    new_err  = np.zeros(len(data))
    
    for i in range(len(data)):
        
        new_data[i] = data[i,1]*L**(beta/nu)
        new_err[i]  = data[i,2]*L**(beta/nu)

    return [new_data,new_err]

def rescale_X(data,nu,hC,L):

    new_X = np.zeros(len(data))

    for i in range(len(data)):

        new_X[i] = ((data[i,0]-hC)/hC)*L**(1.0/nu)

    return new_X



if __name__ == "__main__":

    
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', help='Linear size of the model',type=str)
    args = parser.parse_args()

    if (args.m == "ising"):
        sizes = [500,1000,2000,4000]
        baseName = 'data/Ising/Collapse_Ising_L'
        hC = 1.0
        nu = 1.0
        beta = 0.125

    if (args.m == "random"):
        sizes = [32,64,128,256]
        baseName = 'data/randomIsing/collapse_L'
        hC = 1.0
        nu = 2.0
        beta = (3-m.sqrt(5))/2

    colors = ["ko", "ro", "bo", "go","mo","yo","co","ko"]
    fig = figure(1,figsize=(8,6))

    for i,L in enumerate(sizes):
                
        dataFileName = baseName
        dataFileName += str(L)
        dataFileName += '.dat'

        dataFile = open(dataFileName,'r')
        data = np.loadtxt(dataFileName)
        
        [rDat,rErr] = rescale_Y(data,nu,beta,L)
        x = rescale_X(data,nu,hC,L)

        plot(x[:],rDat,colors[i],marker='^')
        
    show()
        
#p_c=1.0
#nu=2.0
#beta=(3-sqrt(5))/2
#
#
#for i in range(0,10):
#    x16.append(((data16[i,0]-p_c)/p_c)*16**(1/nu)) 
#    x32.append(((data32[i,0]-p_c)/p_c)*32**(1/nu)) 
#    x64.append(((data64[i,0]-p_c)/p_c)*64**(1/nu)) 
#    x128.append(((data128[i,0]-p_c)/p_c)*128**(1/nu))
#    y16.append(data16[i,2]*16**(beta/nu)) 
#    y32.append(data32[i,2]*32**(beta/nu)) 
#    y64.append(data64[i,2]*64**(beta/nu)) 
#    y128.append(data128[i,2]*128**(beta/nu))
     
# 
#
#plt.yscale('linear', nonposy='clip')
#
#plt.xlabel('$L^{1/2} (h-h_c)/h_c$', fontsize=16)
#plt.ylabel('$\Delta\lambda L^{(3-\sqrt{5})/4}$', fontsize=16)
#
#plt.plot(x16[:], y16[:],'yo', label='L=16')
#plt.plot(x32[:], y32[:],'go', label='L=32')
#plt.plot(x64[:], y64[:],'bo', label='L=64')
#plt.plot(x128[:], y128[:],'mo', label='L=128')
#
#
#plt.legend(loc='lower right', numpoints =1)
#
#
#
##plt.plot(x16[:], y16[:],'ro',x32[:], y32[:],'bo',x64[:], y64[:],'go',x128[:], y128[:],'yo')
#
#plt.savefig("FSS")
#plt.clf()

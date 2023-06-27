#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 11:56:30 2022

@author: accdavlo
"""

# We need a couple of packages in this chapter
import numpy as np  
# This is the basic package in python with all the numerical functions

import matplotlib.pyplot as plt 
# This package allows to  plot

from nodepy import rk
#This package already implemented some functions for Runge Kutta and multistep methods
from src import DeC
from src import RungeKutta
from src import ODEproblems
from src.DeC import *
from src.ODEproblems import ODEproblem
import time


from timeit import timeit
import csv
import os

folder = "figures"

try: 
    os.mkdir(folder) 
except OSError as error: 
    print(error)  

#%% Nodepy stuff

colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

max_order=14
    
methods = [DeC, DeC_staggered, DeC_staggered_f, DeC_small_sub, DeC_small_sub_staggered, DeC_small_sub_staggered_f]

distributions = ["equispaced","gaussLobatto"]

def subtimesteps(dist,order):
    if dist == "equispaced":
        return order-1
    elif dist=="gaussLobatto":
        return int((order+1)//2)


for dist in distributions:

    for im, method in enumerate(methods):
        plt.figure(im,figsize=(8,10))
        
    orders = np.zeros((len(methods),max_order),dtype=np.int32)
    stages = np.zeros((len(methods),max_order))
    legends=[None,None,None,None]
    for order in range(3,max_order):
        print("!-----------------------------------!")
        print("!-----Order %d-----------------------!"%order)
        print("!-----------------------------------!")
        
        for im, method in enumerate(methods):
            M_sub=subtimesteps(dist,order)
            dec_method = method(M_sub, order, dist)
            A,b,c=dec_method.compute_RK_from_DeC()
            rkStaggeredDeC = rk.ExplicitRungeKuttaMethod(A,b)
            rkStaggeredDeC.name=dec_method.name+str(order)
            exp_order = rkStaggeredDeC.order()
            print(rkStaggeredDeC.name+" has order "+str(exp_order))
            print(rkStaggeredDeC.name+" has %d stages"%(len(rkStaggeredDeC.c)))
            orders[im, order] = exp_order
            stages[im, order] = len(rkStaggeredDeC.c)

            p1=rkStaggeredDeC.plot_stability_region(bounds=[-10,5,-10,10],filled=False,fignum=im,color=colors[order], linewidth=2,linestyle ='solid')

            if order<=4:
                p,d = rkStaggeredDeC.stability_function()
                print(p)
                
    for im, method in enumerate(methods):
        dec_method = method(2,3,dist)
        plt.figure(im)
        plt.title(None)
        plt.tight_layout()
        plt.grid(True)
        plt.savefig(folder+"/stability_region_"+dec_method.name+".pdf")
            
    plt.show()

    plt.close('all')


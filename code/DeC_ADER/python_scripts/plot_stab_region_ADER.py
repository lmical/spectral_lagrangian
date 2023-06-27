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


import time


from timeit import timeit
import csv
import os
from distutils.dir_util import mkpath

#This package already implemented some functions for Runge Kutta and multistep methods
from src import DeC
from src import RungeKutta
from src import ODEproblems
from src.DeC import *

from src.ADER import ADER, ADER_L2, ADER_u

table_folder = "tables/"
mkpath(table_folder)
figure_folder = "figures/"
mkpath(figure_folder)




#%% Nodepy stuff

colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

max_order=14
    
methods = [ADER, ADER_u, ADER_L2]

distributions = ["equispaced","gaussLobatto","gaussLegendre"]


for dist in distributions:

    for im, method in enumerate(methods):
        plt.figure(im,figsize=(8,10))
        
    orders = np.zeros((len(methods),max_order),dtype=np.int32)
    stages = np.zeros((len(methods),max_order))
    legends=[None,None,None,None]
    for order in range(2,max_order):
        print("!-----------------------------------!")
        print("!-----Order %d-----------------------!"%order)
        print("!-----------------------------------!")
        
        for im, method in enumerate(methods):
            ader_method = method(-1, order, dist)
            A,b,c=ader_method.compute_RK()
            ADERRK = rk.ExplicitRungeKuttaMethod(A,b.flatten())
            ADERRK.name=ader_method.name+str(order)
            exp_order = ADERRK.order()
            print(ADERRK.name+" has order "+str(exp_order))
            print(ADERRK.name+" has %d stages"%(len(ADERRK.c)))
            orders[im, order] = exp_order
            stages[im, order] = len(ADERRK.c)

            p1=ADERRK.plot_stability_region(bounds=[-10,5,-10,10],filled=False,fignum=im,color=colors[order], linewidth=2,linestyle ='solid')

            if order<=4:
                p,d = ADERRK.stability_function()
                print(p)
                
    for im, method in enumerate(methods):
        ader_method = method(2,3,dist)
        plt.figure(im)
        plt.title(None)
        plt.tight_layout()
        plt.grid(True)
        plt.savefig(figure_folder+"stability_region_"+ader_method.name+".pdf")
            
    plt.show()

    plt.close('all')


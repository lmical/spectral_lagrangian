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


#%%
pr=ODEproblem("vibratingDamped")

tt=np.linspace(0,pr.T_fin,10)   #Plot the evolution for order 8
dec2 = DeC_small_sub(2, 3, "equispaced")  #_staggered
tt,uu=dec2.dec(pr.flux, tt, pr.u0)
#%% Nodepy stuff

colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

max_order=14
    
methods = [DeC, DeC_staggered, DeC_staggered_f, DeC_small_sub, DeC_small_sub_staggered, DeC_small_sub_staggered_f]

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
        dec_method = method(order-1, order, "equispaced")
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
    dec_method = method(2,3,"equispaced")
    plt.figure(im)
    plt.title(None)
    plt.tight_layout()
    plt.savefig("stability_region_"+dec_method.name+".pdf")
        
plt.show()

stop()
#%%
dec_small=DeC_small_sub_staggered(2,3,"equispaced")
dec_small.thetas

pr=ODEproblem("linear_system2")

order = 3
tt=np.linspace(0,pr.T_fin,10)   #Plot the evolution for order 8
dec8 = DeC_small_sub_staggered(order-1, order, "equispaced")
order=4
dec_method = DeC_small_sub_staggered(order-1, order, "equispaced")
A,b,c=dec_method.compute_RK_from_DeC()
rkStaggeredDeC = rk.ExplicitRungeKuttaMethod(A,b)
rkStaggeredDeC.name="Stag-DeC-equi"+str(order)
print(rkStaggeredDeC.name+" has order "+str(rkStaggeredDeC.order()))
print(rkStaggeredDeC.name+" has %d stages"%(len(rkStaggeredDeC.c)))
print(A)
print(b)


#%%
dec=DeC_small_sub_staggered_f(4,5,"equispaced")
dec.compute_RK_from_DeC()

#%% Staggered subtimesteps
#Test convergence of DeC for several orders
pr=ODEproblem("nonlinear_scalar")

tt=np.linspace(0,pr.T_fin,10)   #Plot the evolution for order 8
dec8 = DeC_staggered(7, 8, "equispaced")
tt,uu=dec8.dec(pr.flux, tt, pr.u0)
# plt.plot(tt,uu[0,:])
# plt.plot(tt,uu[1,:])
# plt.show()

colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

def compute_integral_error(c,c_exact):  # c is dim x times
    times=np.shape(c)[1]
    error=0.
    for t in range(times):
        error = error + np.linalg.norm(c[:,t]-c_exact[:,t],2)**2.
    error = np.sqrt(error/times) 
    return error

NN=5
dts=[pr.T_fin/2.0**k for k in range(2,2+NN)]
errorsDeCStag=np.zeros(len(dts))
errorsDeC=np.zeros(len(dts))
timesDeC=np.zeros(len(dts))
timesDeCStag=np.zeros(len(dts))

fig1,ax1 = plt.subplots(1,1)
fig2,ax2 = plt.subplots(1,1)


# Compute and plot the errors 
for order in range(2,10):
    myFile = open(f'convergence_{pr.name}_order_{order}.csv', 'w')
    with myFile:
        writer = csv.writer(myFile)
        myData = ["dt", "error DeC", "error Stag DeC", "time DeC", "time Stag DeC"]
        writer.writerow(myData)
        print("order=",order)
        for k in range(NN):
            dt0=dts[k]
            print("dt is ",dt0)
            tt=np.arange(0,pr.T_fin,dt0)
            M_sub = order-1
            K_iter = order

            # Staggered DeC
            dec = DeC_small_sub_staggered_f(M_sub, K_iter, "equispaced")
            tic=time.time()
            dec.dec(pr.flux, tt, pr.u0)
            toc=time.time()
            res=toc-tic
            timesDeCStag[k] = res
            t2,U2=dec.dec(pr.flux, tt, pr.u0)
            u_exact=pr.exact_solution_times(pr.u0,tt)
            errorsDeCStag[k]=compute_integral_error(U2,u_exact)

            #Normal DeC
            dec = DeC_small_sub(M_sub, K_iter, "equispaced")
            tic=time.time()
            dec.dec(pr.flux, tt, pr.u0)
            toc=time.time()
            res=toc-tic
            timesDeC[k] = res
            t2,U2=dec.dec(pr.flux, tt, pr.u0)
            u_exact=pr.exact_solution_times(pr.u0,tt)
            errorsDeC[k]=compute_integral_error(U2,u_exact)
            
            myData = [dt0, errorsDeC[k], errorsDeCStag[k], timesDeC[k], timesDeCStag[k]]
            writer.writerow(myData)
            
            print("speed up factor ", timesDeC[k]/timesDeCStag[k] )
            
        ax1.loglog(dts,errorsDeC,"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errorsDeCStag,"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax2.loglog(timesDeC,errorsDeC,"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(timesDeCStag,errorsDeCStag,"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))


ax1.set_title("DeC error convergence")
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax1.set_ylabel("Error")
ax1.set_xlabel(r"$\Delta t$")
fig1.savefig(f"convergence_DeC_staggered_sub_{pr.name}.pdf")

ax2.set_title("DeC error convergence")
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
ax2.set_ylabel("Error")
ax2.set_xlabel("Computational time")
fig2.savefig(f"convergence_vs_time_DeC_staggered_sub_{pr.name}.pdf")

plt.show()




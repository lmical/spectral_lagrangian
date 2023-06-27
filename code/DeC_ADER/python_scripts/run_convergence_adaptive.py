# We need a couple of packages in this chapter
import numpy as np  
# This is the basic package in python with all the numerical functions

import matplotlib.pyplot as plt 
# This package allows to  plot

from timeit import timeit
import csv
from src import DeC
from src import RungeKutta
from src import ODEproblems
from src.DeC import *
from src.ODEproblems import ODEproblem

import os

folder = "figures/adaptive"

try: 
    os.mkdir(folder) 
except OSError as error: 
    print(error)  

methods = [DeC, DeC_staggered, DeC_staggered_f, DeC_small_sub, DeC_small_sub_staggered, DeC_small_sub_staggered_f]
methods_staggered = [DeC_staggered, DeC_staggered_f, DeC_small_sub_staggered, DeC_small_sub_staggered_f]
N_meth = len(methods)
N_meth_staggered = len(methods_staggered)

#%% 
#Test convergence of DeC for several orders
pr=ODEproblem("linear_system2")

distributions=["equispaced","gaussLobatto"]

def subtimesteps(dist,order):
    if dist == "equispaced":
        return order-1
    elif dist=="gaussLobatto":
        return int((order+1)//2)


colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

def compute_integral_error(c,c_exact):  # c is dim x times
    times=np.shape(c)[1]
    error=0.
    for t in range(times):
        error = error + np.linalg.norm(c[:,t]-c_exact[:,t],2)**2.
    error = np.sqrt(error/times) 
    return error

adaptive_tolerance=1e-8
adaptive_max_order = 16

NN=5
dts=[pr.T_fin/2.0**k for k in range(2,2+NN)]


for dist in distributions:
    print(dist)
    errors=np.zeros((N_meth,len(dts)))
    times =np.zeros((N_meth,len(dts)))

    fig1,ax1 = plt.subplots(1,1)
    fig2,ax2 = plt.subplots(1,1)
    fig3,ax3 = plt.subplots(1,1)
    fig4,ax4 = plt.subplots(1,1)



    # Compute and plot the errors 
    for order in range(2,10):
        myFile = open(folder+"/convergence_DeC_"+dist+f"_{pr.name}_order_{order}.csv", 'w')
        with myFile:
            writer = csv.writer(myFile)
            myData = ["dt"]
            for met in methods:
                myData.append("error "+met(1,2,dist).name)            
            for met in methods:
                myData.append("time "+met(1,2,dist).name)
            writer.writerow(myData)
            print("order=",order)
            for k in range(NN):
                dt0=dts[k]
                print("dt is ",dt0)
                tt=np.arange(0,pr.T_fin,dt0)
                M_sub = subtimesteps(dist,order)
                K_iter = order

                for im, method in enumerate(methods):
                    # Staggered DeC
                    dec = method(M_sub, K_iter, dist)
                    res = %timeit -o dec.dec(pr.flux, tt, pr.u0)
                    times[im,k] = res.average
                    t2,U2=dec.dec(pr.flux, tt, pr.u0)
                    u_exact=pr.exact_solution_times(pr.u0,tt)
                    errors[im,k]=compute_integral_error(U2,u_exact)
                    print(f"error {errors[im,k]}")
                
                
                myData = [dt0]
                for im in range(N_meth):
                    myData.append(errors[im,k])
                for im in range(N_meth):
                    myData.append(times[im,k])
                writer.writerow(myData)
                
                
                for im, method in enumerate(methods):
                    print("speed up factor "+method(1,2,dist).name, times[(im//3)*3,k]/times[im,k] )
                
            ax1.loglog(dts,errors[0,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

            ax2.loglog(times[0,:],errors[0,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax2.loglog(times[1,:],errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax2.loglog(times[2,:],errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                    
            ax3.loglog(dts,errors[3,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

            ax4.loglog(times[3,:],errors[3,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(times[4,:],errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(times[5,:],errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

    myFile = open(folder+"/convergence_DeC_adaptive_"+dist+f"_{pr.name}.csv", 'w')
    with myFile:
        writer = csv.writer(myFile)
        myData = ["dt"]
        for met in methods_staggered:
            myData.append("error "+met(1,2,dist).name)       
            myData.append("time "+met(1,2,dist).name)       
            myData.append("p average"+met(1,2,dist).name)
            myData.append("p std"+met(1,2,dist).name)
        writer.writerow(myData)

        order = adaptive_max_order
        
        times = np.zeros((len(methods_staggered),NN))
        errors = np.zeros((len(methods_staggered),NN))
        p_ave = np.zeros((len(methods_staggered),NN))
        p_std = np.zeros((len(methods_staggered),NN))

        print("order=",order)
        for k in range(NN):
            dt0=dts[k]
            print("dt is ",dt0)
            tt=np.arange(0,pr.T_fin,dt0)
            M_sub = subtimesteps(dist,order)
            K_iter = order

            for im, method in enumerate(methods_staggered):
                # Staggered DeC
                dec = method(M_sub, K_iter, dist)
                res = %timeit -o dec.dec_order_control(pr.flux, tt, pr.u0, adaptive_tolerance)
                times[im,k] = res.average
                t2,U2, p_dist =dec.dec_order_control(pr.flux, tt, pr.u0, adaptive_tolerance)

                #Plot of iterations vs time
                plt.figure(320)
                plt.plot(t2[1:],p_dist[1:])
                plt.xlabel("t")
                plt.ylabel("Iterations")
                plt.savefig(folder+"/p_vs_time_adaptive_"+dist+f"_{pr.name}_{dec.name}_k_{k}.pdf")
                plt.close(320)

                p_ave[im,k] = np.mean(p_dist[1:])
                p_std[im,k] = np.std(p_dist[1:])
                u_exact=pr.exact_solution_times(pr.u0,tt)
                errors[im,k]=compute_integral_error(U2,u_exact)
                print(f"error {errors[im,k]}")
            
            
            myData = [dt0]
            for im in range(N_meth_staggered):
                myData.append(errors[im,k])
                myData.append(times[im,k])
                myData.append(p_ave[im,k])
                myData.append(p_std[im,k])
            writer.writerow(myData)
            
            
        ax1.loglog(dts,errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))

        ax2.loglog(times[0,:],errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(times[1,:],errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                
        ax3.loglog(dts,errors[2,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,errors[3,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
#        ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax4.loglog(times[2,:],errors[2,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax4.loglog(times[3,:],errors[3,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

        

    ax1.set_title("DeC error convergence")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Error")
    ax1.set_xlabel(r"$\Delta t$")
    fig1.set_tight_layout(True)
    fig1.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

    ax2.set_title("DeC error convergence")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel("Error")
    ax2.set_xlabel("Computational time")
    fig2.set_tight_layout(True)
    fig2.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")


    ax3.set_title("DeC error convergence")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylabel("Error")
    ax3.set_xlabel(r"$\Delta t$")
    fig3.set_tight_layout(True)
    fig3.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

    ax4.set_title("DeC error convergence")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_ylabel("Error")
    ax4.set_xlabel("Computational time")
    fig4.set_tight_layout(True)
    fig4.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")


    plt.show(block=False)

    plt.close('all')

#%% 
#Test convergence of DeC for several orders
pr=ODEproblem("vibratingDamped")


colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

def compute_integral_error(c,c_exact):  # c is dim x times
    times=np.shape(c)[1]
    error=0.
    for t in range(times):
        error = error + (c[0,t]-c_exact[0,t])**2.
    error = np.sqrt(error/times) 
    return error

NN=5
dts=[pr.T_fin/2.0**k for k in range(2,2+NN)]

for dist in distributions:

    errors=np.zeros((N_meth,len(dts)))
    times =np.zeros((N_meth,len(dts)))

    fig1,ax1 = plt.subplots(1,1)
    fig2,ax2 = plt.subplots(1,1)
    fig3,ax3 = plt.subplots(1,1)
    fig4,ax4 = plt.subplots(1,1)



    # Compute and plot the errors 
    for order in range(2,10):
        myFile = open(folder+"/convergence_DeC_"+dist+f"_{pr.name}_order_{order}.csv", 'w')
        with myFile:
            writer = csv.writer(myFile)
            myData = ["dt"]
            for met in methods:
                myData.append("error "+met(1,2,dist).name)            
            for met in methods:
                myData.append("time "+met(1,2,dist).name)
            writer.writerow(myData)
            print("order=",order)
            for k in range(NN):
                dt0=dts[k]
                print("dt is ",dt0)
                tt=np.arange(0,pr.T_fin,dt0)
                M_sub = subtimesteps(dist,order)
                K_iter = order

                for im, method in enumerate(methods):
                    # Staggered DeC
                    dec = method(M_sub, K_iter, dist)
                    res = %timeit -o dec.dec(pr.flux, tt, pr.u0)
                    times[im,k] = res.average
                    t2,U2=dec.dec(pr.flux, tt, pr.u0)
                    u_exact=pr.exact_solution_times(pr.u0,tt)
                    errors[im,k]=compute_integral_error(U2,u_exact)
                    print(f"error {errors[im,k]}")
                
                
                myData = [dt0]
                for im in range(N_meth):
                    myData.append(errors[im,k])
                for im in range(N_meth):
                    myData.append(times[im,k])
                writer.writerow(myData)
                
                
                for im, method in enumerate(methods):
                    print("speed up factor "+method(1,2,dist).name, times[(im//3)*3,k]/times[im,k] )
                
            ax1.loglog(dts,errors[0,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

            ax2.loglog(times[0,:],errors[0,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax2.loglog(times[1,:],errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax2.loglog(times[2,:],errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                    
            ax3.loglog(dts,errors[3,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

            ax4.loglog(times[3,:],errors[3,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(times[4,:],errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(times[5,:],errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

    myFile = open(folder+"/convergence_DeC_adaptive_"+dist+f"_{pr.name}.csv", 'w')
    with myFile:
        writer = csv.writer(myFile)
        myData = ["dt"]
        for met in methods_staggered:
            myData.append("error "+met(1,2,dist).name)       
            myData.append("time "+met(1,2,dist).name)       
            myData.append("p average"+met(1,2,dist).name)
            myData.append("p std"+met(1,2,dist).name)
        writer.writerow(myData)

        order = adaptive_max_order
        
        times = np.zeros((len(methods_staggered),NN))
        errors = np.zeros((len(methods_staggered),NN))
        p_ave = np.zeros((len(methods_staggered),NN))
        p_std = np.zeros((len(methods_staggered),NN))

        print("order=",order)
        for k in range(NN):
            dt0=dts[k]
            print("dt is ",dt0)
            tt=np.arange(0,pr.T_fin,dt0)
            M_sub = subtimesteps(dist,order)
            K_iter = order

            for im, method in enumerate(methods_staggered):
                # Staggered DeC
                dec = method(M_sub, K_iter, dist)
                res = %timeit -o dec.dec_order_control(pr.flux, tt, pr.u0, adaptive_tolerance)
                times[im,k] = res.average
                t2,U2, p_dist =dec.dec_order_control(pr.flux, tt, pr.u0, adaptive_tolerance)

                #Plot of iterations vs time
                plt.figure(320)
                plt.plot(t2[1:],p_dist[1:])
                plt.xlabel("t")
                plt.ylabel("Iterations")
                plt.savefig(folder+"/p_vs_time_adaptive_"+dist+f"_{pr.name}_{dec.name}_k_{k}.pdf")
                plt.close(320)

                p_ave[im,k] = np.mean(p_dist[1:])
                p_std[im,k] = np.std(p_dist[1:])
                u_exact=pr.exact_solution_times(pr.u0,tt)
                errors[im,k]=compute_integral_error(U2,u_exact)
                print(f"error {errors[im,k]}")
            
            
            myData = [dt0]
            for im in range(N_meth_staggered):
                myData.append(errors[im,k])
                myData.append(times[im,k])
                myData.append(p_ave[im,k])
                myData.append(p_std[im,k])
            writer.writerow(myData)
            
            
        ax1.loglog(dts,errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))

        ax2.loglog(times[0,:],errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(times[1,:],errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                
        ax3.loglog(dts,errors[2,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,errors[3,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
#        ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax4.loglog(times[2,:],errors[2,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax4.loglog(times[3,:],errors[3,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                

    ax1.set_title("DeC error convergence")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Error")
    ax1.set_xlabel(r"$\Delta t$")
    fig1.set_tight_layout(True)
    fig1.savefig(folder+f"/convergence_DeC_"+dist+f"_staggered_{pr.name}.pdf")

    ax2.set_title("DeC error convergence")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel("Error")
    ax2.set_xlabel("Computational time")
    fig2.set_tight_layout(True)
    fig2.savefig(folder+f"/convergence_vs_time_DeC_"+dist+f"_staggered_{pr.name}.pdf")


    ax3.set_title("DeC error convergence")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylabel("Error")
    ax3.set_xlabel(r"$\Delta t$")
    fig3.set_tight_layout(True)
    fig3.savefig(folder+f"/convergence_DeC_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

    ax4.set_title("DeC error convergence")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_ylabel("Error")
    ax4.set_xlabel("Computational time")
    fig4.set_tight_layout(True)
    fig4.savefig(folder+f"/convergence_vs_time_DeC_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

    plt.show(block=False)
    plt.close('all')


#Test convergence of DeC for several orders
pr=ODEproblem("threeBodies")

tt=np.linspace(0,pr.T_fin,2049)   #Plot the evolution for order 8
dec8 = DeC(8, 9, "gaussLobatto")
tt_exact,uu_exact=dec8.dec(pr.flux, tt, pr.u0)
plt.plot(tt_exact,uu_exact[0,:])
plt.plot(tt_exact,uu_exact[1,:])
plt.show(block=False)



colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


def compute_integral_error(tt,c,tt_exact,c_exact):  # c is dim x times
    error=0.
    iex=0
    for it, t in enumerate(tt):
        while (t> tt_exact[iex]+1e-5):
            iex+=1
        error = error + np.linalg.norm(c[:,it]-c_exact[:,iex],2)**2.
    error = np.sqrt(error/len(tt)) 
    return error

NN=5
dts=[pr.T_fin/2.0**k for k in range(4,4+NN)]

for dist in distributions:

    errors=np.zeros((N_meth,len(dts)))
    times =np.zeros((N_meth,len(dts)))

    fig1,ax1 = plt.subplots(1,1)
    fig2,ax2 = plt.subplots(1,1)
    fig3,ax3 = plt.subplots(1,1)
    fig4,ax4 = plt.subplots(1,1)



    # Compute and plot the errors 
    for order in range(2,10):
        myFile = open(folder+"/convergence_DeC_"+dist+f"_{pr.name}_order_{order}.csv", 'w')
        with myFile:
            writer = csv.writer(myFile)
            myData = ["dt"]
            for met in methods:
                myData.append("error "+met(1,2,dist).name)            
            for met in methods:
                myData.append("time "+met(1,2,dist).name)
            writer.writerow(myData)
            print("order=",order)
            for k in range(NN):
                dt0=dts[k]
                print("dt is ",dt0)
                tt=np.arange(0,pr.T_fin,dt0)
                M_sub = subtimesteps(dist,order)
                K_iter = order

                for im, method in enumerate(methods):
                    # Staggered DeC
                    dec = method(M_sub, K_iter, dist)
                    res = %timeit -o dec.dec(pr.flux, tt, pr.u0)
                    times[im,k] = res.average
                    t2,U2=dec.dec(pr.flux, tt, pr.u0)
                    errors[im,k]=compute_integral_error(t2,U2,tt_exact,uu_exact)
                    print(f"error {errors[im,k]}")
                
                
                myData = [dt0]
                for im in range(N_meth):
                    myData.append(errors[im,k])
                for im in range(N_meth):
                    myData.append(times[im,k])
                writer.writerow(myData)
                
                
                for im, method in enumerate(methods):
                    print("speed up factor "+method(1,2,dist).name, times[(im//3)*3,k]/times[im,k] )
                
            ax1.loglog(dts,errors[0,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

            ax2.loglog(times[0,:],errors[0,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax2.loglog(times[1,:],errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax2.loglog(times[2,:],errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                    
            ax3.loglog(dts,errors[3,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

            ax4.loglog(times[3,:],errors[3,:],"-",color=colors[order],label="DeC(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(times[4,:],errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(times[5,:],errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

    myFile = open(folder+"/convergence_DeC_adaptive_"+dist+f"_{pr.name}.csv", 'w')
    with myFile:
        writer = csv.writer(myFile)
        myData = ["dt"]
        for met in methods_staggered:
            myData.append("error "+met(1,2,dist).name)       
            myData.append("time "+met(1,2,dist).name)       
            myData.append("p average"+met(1,2,dist).name)
            myData.append("p std"+met(1,2,dist).name)
        writer.writerow(myData)

        order = adaptive_max_order
        
        times = np.zeros((len(methods_staggered),NN))
        errors = np.zeros((len(methods_staggered),NN))
        p_ave = np.zeros((len(methods_staggered),NN))
        p_std = np.zeros((len(methods_staggered),NN))

        print("order=",order)
        for k in range(NN):
            dt0=dts[k]
            print("dt is ",dt0)
            tt=np.arange(0,pr.T_fin,dt0)
            M_sub = subtimesteps(dist,order)
            K_iter = order

            for im, method in enumerate(methods_staggered):
                # Staggered DeC
                dec = method(M_sub, K_iter, dist)
                res = %timeit -o dec.dec_order_control(pr.flux, tt, pr.u0, adaptive_tolerance)
                times[im,k] = res.average
                t2,U2, p_dist =dec.dec_order_control(pr.flux, tt, pr.u0, adaptive_tolerance)

                #Plot of iterations vs time
                plt.figure(320)
                plt.plot(t2[1:],p_dist[1:])
                plt.xlabel("t")
                plt.ylabel("Iterations")
                plt.savefig(folder+"/p_vs_time_adaptive_"+dist+f"_{pr.name}_{dec.name}_k_{k}.pdf")
                plt.close(320)

                p_ave[im,k] = np.mean(p_dist[1:])
                p_std[im,k] = np.std(p_dist[1:])
                errors[im,k]=compute_integral_error(t2,U2,tt_exact,uu_exact)
                print(f"error {errors[im,k]}")
            
            
            myData = [dt0]
            for im in range(N_meth_staggered):
                myData.append(errors[im,k])
                myData.append(times[im,k])
                myData.append(p_ave[im,k])
                myData.append(p_std[im,k])
            writer.writerow(myData)
            
            
        ax1.loglog(dts,errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))

        ax2.loglog(times[0,:],errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(times[1,:],errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                
        ax3.loglog(dts,errors[2,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,errors[3,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
#        ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax4.loglog(times[2,:],errors[2,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax4.loglog(times[3,:],errors[3,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                

    ax1.set_title("DeC error convergence")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Error")
    ax1.set_xlabel(r"$\Delta t$")
    fig1.set_tight_layout(True)
    fig1.savefig(folder+f"/convergence_DeC_"+dist+f"_staggered_{pr.name}.pdf")

    ax2.set_title("DeC error convergence")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel("Error")
    ax2.set_xlabel("Computational time")
    fig2.set_tight_layout(True)
    fig2.savefig(folder+f"/convergence_vs_time_DeC_"+dist+f"_staggered_{pr.name}.pdf")


    ax3.set_title("DeC error convergence")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylabel("Error")
    ax3.set_xlabel(r"$\Delta t$")
    fig3.set_tight_layout(True)
    fig3.savefig(folder+f"/convergence_DeC_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

    ax4.set_title("DeC error convergence")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_ylabel("Error")
    ax4.set_xlabel("Computational time")
    fig4.set_tight_layout(True)
    fig4.savefig(folder+f"/convergence_vs_time_DeC_"+dist+f"_staggered_small_sub_{pr.name}.pdf")


    plt.show(block=False)


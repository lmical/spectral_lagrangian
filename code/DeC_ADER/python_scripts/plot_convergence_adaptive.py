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

print("Here")

#adaptive_tolerance=1e-8
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

    print("There")


    # Compute and plot the errors 
    for order in range(3,10):
        fileName = folder+"/convergence_DeC_"+dist+f"_{pr.name}_order_{order}.csv"
        print("Evverywhere")

        errors = np.zeros((len(methods),len(dts)))
        times  = np.zeros((len(methods),len(dts)))

        with open(fileName) as myFile:
            reader = csv.reader(myFile)
            for ir, row in enumerate(reader):
                if ir!=0:
                    dt0 = float(row[0])
                    if abs(dt0 - dts[ir-1])>1e-14:
                        raise ValueError("Wrong dt!!!!") 
                    errors_all_methods = [float(z) for z in row[1:len(methods)+1]]
                    times_all_methods = [float(z) for z in row[len(methods)+1:2*len(methods)+1]]
                    for im, method in enumerate(methods):
                        times[im,ir-1] = times_all_methods[im]
                        errors[im,ir-1]= errors_all_methods[im] 
                        print(f"error {errors[im,ir-1]}")
                
               
        ax1.loglog(dts,errors[0,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax1.loglog(dts,errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax2.loglog(times[0,:],errors[0,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax2.loglog(times[1,:],errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(times[2,:],errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                    
        ax3.loglog(dts,errors[3,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax3.loglog(dts,errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax4.loglog(times[3,:],errors[3,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax4.loglog(times[4,:],errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax4.loglog(times[5,:],errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))



    times = np.zeros((len(methods_staggered),NN))
    errors = np.zeros((len(methods_staggered),NN))
    p_ave = np.zeros((len(methods_staggered),NN))
    p_std = np.zeros((len(methods_staggered),NN))

    
    fileName = folder+"/convergence_DeC_adaptive_"+dist+f"_{pr.name}.csv"
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile)
        for ir, row in enumerate(reader):
            if ir!=0:
                dt0 = float(row[0])
                if abs(dt0 - dts[ir-1])>1e-14:
                    raise ValueError("Wrong dt!!!!") 

                for im, method in enumerate(methods_staggered):
                    errors[im,ir-1] = row[1+4*im]                    
                    times[im,ir-1]  = row[2+4*im]                    
                    p_ave[im,ir-1]  = row[3+4*im]                  
                    p_std[im,ir-1]  = row[4+4*im]                    
            
            
            
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

        

#    ax1.set_title("DeC error convergence")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Error")
    ax1.set_xlabel(r"$\Delta t$")
    ax1.grid(True,which="both")
    fig1.set_tight_layout(True)
    fig1.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

#    ax2.set_title("DeC error convergence")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel("Error")
    ax2.set_xlabel("Computational time")
    ax2.grid(True,which="both")
    fig2.set_tight_layout(True)
    fig2.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")


 #   ax3.set_title("DeC error convergence")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylabel("Error")
    ax3.set_xlabel(r"$\Delta t$")
    ax3.grid(True,which="both")
    fig3.set_tight_layout(True)
    fig3.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

#    ax4.set_title("DeC error convergence")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_ylabel("Error")
    ax4.set_xlabel("Computational time")
    ax4.grid(True,which="both")
    fig4.set_tight_layout(True)
    fig4.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")


    plt.show(block=False)

    plt.close('all')

#%% 
#Test convergence of DeC for several orders
pr=ODEproblem("vibratingDamped")


colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


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

    print("There")


    # Compute and plot the errors 
    for order in range(3,10):
        fileName = folder+"/convergence_DeC_"+dist+f"_{pr.name}_order_{order}.csv"
        print("Evverywhere")

        errors = np.zeros((len(methods),len(dts)))
        times  = np.zeros((len(methods),len(dts)))

        with open(fileName) as myFile:
            reader = csv.reader(myFile)
            for ir, row in enumerate(reader):
                if ir!=0:
                    dt0 = float(row[0])
                    if abs(dt0 - dts[ir-1])>1e-14:
                        raise ValueError("Wrong dt!!!!") 
                    errors_all_methods = [float(z) for z in row[1:len(methods)+1]]
                    times_all_methods = [float(z) for z in row[len(methods)+1:2*len(methods)+1]]
                    for im, method in enumerate(methods):
                        times[im,ir-1] = times_all_methods[im]
                        errors[im,ir-1]= errors_all_methods[im] 
                        print(f"error {errors[im,ir-1]}")
                
               
        ax1.loglog(dts,errors[0,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax1.loglog(dts,errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax2.loglog(times[0,:],errors[0,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax2.loglog(times[1,:],errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(times[2,:],errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                    
        ax3.loglog(dts,errors[3,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax3.loglog(dts,errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax4.loglog(times[3,:],errors[3,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax4.loglog(times[4,:],errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax4.loglog(times[5,:],errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))



    times = np.zeros((len(methods_staggered),NN))
    errors = np.zeros((len(methods_staggered),NN))
    p_ave = np.zeros((len(methods_staggered),NN))
    p_std = np.zeros((len(methods_staggered),NN))

    
    fileName = folder+"/convergence_DeC_adaptive_"+dist+f"_{pr.name}.csv"
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile)
        for ir, row in enumerate(reader):
            if ir!=0:
                dt0 = float(row[0])
                if abs(dt0 - dts[ir-1])>1e-14:
                    raise ValueError("Wrong dt!!!!") 

                for im, method in enumerate(methods_staggered):
                    errors[im,ir-1] = row[1+4*im]                    
                    times[im,ir-1]  = row[2+4*im]                    
                    p_ave[im,ir-1]  = row[3+4*im]                  
                    p_std[im,ir-1]  = row[4+4*im]                    
            
            
            
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

        

#    ax1.set_title("DeC error convergence")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Error")
    ax1.set_xlabel(r"$\Delta t$")
    ax1.grid(True,which="both")
    fig1.set_tight_layout(True)
    fig1.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

#    ax2.set_title("DeC error convergence")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel("Error")
    ax2.set_xlabel("Computational time")
    ax2.grid(True,which="both")
    fig2.set_tight_layout(True)
    fig2.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")


 #   ax3.set_title("DeC error convergence")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylabel("Error")
    ax3.set_xlabel(r"$\Delta t$")
    ax3.grid(True,which="both")
    fig3.set_tight_layout(True)
    fig3.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

#    ax4.set_title("DeC error convergence")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_ylabel("Error")
    ax4.set_xlabel("Computational time")
    ax4.grid(True,which="both")
    fig4.set_tight_layout(True)
    fig4.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")


    plt.show(block=False)

    plt.close('all')


#Test convergence of DeC for several orders
pr=ODEproblem("threeBodies")


colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]



NN=5
dts=[pr.T_fin/2.0**k for k in range(4,4+NN)]


for dist in distributions:
    print(dist)
    errors=np.zeros((N_meth,len(dts)))
    times =np.zeros((N_meth,len(dts)))

    fig1,ax1 = plt.subplots(1,1)
    fig2,ax2 = plt.subplots(1,1)
    fig3,ax3 = plt.subplots(1,1)
    fig4,ax4 = plt.subplots(1,1)

    print("There")


    # Compute and plot the errors 
    for order in range(3,10):
        fileName = folder+"/convergence_DeC_"+dist+f"_{pr.name}_order_{order}.csv"
        print("Evverywhere")

        errors = np.zeros((len(methods),len(dts)))
        times  = np.zeros((len(methods),len(dts)))

        with open(fileName) as myFile:
            reader = csv.reader(myFile)
            for ir, row in enumerate(reader):
                if ir!=0:
                    dt0 = float(row[0])
                    if abs(dt0 - dts[ir-1])>1e-14:
                        raise ValueError("Wrong dt!!!!") 
                    errors_all_methods = [float(z) for z in row[1:len(methods)+1]]
                    times_all_methods = [float(z) for z in row[len(methods)+1:2*len(methods)+1]]
                    for im, method in enumerate(methods):
                        times[im,ir-1] = times_all_methods[im]
                        errors[im,ir-1]= errors_all_methods[im] 
                        print(f"error {errors[im,ir-1]}")
                
               
        ax1.loglog(dts,errors[0,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax1.loglog(dts,errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax2.loglog(times[0,:],errors[0,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax2.loglog(times[1,:],errors[1,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax2.loglog(times[2,:],errors[2,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

                    
        ax3.loglog(dts,errors[3,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax3.loglog(dts,errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax3.loglog(dts,[dt**(order)*errors[3,2]/dts[2]**(order) for dt in dts],":",color=colors[order])#,label="ref %d"%(order))

        ax4.loglog(times[3,:],errors[3,:],"-",color=colors[order],label="DeC(%d)"%(order))
        ax4.loglog(times[4,:],errors[4,:],"--",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        ax4.loglog(times[5,:],errors[5,:],"-.",color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))



    times = np.zeros((len(methods_staggered),NN))
    errors = np.zeros((len(methods_staggered),NN))
    p_ave = np.zeros((len(methods_staggered),NN))
    p_std = np.zeros((len(methods_staggered),NN))

    
    fileName = folder+"/convergence_DeC_adaptive_"+dist+f"_{pr.name}.csv"
    with open(fileName) as csvfile:
        reader = csv.reader(csvfile)
        for ir, row in enumerate(reader):
            if ir!=0:
                dt0 = float(row[0])
                if abs(dt0 - dts[ir-1])>1e-14:
                    raise ValueError("Wrong dt!!!!") 

                for im, method in enumerate(methods_staggered):
                    errors[im,ir-1] = row[1+4*im]                    
                    times[im,ir-1]  = row[2+4*im]                    
                    p_ave[im,ir-1]  = row[3+4*im]                  
                    p_std[im,ir-1]  = row[4+4*im]                    
            
            
            
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

        

#    ax1.set_title("DeC error convergence")
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax1.set_ylabel("Error")
    ax1.set_xlabel(r"$\Delta t$")
    ax1.grid(True,which="both")
    fig1.set_tight_layout(True)
    fig1.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

#    ax2.set_title("DeC error convergence")
    ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax2.set_ylabel("Error")
    ax2.set_xlabel("Computational time")
    ax2.grid(True,which="both")
    fig2.set_tight_layout(True)
    fig2.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_{pr.name}.pdf")


 #   ax3.set_title("DeC error convergence")
    ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax3.set_ylabel("Error")
    ax3.set_xlabel(r"$\Delta t$")
    ax3.grid(True,which="both")
    fig3.set_tight_layout(True)
    fig3.savefig(folder+f"/convergence_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")

#    ax4.set_title("DeC error convergence")
    ax4.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax4.set_ylabel("Error")
    ax4.set_xlabel("Computational time")
    ax4.grid(True,which="both")
    fig4.set_tight_layout(True)
    fig4.savefig(folder+f"/convergence_vs_time_DeC_adaptive_"+dist+f"_staggered_small_sub_{pr.name}.pdf")


    plt.show(block=False)

    plt.close('all')


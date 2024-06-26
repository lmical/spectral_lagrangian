#%%
# We need a couple of packages in this chapter
import numpy as np  
# This is the basic package in python with all the numerical functions

import matplotlib.pyplot as plt 
# This package allows to  plot

from timeit import timeit
import csv
from src import RungeKutta
from src import ODEproblems
from src.ADER import ADER, ADER_L2, ADER_u,ClassicalADER
from src.ODEproblems import ODEproblem
from src.utils import compute_integral_error, compute_integral_error_approx
import matplotlib

matplotlib.rcParams.update({'font.size': 13})

import os

folder = "figures/adaptive"

try: 
    os.mkdir(folder) 
except OSError as error: 
    print(error)  

methods = [ADER, ADER_L2, ADER_u, ClassicalADER]
methods_string = ["ADER","ADERdu", "ADERu", "Classic ADER"]
linestyles_methods = ["-","-.","--",":"]
markers_methods = ["o","x","+","."]

methods_string_ref = ["ADER","ADERdu", "ADERu", "cADER","ref"]
linestyles_methods_ref = ["-","-.","--",":",(0,(1,10))]

methods_staggered = [ADER_L2, ADER_u]
linestyles_methods_staggered = ["-.","--"]
markers_methods_staggered  = ["x","+"]

methods_improved = [ADER, ADER_L2, ADER_u]
linestyles_methods_improved = ["-","-.","--"]
markers_methods_improved = ["o","x","+"]

N_meth = len(methods)
N_meth_staggered = len(methods_staggered)

N_meth_improved = len(methods_improved)

distributions=["equispaced","gaussLobatto","gaussLegendre"]

colors=["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",\
       "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


adaptive_tolerance=1e-8
adaptive_max_order = 16


dist_names = {"equispaced":"Equispaced",
              "gaussLobatto":"Gauss-Lobatto",
              "gaussLegendre":"Gauss-Legendre"}

problems = ["linear_system2", "vibratingDamped", "threeBodies", "C5"]
N0s = [2, 2, 4 , 3]

#%% 
#Test convergence of DeC for several orders

#which type of plot 4 5 6 # error  # error vs time # speed up
iiii=6

for i_prob , prob in enumerate(problems):
    N0 = N0s[i_prob]

    pr=ODEproblem(prob)

    NN=5
    dts=[pr.T_fin/2.0**k for k in range(N0,N0+NN)]


    for dist in distributions:
        print(dist)
        errors=np.zeros((N_meth,len(dts)))
        times =np.zeros((N_meth,len(dts)))

        # fig1,ax1 = plt.subplots(1,1) # error
        # fig2,ax2 = plt.subplots(1,1) # error vs time
        # fig3,ax3 = plt.subplots(1,1) # speed-up

        

        if iiii==4:
            fig4,ax4 = plt.subplots(1,1) # error
        elif iiii==5:
            fig5,ax5 = plt.subplots(1,1) # error vs time
        elif iiii==6:
            fig6,ax6 = plt.subplots(1,1) # speed-up

        # fig7,ax7 = plt.subplots(1,1) # error
        # fig8,ax8 = plt.subplots(1,1) # error vs time
        # fig9,ax9 = plt.subplots(1,1) # speed-up

        

        # Compute and plot the errors 
        for order in range(3,10):
            fileName = folder+"/convergence_ADER_"+dist+f"_{pr.name}_order_{order}.csv"

            errors = np.zeros((len(methods),len(dts)))
            times  = np.zeros((len(methods),len(dts)))

            with open(fileName) as myFile:
                reader = csv.reader(myFile)
                for ir, row in enumerate(reader):
                    print(row)
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
                    
                
        #     ax1.loglog(dts,errors[0,:],\
        #             linestyle = linestyles_methods[0],\
        #             color=colors[order],label="ADER(%d)"%(order))
        #     ax1.loglog(dts,errors[1,:],\
        #             linestyle = linestyles_methods[1],\
        #             color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #     ax1.loglog(dts,errors[2,:],\
        #             linestyle = linestyles_methods[2],\
        #             color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #     ax1.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],\
        #             linestyle=(0,(1,10)),color=colors[order])#,label="ref %d"%(order))

        #     ax2.loglog(times[0,:],errors[0,:],\
        #             linestyle = linestyles_methods[0],\
        #             color=colors[order],label="ADER(%d)"%(order))
        #     ax2.loglog(times[1,:],errors[1,:],\
        #             linestyle = linestyles_methods[1],\
        #             color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #     ax2.loglog(times[2,:],errors[2,:],\
        #             linestyle = linestyles_methods[2],\
        #             color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        # #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

        #     ax3.semilogx(dts,times[0,:]/times[1,:],\
        #             linestyle = linestyles_methods[1],\
        #             marker = markers_methods[1],\
        #             color=colors[order])#,label="ADERdu(%d)"%(order))
        #     ax3.semilogx(dts,times[0,:]/times[2,:],\
        #             linestyle = linestyles_methods[2],\
        #             marker = markers_methods[2],\
        #             color=colors[order])#,label="ADERu(%d)"%(order))
        #     ax3.semilogx([], [],color=colors[order], label="ADER(%d)"%(order))
        # #   
            if iiii==4:    
                ax4.loglog(dts,errors[0,:],\
                        linestyle = linestyles_methods[0],\
                        color=colors[order],label="ADER(%d)"%(order))
                ax4.loglog(dts,errors[1,:],\
                        linestyle = linestyles_methods[1],\
                        color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
                ax4.loglog(dts,errors[2,:],\
                        linestyle = linestyles_methods[2],\
                        color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
                ax4.loglog(dts,errors[3,:],\
                        linestyle = linestyles_methods[3],\
                        color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
                ax4.loglog(dts,[max(dt**(order)*errors[0,2]/dts[2]**(order),1e-15) for dt in dts],\
                        linestyle=(0,(1,10)),color=colors[order])#,label="ref %d"%(order))

            if iiii==5:
                ax5.loglog(times[0,:],errors[0,:],\
                        linestyle = linestyles_methods[0],\
                            color=colors[order],label="ADER(%d)"%(order))
                ax5.loglog(times[1,:],errors[1,:],\
                        linestyle = linestyles_methods[1],\
                            color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
                ax5.loglog(times[2,:],errors[2,:],\
                        linestyle = linestyles_methods[2],\
                            color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
                ax5.loglog(times[3,:],errors[3,:],\
                        linestyle = linestyles_methods[3],\
                            color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            #    ax5.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

            if iiii==6:
                ax6.semilogx(dts,times[-1,:]/times[0,:],\
                        linestyle = linestyles_methods[0],\
                        marker = markers_methods[0],\
                        color=colors[order],label="ADER(%d)"%(order))
                ax6.semilogx(dts,times[-1,:]/times[1,:],\
                        linestyle = linestyles_methods[1],\
                        marker = markers_methods[1],\
                        color=colors[order])#,label="ADERdu(%d)"%(order))
                ax6.semilogx(dts,times[-1,:]/times[2,:],\
                        linestyle = linestyles_methods[2],\
                        marker = markers_methods[2],\
                        color=colors[order])#,label="ADERu(%d)"%(order))
        #   

                
        #     ax7.loglog(dts,errors[0,:],\
        #             linestyle = linestyles_methods[0],\
        #             color=colors[order],label="ADER(%d)"%(order))
        #     ax7.loglog(dts,errors[3,:],\
        #             linestyle = linestyles_methods[3],\
        #             color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        #     ax7.loglog(dts,[dt**(order)*errors[0,2]/dts[2]**(order) for dt in dts],\
        #             linestyle=(0,(1,10)),color=colors[order])#,label="ref %d"%(order))

        #     ax8.loglog(times[0,:],errors[0,:],\
        #             linestyle = linestyles_methods[0],\
        #             color=colors[order],label="ADER(%d)"%(order))
        #     ax8.loglog(times[3,:],errors[3,:],\
        #             linestyle = linestyles_methods[3],\
        #             color=colors[order])#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        # #    ax8.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

        #     ax9.semilogx(dts,times[-1,:]/times[0,:],\
        #             linestyle = linestyles_methods[0],\
        #             marker = markers_methods[0],\
        #             color=colors[order],label="ADER(%d)"%(order))
        # #   



        times = np.zeros((len(methods_staggered),NN))
        errors = np.zeros((len(methods_staggered),NN))
        p_ave = np.zeros((len(methods_staggered),NN))
        p_std = np.zeros((len(methods_staggered),NN))

        
        fileName = folder+"/convergence_ADER_adaptive_"+dist+f"_{pr.name}.csv"
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
                
                
                
    #     ax1.loglog(dts,errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #     ax1.loglog(dts,errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))

    #     ax2.loglog(times[0,:],errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #     ax2.loglog(times[1,:],errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    # #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

        if iiii==4:        
            ax4.loglog(dts,errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax4.loglog(dts,errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
        if iiii==5:
            ax5.loglog(times[0,:],errors[0,:],"--",color='black',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
            ax5.loglog(times[1,:],errors[1,:],"-.",color='grey',linewidth=3)#,label="DeCStag(%d,%d)"%(M_sub,K_iter))
    #    ax2.loglog(dts,[dt**(order)*errorsDeC[2]/dts[2]**(order) for dt in dts],":",label="ref %d"%(order))

        

    # #    ax1.set_title("DeC error convergence")
    #     index_set = [0,1,2,4]
    #     ax1.set_ylabel("Error")
    #     ax1.set_xlabel(r"$\Delta t$")
    #     ax1.grid(True,which="both")
    #     leg1=plt.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    #     lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
    #             linestyle=linestyles_methods_ref[i]) \
    #             for i in index_set]
    #     leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
    #                       loc='center left', bbox_to_anchor=(1, 0.3))
    #     ax1.add_artist(leg1)
    #     ax1.add_artist(leg2)
    #     fig1.set_tight_layout(True)
    #     fig1.savefig(folder+f"/convergence_ADER_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

    # #    ax2.set_title("DeC error convergence")
    #     index_set = [0,1,2]
    #     leg1=ax2.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    #     ax2.set_ylabel("Error")
    #     ax2.set_xlabel("Computational time")
    #     ax2.grid(True,which="both")
    #     lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
    #             linestyle=linestyles_methods_ref[i]) \
    #             for i in index_set]
    #     leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
    #                       loc='center left', bbox_to_anchor=(1, 0.3))
    #     ax2.add_artist(leg1)
    #     ax2.add_artist(leg2)
    #     fig2.set_tight_layout(True)
    #     fig2.savefig(folder+f"/convergence_vs_time_ADER_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

    # #    ax2.set_title("DeC error convergence")
    #     index_set =  [1,2]
    #     leg1 = ax3.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    #     ax3.set_ylabel("Speed-up")
    #     ax3.set_xlabel(r"$\Delta t$")
    #     lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
    #             linestyle=linestyles_methods_ref[i]) \
    #             for i in index_set]
    #     leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
    #                       loc='center left', bbox_to_anchor=(1, 0.3))
    #     ax3.add_artist(leg1)
    #     ax3.add_artist(leg2)
    #     fig3.set_tight_layout(True)
    #     fig3.savefig(folder+f"/speed_up_ADER_"+dist+f"_staggered_{pr.name}.pdf")


        if iiii==4:
            index_set = [3,0,2,1,4]#[0,1,2,3,4]
            leg1 = ax4.legend(loc='center left', bbox_to_anchor=(1, 0.7))
            ax4.set_ylabel("Error")
            ax4.set_xlabel(r"$\Delta t$")
            ax4.grid(True,which="both")
            lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
                    linestyle=linestyles_methods_ref[i]) \
                    for i in index_set]
            leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
                            loc='center left', bbox_to_anchor=(1, 0.2))
            ax4.add_artist(leg1)
            ax4.add_artist(leg2)
            ax4.set_title(dist_names[dist])
            fig4.set_tight_layout(True)
            fig4.savefig(folder+f"/convergence_ADER_wc_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

        if iiii==5:
    #    ax2.set_title("DeC error convergence")
            index_set = [3,0,2,1]# [0,1,2,3]
            leg1=ax5.legend(loc='center left', bbox_to_anchor=(1, 0.7))
            ax5.set_ylabel("Error")
            ax5.set_xlabel("Computational time")
            ax5.grid(True,which="both")
            lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
                    linestyle=linestyles_methods_ref[i]) \
                    for i in index_set]
            leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
                            loc='center left', bbox_to_anchor=(1, 0.2))
            ax5.add_artist(leg1)
            ax5.add_artist(leg2)
            ax5.set_title(dist_names[dist])
            fig5.set_tight_layout(True)
            fig5.savefig(folder+f"/convergence_vs_time_ADER_wc_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

        if iiii==6:
        #    ax2.set_title("DeC error convergence")
            index_set = [0,2,1]#[0,1,2]
            leg1=ax6.legend(loc='center left', bbox_to_anchor=(1, 0.7))
            ax6.set_ylabel("Speed-up")
            ax6.set_xlabel(r"$\Delta t$")
            #ax6.grid(True,which="both")        
            lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
                    linestyle=linestyles_methods_ref[i], marker = markers_methods[i]) \
                    for i in index_set]
            leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
                            loc='center left', bbox_to_anchor=(1, 0.2))
            ax6.add_artist(leg1)
            ax6.add_artist(leg2)
            ax6.set_ylim([0.9,2.1])
            ax6.grid(True)
            ax6.set_title(dist_names[dist])
            fig6.set_tight_layout(True)
            fig6.savefig(folder+f"/speed_up_ADER_wc_"+dist+f"_staggered_{pr.name}.pdf")


    # #    ax1.set_title("DeC error convergence")
    #     index_set = [0,3,4]
    #     leg1 = ax7.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    #     ax7.set_ylabel("Error")
    #     ax7.set_xlabel(r"$\Delta t$")
    #     ax7.grid(True,which="both")       
    #     lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
    #             linestyle=linestyles_methods_ref[i]) \
    #             for i in index_set]
    #     leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
    #                       loc='center left', bbox_to_anchor=(1, 0.3))
    #     ax7.add_artist(leg1)
    #     ax7.add_artist(leg2)        
    #     fig7.set_tight_layout(True)
    #     fig7.savefig(folder+f"/convergence_ADER_oc_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

    # #    ax2.set_title("DeC error convergence")
    #     index_set = [0,3]
    #     ax8.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    #     ax8.set_ylabel("Error")
    #     ax8.set_xlabel("Computational time")
    #     ax8.grid(True,which="both")      
    #     lines = [matplotlib.lines.Line2D([0], [0], color='black', linewidth=1,\
    #             linestyle=linestyles_methods_ref[i]) \
    #             for i in index_set]
    #     leg2 = plt.legend(lines, [methods_string_ref[i] for i in index_set],\
    #                       loc='center left', bbox_to_anchor=(1, 0.3))
    #     ax8.add_artist(leg1)
    #     ax8.add_artist(leg2)        
    #     fig8.set_tight_layout(True)
    #     fig8.savefig(folder+f"/convergence_vs_time_ADER_oc_adaptive_"+dist+f"_staggered_{pr.name}.pdf")

    # #    ax2.set_title("DeC error convergence")
    #     ax9.legend(loc='center left', bbox_to_anchor=(1, 0.7))
    #     ax9.set_ylabel("Speed-up")
    #     ax9.set_xlabel(r"$\Delta t$")
    #     #ax6.grid(True,which="both")
    #     fig9.set_tight_layout(True)
    #     fig9.savefig(folder+f"/speed_up_ADER_oc_"+dist+f"_staggered_{pr.name}.pdf")




        plt.show(block=False)

        plt.close('all')
        plt.cla()

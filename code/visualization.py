import numpy as np
import matplotlib.pyplot as plt

#==============================================================
# Function to print some infos
# Timestep
# Time
# Max and Min of the variables
#==============================================================
def printing_function(indt,t,H_field,v_field):
    print()
    print("------------------------------------------")
    print("Timestep:",indt)
    print("Time:",       t)
    print("Maxima H and v:",np.max(H_field),np.max(v_field))
    print("Minima H and v:",np.min(H_field),np.min(v_field))
    print("------------------------------------------")
    print()
#==============================================================
#
#
#
#==============================================================
# Function to plot H_field, v_field
#==============================================================
def plotting_function(indt,t,x_H,H_field,B_field,x_v,v_field,H_in_x_v,DATA,storing_info):
    N_el, N_local_nodes_H = x_H.shape
    degree_H=N_local_nodes_H-1
    degree_v=degree_H+1

    fig, axs = plt.subplots(2,2) #Array of subplots
    fig.suptitle(DATA.test)

    #H
    for inde in range(N_el):
        axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker="*")
    axs[0,0].set_title("H")
    axs[0,0].set_xlabel("x")
    axs[0,0].grid()
    # axs[0,0].set_ylabel("y")
    # axs[0,0].legend("H")

    #v
    axs[0,1].plot(x_v,v_field)
    axs[0,1].set_title("v")
    axs[0,1].set_xlabel("x")
    axs[0,1].grid()

    #eta
    for inde in range(N_el):
        axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker="*")
    axs[1,0].set_title("\eta")
    axs[1,0].set_xlabel("x")
    axs[1,0].grid()

    #q
    axs[1,1].plot(x_v,v_field*H_in_x_v)
    axs[1,1].set_title("q")
    axs[1,1].set_xlabel("x")
    axs[1,1].grid()

    fig.tight_layout()

    if DATA.storing==True and storing_info==True:
        plt.savefig(DATA.folder+"/"+DATA.test+"/"+DATA.test+"_perturbation"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+str(N_el)+"_CFL_"+str(DATA.CFL)+".pdf", format="pdf", bbox_inches="tight")


    plt.show()

#==============================================================
#
#
#
#==============================================================
# Function to store in a file
#==============================================================
def storing(H_field, v_field, x_v, B_field, H_in_x_v, B_in_x_v, M_Local_to_Global, DATA):
    f=open(DATA.folder+"/"+DATA.test+"/values_perturbation"+str(DATA.perturbation)+"_LxF"+str(DATA.LaxFriedrichs)+"_P"+str(DATA.order_space-1)+"P"+str(DATA.order_space)+"_"+"{:05d}".format(DATA.N_el)+"_CFL_"+str(DATA.CFL)+".dat","w+")
    f.write(" indi,                  x_v,                   v,                   H,                   q,                 eta\n")
    for indi in range(len(v_field)):
        towrite=format(indi, '5d')+",     "+format(x_v[indi], '.14f')+",    "+format(v_field[indi], '.14f')+",    "+format(H_in_x_v[indi], '.14f')+",    "+format(H_in_x_v[indi]*v_field[indi], '.14f')+",    "+format(H_in_x_v[indi]+B_in_x_v[indi], '.14f')+",\n"
        f.write(towrite)
    f.close()

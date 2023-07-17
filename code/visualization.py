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
def plotting_function(indt,t,x_H,H_field,x_v,v_field):
    N_el, N_local_nodes_H = x_H.shape

    #H
    plt.figure()
    plt.title("H")
    for inde in range(N_el):
        plt.plot(x_H[inde,:],H_field[inde,:], marker="*")
    #plt.legend()
    plt.xlabel("x")
    #plt.ylabel("y")
    plt.grid()
    plt.show()    

    #v
    plt.figure()
    plt.title("v")
    plt.plot(x_v,v_field, marker="*")
    #plt.legend()
    plt.xlabel("x")
    #plt.ylabel("y")
    plt.grid()
    plt.show()    
#==============================================================
#
#
#
#==============================================================
# Function to store in a file
#==============================================================
def storing(H_field, v_field, x_v, B_field, local_values_v_in_H, M_Local_to_Global, DATA):
    f=open(DATA.folder+"/"+DATA.test+"/"+DATA.test+"_"+"P"+str(DATA.order_space-1)+"P"+str(DATA.order_space)+"_"+str(DATA.N_el)+"_CFL_"+str(DATA.CFL)+".dat","w+")
    f.write("  indi    x_v    v\n")
    for indi in range(len(v_field)):
        towrite=" "+str(indi+1)+"     "+format(x_v[indi], '.14f')+"    "+format(v_field[indi], '.14f')+"\n"
        f.write(towrite)
    f.close()

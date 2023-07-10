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
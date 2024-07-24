import numpy as np
import matplotlib.pyplot as plt
import test_dependent
import lagrangian_scheme

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
def plotting_function(indt,x_H,H_field,B_field,x_v,v_field,H_in_x_v,DATA,storing_info):
    N_el, N_local_nodes_H = x_H.shape
    degree_H=N_local_nodes_H-1
    degree_v=degree_H+1

    fig, axs = plt.subplots(2,2) #Array of subplots
    fig.suptitle(DATA.test)

    #H
    for inde in range(N_el):
        axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker=".")
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
        axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker=".")
    axs[1,0].set_title("\eta")
    axs[1,0].set_xlabel("x")
    axs[1,0].grid()

    #q
    axs[1,1].plot(x_v,v_field*H_in_x_v)
    axs[1,1].set_title("q")
    axs[1,1].set_xlabel("x")
    axs[1,1].grid()

    fig.tight_layout()

    #SAVE
    if DATA.storing==True and storing_info==True:
        if DATA.LaxFriedrichs=="Active" or DATA.LaxFriedrichs=="Disabled":
            plt.savefig(DATA.folder+"/"+DATA.test+"/pic_values_pert"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".pdf", format="pdf", bbox_inches="tight")
        elif DATA.LaxFriedrichs=="ShockDetector_divV" or DATA.LaxFriedrichs=="ShockDetector_divV_tn" or DATA.LaxFriedrichs=="ShockDetector_divV_tn_memory":
            plt.savefig(DATA.folder+"/"+DATA.test+"/pic_values_pert"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_K"+str(DATA.K_limiter_divV)+"_NLimitedNeighbours"+str(DATA.N_limited_neighbours)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".pdf", format="pdf", bbox_inches="tight")
        else:
            print("Invalid LxF option in plotting_function")
            print(DATA.LaxFriedrichs)
            quit()

    # plt.show()
#==============================================================
#
#
#
#==============================================================
# Function to plot Shock Detector
#==============================================================
def plotting_ShockDetector(indt,x_H,H_field,B_field,x_v,v_field,H_in_x_v,M_Local_to_Global,w_H,local_derivatives_v_in_H,TroubledCells_memory,DATA,storing_info):


    TroubledCells=lagrangian_scheme.ShockDetector(v_field,M_Local_to_Global,H_field,x_v,w_H,local_derivatives_v_in_H,DATA)
    if DATA.LaxFriedrichs=="ShockDetector_divV_tn_memory":
        TroubledCells=lagrangian_scheme.MergeTroubledCells(TroubledCells,TroubledCells_memory,DATA)


    N_el, N_local_nodes_H = x_H.shape
    degree_H=N_local_nodes_H-1
    degree_v=degree_H+1

    fig, axs = plt.subplots(2,2) #Array of subplots
    fig.suptitle(DATA.test)


    if degree_H==0:
        marker_choice="." #Otherwise H is not plotted
    else:
        marker_choice=" "

    #H
    for inde in range(N_el):
        if TroubledCells[inde]==0:
            axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker=marker_choice,color="k")
        elif TroubledCells[inde]==1:
            axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker=marker_choice,color="r")
        elif TroubledCells[inde]>1:
            axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker=marker_choice,color='orange')
        elif TroubledCells[inde]==-1:
            axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker=marker_choice,color='blue')
        elif TroubledCells[inde]<-1:
            axs[0,0].plot(x_H[inde,:],H_field[inde,:], marker=marker_choice,color='cyan')
        else:
            print("Invalid adimssible value for TroubledCells")
            quit()
    axs[0,0].set_title("H")
    axs[0,0].set_xlabel("x")
    axs[0,0].grid()
    # axs[0,0].set_ylabel("y")
    # axs[0,0].legend("H")

    #v
    for inde in range(N_el):
        global_indices_v=M_Local_to_Global[inde,:]
        x_v_local=x_v[global_indices_v]
        v_local=v_field[global_indices_v]
        if TroubledCells[inde]==0:
            axs[0,1].plot(x_v_local,v_local, marker=marker_choice,color="k")
        elif TroubledCells[inde]==1:
            axs[0,1].plot(x_v_local,v_local, marker=marker_choice,color="r")
        elif TroubledCells[inde]>1:
            axs[0,1].plot(x_v_local,v_local, marker=marker_choice,color='orange')
        elif TroubledCells[inde]==-1:
            axs[0,1].plot(x_v_local,v_local, marker=marker_choice,color='blue')
        elif TroubledCells[inde]<-1:
            axs[0,1].plot(x_v_local,v_local, marker=marker_choice,color='cyan')
        else:
            print("Invalid adimssible value for TroubledCells")
            quit()
    axs[0,1].set_title("v")
    axs[0,1].set_xlabel("x")
    axs[0,1].grid()

    #eta
    for inde in range(N_el):
        if TroubledCells[inde]==0:
            axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker=marker_choice,color="k")
        elif TroubledCells[inde]==1:
            axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker=marker_choice,color="r")
        elif TroubledCells[inde]>1:
            axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker=marker_choice,color='orange')
        elif TroubledCells[inde]==-1:
            axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker=marker_choice,color='blue')
        elif TroubledCells[inde]<-1:
            axs[1,0].plot(x_H[inde,:],H_field[inde,:]+B_field[inde,:], marker=marker_choice,color='cyan')
        else:
            print("Invalid adimssible value for TroubledCells")
            quit()
    axs[1,0].set_title("\eta")
    axs[1,0].set_xlabel("x")
    axs[1,0].grid()

    #q
    q_field=v_field*H_in_x_v
    for inde in range(N_el):
        global_indices_v=M_Local_to_Global[inde,:]
        x_v_local=x_v[global_indices_v]
        q_local=q_field[global_indices_v]
        if TroubledCells[inde]==0:
            axs[1,1].plot(x_v_local,q_local, marker=marker_choice,color="k")
        elif TroubledCells[inde]==1:
            axs[1,1].plot(x_v_local,q_local, marker=marker_choice,color="r")
        elif TroubledCells[inde]>1:
            axs[1,1].plot(x_v_local,q_local, marker=marker_choice,color='orange')
        elif TroubledCells[inde]==-1:
            axs[1,1].plot(x_v_local,q_local, marker=marker_choice,color='blue')
        elif TroubledCells[inde]<-1:
            axs[1,1].plot(x_v_local,q_local, marker=marker_choice,color='cyan')
        else:
            print("Invalid adimssible value for TroubledCells")
            quit()
    axs[1,1].set_title("q")
    axs[1,1].set_xlabel("x")
    axs[1,1].grid()

    fig.tight_layout()

    #SAVE
    if DATA.storing==True and storing_info==True:
        if DATA.LaxFriedrichs=="Active" or DATA.LaxFriedrichs=="Disabled":
            plt.savefig(DATA.folder+"/"+DATA.test+"/pic_ShockDetector_pert"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".pdf", format="pdf", bbox_inches="tight")
        elif DATA.LaxFriedrichs=="ShockDetector_divV" or DATA.LaxFriedrichs=="ShockDetector_divV_tn" or DATA.LaxFriedrichs=="ShockDetector_divV_tn_memory":
            plt.savefig(DATA.folder+"/"+DATA.test+"/pic_ShockDetector_pert"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_K"+str(DATA.K_limiter_divV)+"_NLimitedNeighbours"+str(DATA.N_limited_neighbours)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".pdf", format="pdf", bbox_inches="tight")
        else:
            print("Invalid LxF option in ShockDetector")
            print(DATA.LaxFriedrichs)
            quit()



    # plt.show()
#==============================================================
#
#
#
#==============================================================
# Function to store in a file
# indi, x_v, v, H, q, eta
#==============================================================
def storing(H_field, v_field, x_v, B_field, H_in_x_v, B_in_x_v, M_Local_to_Global, DATA):

    #SAVE
    if DATA.LaxFriedrichs=="Active" or DATA.LaxFriedrichs=="Disabled":
        f=open(DATA.folder+"/"+DATA.test+"/values_pert"+str(DATA.perturbation)+"_"+"P"+str(DATA.order_space-1)+"P"+str(DATA.order_space)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".dat","w+")
    elif DATA.LaxFriedrichs=="ShockDetector_divV" or DATA.LaxFriedrichs=="ShockDetector_divV_tn" or DATA.LaxFriedrichs=="ShockDetector_divV_tn_memory":
        f=open(DATA.folder+"/"+DATA.test+"/values_pert"+str(DATA.perturbation)+"_"+"P"+str(DATA.order_space-1)+"P"+str(DATA.order_space)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_K"+str(DATA.K_limiter_divV)+"_NLimitedNeighbours"+str(DATA.N_limited_neighbours)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".dat","w+")
    else:
        print("Invalid LxF option in storing")
        print(DATA.LaxFriedrichs)
        quit()



    f.write(" indi,                  x_v,                   v,                   H,                   q,                 eta\n")
    for indi in range(len(v_field)):
        towrite=format(indi, '5d')+",     "+format(x_v[indi], '.14f')+",    "+format(v_field[indi], '.14f')+",    "+format(H_in_x_v[indi], '.14f')+",    "+format(H_in_x_v[indi]*v_field[indi], '.14f')+",    "+format(H_in_x_v[indi]+B_in_x_v[indi], '.14f')+",\n"
        f.write(towrite)
    f.close()
#==============================================================
#
#
#
#==============================================================
# Function to compute the error
# v, H, q
#==============================================================
def compute_error(H_field, v_field, x_v, x_H, H_in_x_v, M_Local_to_Global, w_H, w_v, local_derivatives_v_in_H, local_derivatives_v, DATA):

    #Again, just be sure
    if DATA.analytical_solution==False or DATA.perturbation!=0:
        print("You are trying to compute the error without analytical solution")
        quit()


    N_el, N_local_nodes_H = H_field.shape
    N_local_nodes_v=N_local_nodes_H+1
    N_global_nodes_v=len(x_v)


    #---------------------------------------
    # L^1 error
    #---------------------------------------

    #In local_derivatives_v we have
    #Rows basis functions
    #Columns x_j_v
    #I transpose it to have in each row a specifc x_i_v and in the columns the basis functions
    in_v_local_derivatives_v=local_derivatives_v.transpose()

    #Same for H
    in_H_local_derivatives_v=local_derivatives_v_in_H.transpose()

    #Initialization
    errorH=0
    errorv=0
    errorq=0


    for inde in range(N_el):

        global_indices_v=M_Local_to_Global[inde,:]
        x_v_local=x_v[global_indices_v]

        #H
        H_local=H_field[inde,:]
        for indi in range(N_local_nodes_H):
            H_ex   = test_dependent.Analytical_State(x_H[inde,indi],DATA.T,DATA)[0]
            errorH = errorH+np.abs(H_local[indi]-H_ex)*w_H[indi]*np.abs(np.sum(in_H_local_derivatives_v[indi,:]*x_v_local)) #NB: This last term is dx/dxi

        #v
        v_local=v_field[global_indices_v]
        for indi in range(N_local_nodes_v):
            v_ex   = test_dependent.Analytical_State(x_v_local[indi],DATA.T,DATA)[1]
            errorv = errorv+np.abs(v_local[indi]-v_ex)*w_v[indi]*np.abs(np.sum(in_v_local_derivatives_v[indi,:]*x_v_local)) #NB: This last term is dx/dxi

        #q
        H_local_in_v=H_in_x_v[global_indices_v]
        q_local=v_local*H_local_in_v
        for indi in range(N_local_nodes_v):
            H_ex, v_ex=test_dependent.Analytical_State(x_v_local[indi],DATA.T,DATA)
            q_ex   = H_ex*v_ex
            errorq = errorq+np.abs(q_local[indi]-q_ex)*w_v[indi]*np.abs(np.sum(in_v_local_derivatives_v[indi,:]*x_v_local)) #NB: This last term is dx/dxi
        


    #---------------------------------------
    # 2-error
    #---------------------------------------

    errorHbis = 0.
    errorvbis = 0.

    H_sol=np.zeros(H_field.shape)
    v_sol=np.zeros(len(v_field))

    for inde in range(N_el):
        for indi in range(N_local_nodes_H):
            H_sol[inde,indi] = test_dependent.Analytical_State(x_H[inde,indi],DATA.T,DATA)[0]

    for indi in range(N_global_nodes_v):
        v_sol[indi] = test_dependent.Analytical_State(x_v[indi],DATA.T,DATA)[1]



    errorHbis=np.linalg.norm(H_sol-H_field)/np.sqrt(H_field.size)
    errorvbis=np.linalg.norm(v_sol-v_field)/np.sqrt(len(v_field))

    #SAVE
    if DATA.LaxFriedrichs=="Active" or DATA.LaxFriedrichs=="Disabled":
        f=open(DATA.folder+"/"+DATA.test+"/errors_pert"+str(DATA.perturbation)+"_"+"P"+str(DATA.order_space-1)+"P"+str(DATA.order_space)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".dat","w+")
    elif DATA.LaxFriedrichs=="ShockDetector_divV" or DATA.LaxFriedrichs=="ShockDetector_divV_tn" or DATA.LaxFriedrichs=="ShockDetector_divV_tn_memory":
        f=open(DATA.folder+"/"+DATA.test+"/errors_pert"+str(DATA.perturbation)+"_"+"P"+str(DATA.order_space-1)+"P"+str(DATA.order_space)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_K"+str(DATA.K_limiter_divV)+"_NLimitedNeighbours"+str(DATA.N_limited_neighbours)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".dat","w+")
    else:
        print("Invalid LxF option in storing")
        print(DATA.LaxFriedrichs)
        quit()


    f.write(" N_el,       L^1(v),       L^1(H),       L^1(q),      ||v||_2,      ||H||_2,\n")
    towrite=format(N_el, '5d')+",    "+format(errorv, '.3e')+",    "+format(errorH, '.3e')+",    "+format(errorq, '.3e')+",    "+format(errorvbis, '.3e')+",    "+format(errorHbis, '.3e')+",\n"
    f.write(towrite)
    f.close()
#==============================================================
#
#
#
#==============================================================
# Function to plot the error
# v, H, q
#==============================================================
def plot_error(H_field, v_field, x_v, x_H, H_in_x_v, DATA):

    #Again, just be sure
    if DATA.analytical_solution==False:
        print("You are trying to plot the error without analytical solution")
        quit()


    N_el, N_local_nodes_H = H_field.shape
    N_local_nodes_v=N_local_nodes_H+1
    N_global_nodes_v=len(x_v)
    degree_H=N_local_nodes_H-1
    degree_v=N_local_nodes_v-1

    q_field=np.zeros(len(v_field))

    for indi in range(N_global_nodes_v):
        q_field[indi]=H_in_x_v[indi]*v_field[indi]

    H_sol=np.zeros(H_field.shape)
    v_sol=np.zeros(len(v_field))
    q_sol=np.zeros(len(v_field))


    for inde in range(N_el):
        for indi in range(N_local_nodes_H):
            H_sol[inde,indi] = test_dependent.Analytical_State(x_H[inde,indi],DATA.T,DATA)[0]

    for indi in range(N_global_nodes_v):
        H_ex, v_ex = test_dependent.Analytical_State(x_v[indi],DATA.T,DATA)
        v_sol[indi] = v_ex
        q_sol[indi] = H_ex*v_ex



    errorH=H_sol-H_field
    errorv=v_sol-v_field
    errorq=q_sol-q_field

    fig, axs = plt.subplots(1,3) #Array of subplots
    fig.suptitle(DATA.test+" Error")

    #H
    for inde in range(N_el):
        axs[0].plot(x_H[inde,:],errorH[inde,:], marker=".")
    axs[0].set_title("H")
    axs[0].set_xlabel("x")
    axs[0].grid()
    # axs[0].set_ylabel("y")
    # axs[0].legend("H")

    #v
    axs[1].plot(x_v,errorv)
    axs[1].set_title("v")
    axs[1].set_xlabel("x")
    axs[1].grid()


    #q
    axs[2].plot(x_v,errorq)
    axs[2].set_title("q")
    axs[2].set_xlabel("x")
    axs[2].grid()


    fig.tight_layout()

    if DATA.storing==True:
        if DATA.LaxFriedrichs=="Active" or DATA.LaxFriedrichs=="Disabled":
            plt.savefig(DATA.folder+"/"+DATA.test+"/pic_errors_pert"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".pdf", format="pdf", bbox_inches="tight")
        elif DATA.LaxFriedrichs=="ShockDetector_divV" or DATA.LaxFriedrichs=="ShockDetector_divV_tn" or DATA.LaxFriedrichs=="ShockDetector_divV_tn_memory":
            plt.savefig(DATA.folder+"/"+DATA.test+"/pic_errors_pert"+str(DATA.perturbation)+"_"+"P"+str(degree_H)+"P"+str(degree_v)+"_"+DATA.time_scheme+"_LxF"+str(DATA.LaxFriedrichs)+"_K"+str(DATA.K_limiter_divV)+"_NLimitedNeighbours"+str(DATA.N_limited_neighbours)+"_"+DATA.jump_CIP_in_v+"_jeta"+str(DATA.jump_eta_in_x)+"_CFL"+str(DATA.CFL)+"_N_el"+"{:05d}".format(DATA.N_el)+".pdf", format="pdf", bbox_inches="tight")
        else:
            print("Invalid LxF option in plot_error")
            print(DATA.LaxFriedrichs)
            quit()

    # plt.show()

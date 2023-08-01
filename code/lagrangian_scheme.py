import numpy as np
import Riemann_solver
import test_dependent

#==============================================================
# Obtaining Hhat_field at the beginning of the simulation for string mass conservation
#==============================================================
def get_Hhat_on_reference_element(H_field,x_v,local_derivatives_v_in_H,M_Local_to_Global):
    """
    Getting the field on the reference element for strong mass conservation
     Hhat_i = H_i(0)*det J(x_i,0)
     J(xi,0) = grad_xi x (xi,0) = sum_j x_j(0) grad_xi phi_j (xi)
    """
    Hhat_field = np.zeros(H_field.shape)

    N_el, local_nodes_H = H_field.shape

    for inde in range(N_el):
        global_indices_v=M_Local_to_Global[inde,:]
        x_v_local=x_v[global_indices_v]
        in_H_local_derivatives_v=local_derivatives_v_in_H.transpose()
        #----------------------------------
        # print(x_v_local)
        # print(local_derivatives_v_in_H)
        # print(in_H_local_derivatives_v)
        #----------------------------------
        for indi in range(local_nodes_H):
            # d_xi x(xi,0)
            d_xi_x=sum(in_H_local_derivatives_v[indi,:]*x_v_local)
            Hhat_field[inde,indi] = H_field[inde,indi]*d_xi_x

    #----------------------------------
    # print(H_field)
    # print(Hhat_field)
    #----------------------------------
    return Hhat_field
#==============================================================
#
#
#
#==============================================================
# Strong mass conservation
# Compute the H_field out of the Hhat_field
# NB: It is kind of the inverse as get_Hhat_on_reference_element
#==============================================================
def strong_mass_conservation(Hhat_field,x_v,local_derivatives_v_in_H,M_Local_to_Global):
    """
    Getting the H_field from Hhat_field
     H_i(t) = Hhat_i/det J(x_i,t)
     J(xi,t) = grad_xi x (xi,t) = sum_j x_j(t) grad_xi phi_j (xi)
    """
    H_field = np.zeros(Hhat_field.shape)

    N_el, local_nodes_H = Hhat_field.shape

    for inde in range(N_el):
        global_indices_v=M_Local_to_Global[inde,:]
        x_v_local=x_v[global_indices_v]
        in_H_local_derivatives_v=local_derivatives_v_in_H.transpose()
        #----------------------------------
        # print(x_v_local)
        # print(local_derivatives_v_in_H)
        # print(in_H_local_derivatives_v)
        #----------------------------------
        for indi in range(local_nodes_H):
            # d_xi x(xi,0)
            d_xi_x=sum(in_H_local_derivatives_v[indi,:]*x_v_local)
            H_field[inde,indi] = Hhat_field[inde,indi]/d_xi_x


    #----------------------------------
    # print(H_field)
    # print(Hhat_field)
    #----------------------------------
    return H_field
#==============================================================
#
#
#
#==============================================================
# Local computation of the lumped mass matrix
# NB: It is time-dependent because of the lumping
# NB: there must be a global assembling later
#==============================================================
def Lumped_Mass_Matrix_K(w_v,in_v_local_derivatives_v,x_v_local):
    """
    Local computation of the lumped mass matrix (indeed, there must be a global assembling later)
    Computation of
    w_i^K det J(x_i,t)|_Khat
    J(xi,t) = grad_xi x (xi,t) = sum_j x_j(t) grad_xi phi_j (xi)
    """
    N_local_nodes_v=len(w_v)
    M_v_in_K=np.zeros(N_local_nodes_v)
    for indi in range(N_local_nodes_v):
        #in_v_local_derivatives has
        #rows x_i_v
        #column basis functions
        d_xi_x=sum(in_v_local_derivatives_v[indi,:]*x_v_local)
        #print(d_xi_x)
        M_v_in_K[indi]=w_v[indi]*d_xi_x
    return M_v_in_K
#==============================================================
#
#
#
#==============================================================
# Computation of the lumped mass matrix
# NB: It is time-dependent because of the lumping
# NB: Uses Lumped_Mass_Matrix_K
#==============================================================
def Lumped_Mass_Matrix(w_v,x_v,M_Local_to_Global,local_derivatives_v,DATA):
    """
    With the assumed lumping, the mass matrix is diagonal and reads
    M_ii=sum_{K in K_i} w_i^K det J(x_i,t)|_Khat
    J(xi,t) = grad_xi x (xi,t) = sum_j x_j(t) grad_xi phi_j (xi)
    """

    #Initialize mass matrix
    M_v=np.zeros(len(x_v))

    N_el, N_local_nodes_v = M_Local_to_Global.shape

    #In local_derivatives_v we have
    #Rows basis functions
    #Columns x_j_v
    #I transpose it to have in each row a specifc x_i_v and in the columns the basis functions
    in_v_local_derivatives_v=local_derivatives_v.transpose()

    # print("Rows basis functions, columns DoFs")
    # print(local_derivatives_v) 
    # print("Columns basis functions, rows DoFs")
    # print(in_v_local_derivatives_v) 


    #-----------------------------------------------
    # #NB: This should be zero
    # for indi in range(N_local_nodes_v):
    #     print(sum(in_v_local_derivatives_v[indi,:]))
    #-----------------------------------------------

    for inde in range(N_el):
        global_indices_v=M_Local_to_Global[inde,:]
        x_v_local=x_v[global_indices_v]


        # Computation of
        # w_i^K det J(x_i,t)|_Khat
        # J(xi,t) = grad_xi x (xi,t) = sum_j x_j(t) grad_xi phi_j (xi)
        M_v_in_K=Lumped_Mass_Matrix_K(w_v,in_v_local_derivatives_v,x_v_local)


        #Assembling
        for indi in range(N_local_nodes_v):
            
            #-------------------------------------------------------
            #NB: This is what is in M_v_in_K[indi] 
            #d_xi_x=sum(in_v_local_derivatives_v[indi,:]*x_v_local)
            #w_v[indi]*d_xi_x
            #-------------------------------------------------------
            M_v[global_indices_v[indi]]=M_v[global_indices_v[indi]]+M_v_in_K[indi] 
    

    if DATA.periodic==True:
        M_v[0]=M_v[0]+M_v[-1]
        M_v[-1]=M_v[0]


    return M_v
#==============================================================
#
#
#
#==============================================================
# Local computation of the space residuals v
# NB: there must be a global assembling later
# NB: To be multiplied by g
#==============================================================
def Space_Residuals_v_K(w_v,in_v_local_derivatives_H,H_local,B_local):
    """
    Computation of
    w_i^K sum_{x_j_H in K} (H_j+B_j) grad_xi psi_j(xi_i)
    """

    N_local_nodes_v=len(w_v)
    phi_i_v_in_K=np.zeros(N_local_nodes_v)
    for indi in range(N_local_nodes_v):
        d_H_B=sum(in_v_local_derivatives_H[indi,:]*(H_local+B_local))
        phi_i_v_in_K[indi]=w_v[indi]*d_H_B
    
    return phi_i_v_in_K
#==============================================================
#
#
#
#==============================================================
# Space residuals v
# NB: Uses Space_residuals_v_K
# NB: To be multiplied by g
#==============================================================
def Space_Residuals_v(H_field, B_field, w_v, local_derivatives_H_in_v, M_Local_to_Global,DATA):
    """
    With the assumed lumping, the space residuals read
    phi_i=sum_{K in K_i} sum_{x_j_H in K} w_i^K grad_xi psi_j(xi_i)
    # NB: According to my computations, the Jacobian should go away
    # NB: To be multiplied by g
    # NB: One must add some boundary terms to guarantee coupling
    """

    N_local_nodes_H=H_field.shape[1]
    N_el, N_local_nodes_v = M_Local_to_Global.shape
    N_global_nodes_v=(N_local_nodes_v-1)*N_el+1

    #In local_derivatives_H_in_v we have
    #Rows basis functions
    #Columns x_j_v
    #I transpose it to have in each row a specifc x_i_v and in the columns the basis functions
    in_v_local_derivatives_H=local_derivatives_H_in_v.transpose()
    #----------------------------------------------
    # This should be 0
    # for indi in range(N_local_nodes_v):
    #     print(sum(in_v_local_derivatives_H[indi,:]))
    # print(in_v_local_derivatives_H)
    # quit()
    #----------------------------------------------


    phi_v=np.zeros(N_global_nodes_v)

    for inde in range(N_el):
        global_indices_v=M_Local_to_Global[inde,:]
 
        H_local=H_field[inde,:]
        B_local=B_field[inde,:]

        # Computation of
        # w_i^K sum_{x_j_H in K} (H_j+B_j) grad_xi psi_j(xi_i)

        phi_i_v_in_K=Space_Residuals_v_K(w_v,in_v_local_derivatives_H,H_local,B_local)

        #Assembling
        for indi in range(N_local_nodes_v):
            
            #-------------------------------------------------------
            #NB: This is what is in phi_i_v_in_K[indi] 
            #d_H_and_B=sum(in_v_local_derivatives_v[indi,:]*x_v_local)
            #w_v[indi]*d_H_and_B
            #-------------------------------------------------------
            phi_v[global_indices_v[indi]]=phi_v[global_indices_v[indi]]+phi_i_v_in_K[indi] 


    if DATA.periodic==True:
        phi_v[0]=phi_v[0]+phi_v[-1]
        phi_v[-1]=phi_v[0]


    return phi_v
#==============================================================
# Coupling terms in the space residuals v
# NB: To be multiplied by g
#==============================================================
def Coupling_Terms_Space_Residuals_v(H_field, B_field, v_field, M_Local_to_Global, M_faces, x_v, DATA):
    """
    # Coupling boundary terms to guarantee coupling
    CT_phi_i^K=int_{partial K} (etahat-eta|_K) 
    """

    N_local_nodes_H=H_field.shape[1]
    N_el, N_local_nodes_v = M_Local_to_Global.shape
    N_global_nodes_v=(N_local_nodes_v-1)*N_el+1

    
    CT_phi_v=np.zeros(N_global_nodes_v)

    #Loop on the internal elements
    for inde in range(N_el):

        global_indices_v=M_Local_to_Global[inde,:]


        #Riemann Problem Left
        #outside  L
        #inside   R
        if inde==0: #Care for the BC
            H_outside, B_outside, v_outside = test_dependent.BC_state(DATA, x_v[0], H_field[0], B_field[0], v_field[0], H_field[-1,-1], B_field[-1,-1], v_field[-1],"L")
            H_inside  = H_field[inde,0]                            #First from the current cell
            B_inside  = B_field[inde,0]                            #First from the current cell
            v_inside    = v_field[global_indices_v[0]]             #v_continuous from first node of the current cell
            q_outside = np.array([H_outside,H_outside*v_outside])  #L state
            q_inside  = np.array([H_inside ,H_inside *v_inside])   #R state
            if DATA.WB==True:
                Hhat, HUhat, Bhat=Riemann_solver.shallow_water_hll_WB(q_outside, q_inside, B_outside, B_inside, DATA.g)
            else:
                Hhat, HUhat = Riemann_solver.shallow_water_hll(q_outside, q_inside, DATA.g)
                Bhat = test_dependent.Bathymetry(x_v[global_indices_v[0]],DATA)
            CT_phi_v[global_indices_v[0]]=CT_phi_v[global_indices_v[0]]- ( (Hhat+Bhat)-(H_inside+B_inside) ) #IMPORTANT: - sign because of normal
        else:
            H_outside = H_field[inde-1,-1]                       #Last from the left cell
            B_outside = B_field[inde-1,-1]                       #Last from the left cell
            H_inside  = H_field[inde,0]                          #First from the current cell
            B_inside  = B_field[inde,0]                          #First from the current cell
            v_both    = v_field[global_indices_v[0]]             #v_continuous from first node of the current cell
            q_outside = np.array([H_outside,H_outside*v_both])   #L state
            q_inside  = np.array([H_inside ,H_inside *v_both])   #R state
            if DATA.WB==True:
                Hhat, HUhat, Bhat=Riemann_solver.shallow_water_hll_WB(q_outside, q_inside, B_outside, B_inside, DATA.g)
            else:
                Hhat, HUhat = Riemann_solver.shallow_water_hll(q_outside, q_inside, DATA.g)
                Bhat = test_dependent.Bathymetry(x_v[global_indices_v[0]],DATA)

            #Add contribution
            CT_phi_v[global_indices_v[0]]=CT_phi_v[global_indices_v[0]]- ( (Hhat+Bhat)-(H_inside+B_inside) ) #IMPORTANT: - sign because of normal

        #Riemann Problem Right
        #inside    L
        #outside   R
        if inde==N_el-1: #Care for the BC
            #Riemann Problem Right
            H_outside, B_outside, v_outside = test_dependent.BC_state(DATA, x_v[-1], H_field[-1,-1], B_field[-1,-1], v_field[-1] ,H_field[0,0], B_field[0,0], v_field[0],"R")
            H_inside  = H_field[inde,-1]                          #Last from the current cell
            B_inside  = B_field[inde,-1]                          #Last from the current cell
            v_inside  = v_field[global_indices_v[-1]]             #v_continuous from last node of the current cell
            q_outside = np.array([H_outside,H_outside*v_outside]) #L state
            q_inside  = np.array([H_inside ,H_inside *v_inside])  #R state
            if DATA.WB==True:
                # Hhat, HUhat, Bhat=Riemann_solver.shallow_water_hll_WB(q_inside, q_outside, B_inside, B_outside, DATA.g)
                Hhat, HUhat, Bhat=Riemann_solver.shallow_water_hll_WB(q_outside, q_inside, B_outside, B_inside, DATA.g)
            else:
                # Hhat, HUhat = Riemann_solver.shallow_water_hll(q_inside, q_outside, DATA.g)
                Hhat, HUhat = Riemann_solver.shallow_water_hll(q_outside, q_inside, DATA.g)

                Bhat = test_dependent.Bathymetry(x_v[global_indices_v[-1]],DATA)

            CT_phi_v[global_indices_v[-1]]=CT_phi_v[global_indices_v[-1]]+ ( (Hhat+Bhat)-(H_inside+B_inside) ) 
        else:
            #Riemann Problem Right
            H_outside = H_field[inde+1,0]                        #First from the right cell
            H_inside  = H_field[inde,-1]                         #Last from the current cell
            B_outside = B_field[inde+1,0]                        #First from the right cell
            B_inside  = B_field[inde,-1]                         #Last from the current cell
            v_both    = v_field[global_indices_v[-1]]            #v_continuous from last node of the current cell
            q_outside = np.array([H_outside,H_outside*v_both])   #L state
            q_inside  = np.array([H_inside ,H_inside *v_both])   #R state
            if DATA.WB==True:
                # Hhat, HUhat, Bhat=Riemann_solver.shallow_water_hll_WB(q_inside, q_outside, B_inside, B_outside, DATA.g)
                Hhat, HUhat, Bhat=Riemann_solver.shallow_water_hll_WB(q_outside, q_inside, B_outside, B_inside, DATA.g)
            else:
                # Hhat, HUhat = Riemann_solver.shallow_water_hll(q_inside, q_outside, DATA.g)
                Hhat, HUhat = Riemann_solver.shallow_water_hll(q_outside, q_inside, DATA.g)

                Bhat = test_dependent.Bathymetry(x_v[global_indices_v[-1]],DATA)


            #Add contribution
            CT_phi_v[global_indices_v[-1]]=CT_phi_v[global_indices_v[-1]]+ ( (Hhat+Bhat)-(H_inside+B_inside) ) 

    if DATA.periodic==True:
        CT_phi_v[0]=CT_phi_v[0]+CT_phi_v[-1]
        CT_phi_v[-1]=CT_phi_v[0]

    return CT_phi_v
#==============================================================
#
#
#
#==============================================================
# Local computation for the Lax-Friedrichs stabilization
# NB: Uses Lax_Friedrichs_K
#==============================================================
def Lax_Friedrichs_K(v_local,Hmax):
    """
    Computation
    phi^K_i=alpha (c_i-cbar_K)
    """
    #N_local_nodes_v=len(v_local)
    alpha=np.max(np.absolute(v_local)+1e-6)+np.sqrt(9.81*Hmax)
    vbar=np.average(v_local)
    ST_i_K=alpha*(v_local-vbar)
    return ST_i_K
#==============================================================
#
#
#
#==============================================================
# Lax-Friedrichs stabilization
# NB: Uses Lax_Friedrichs_K
#==============================================================
def Lax_Friedrichs(v_field,M_Local_to_Global,H_field,x_v,DATA):
    """
    Computation
    phi_i=sum_{K in K_i} alpha (c_i-cbar_K)
    """

    N_global_nodes_v=len(v_field)
    N_el, N_local_nodes_v = M_Local_to_Global.shape

    ST_i=np.zeros(N_global_nodes_v)

    for inde in range(N_el):

        Hmax=max(H_field[inde,:])
        Hmax=0
        global_indices_v=M_Local_to_Global[inde,:]
        v_local=v_field[global_indices_v]
        ST_i_K=Lax_Friedrichs_K(v_local,Hmax)
        for indi in range(N_local_nodes_v):
            ST_i[global_indices_v[indi]]=ST_i[global_indices_v[indi]]+ST_i_K[indi]

    if DATA.periodic==True:
        ST_i[0]=ST_i[0]+ST_i[-1]
        ST_i[-1]=ST_i[0]


    return ST_i
#==============================================================
#
#
#
#==============================================================
# Computation of the time step
#==============================================================
def Compute_Time_Step(H_field,v_field,x_v,M_Local_to_Global,DATA,degree_v,CFL):
    """
    Computation of the time step
    NB: I divide by (2*degree_v+1) to compute the CFL
    """
    N_local_nodes_H=H_field.shape[1]
    N_el, N_local_nodes_v = M_Local_to_Global.shape
    N_global_nodes_v=(N_local_nodes_v-1)*N_el+1

    dt_max=DATA.T
    for inde in range(N_el):
        H_local=H_field[inde,:]
        global_indices_v=M_Local_to_Global[inde,:]
        v_local=v_field[global_indices_v]
        x_local=x_v[global_indices_v]

        max_sound_speed=0
        max_H=0
        max_abs_v=0
        for indi in range(N_local_nodes_H):
            max_H=max(max_H,H_local[indi])
        for indi in range(N_local_nodes_v):
            max_abs_v=max(max_abs_v,abs(v_local[indi]))

        max_sound_speed=np.sqrt(DATA.g*max_H)
        dx=x_local[-1]-x_local[0]

        dt_max=min(dt_max,dx/(max_sound_speed+max_abs_v))

    dt_max=dt_max/(2*degree_v+1)*CFL

    return dt_max
#==============================================================
#
#
#
#==============================================================
# Getting location of the x_H DoF from the x_v field
#==============================================================
def get_x_H(x_v,local_values_v_in_H,M_Local_to_Global):
    """
    Getting location of the x_H DoF from the x_v field
    """
    N_global_nodes_v=len(x_v)
    N_local_nodes_v, N_local_nodes_H = local_values_v_in_H.shape
    N_el=int((N_global_nodes_v-1)/(N_local_nodes_v-1))

    x_H=np.zeros((N_el,N_local_nodes_H))

    in_H_local_values_v=local_values_v_in_H.transpose()

    for inde in range(N_el):    

        global_indices_v=M_Local_to_Global[inde,:]
        x_local=x_v[global_indices_v]

        for indi in range(N_local_nodes_H):
            x_H[inde,indi]=sum(x_local*in_H_local_values_v[indi,:])

    return x_H
#==============================================================
#
#
#
#==============================================================
# Getting bathymetry from x_H 
#==============================================================
def get_B(x_H,DATA):
    """
    Getting bathymetry from x_H 
    """

    N_el, N_local_nodes_H=x_H.shape
    B_field=np.zeros((N_el,N_local_nodes_H))

    for inde in range(N_el):    
        for indi in range(N_local_nodes_H):
            B_field[inde,indi]=test_dependent.Bathymetry(x_H[inde,indi],DATA)

    return B_field
#==============================================================
#
#
#
#==============================================================
# Getting v at x_H 
#==============================================================
def get_v_in_x_H(v_field,x_H,local_values_v_in_H,M_Local_to_Global):
    """
    Getting v on x_H 
    """

    N_el, N_local_nodes_H=x_H.shape
    v_in_x_H=np.zeros((N_el,N_local_nodes_H))
    in_H_local_values_v=local_values_v_in_H.transpose()

    for inde in range(N_el): 
        global_indices_v = M_Local_to_Global[inde,:]    
        v_local          = v_field[global_indices_v]
        v_in_x_H[inde,:]   = in_H_local_values_v @ v_local

    return v_in_x_H
#==============================================================
#
#
#
#==============================================================
# Getting H at x_v
#==============================================================
def get_H_in_x_v(H_field,x_v,local_values_H_in_v,M_Local_to_Global):
    """
    Getting H on x_v 
    """

    N_global_nodes_v=len(x_v)
    H_in_x_v=np.zeros((N_global_nodes_v))
    in_v_local_values_H=local_values_H_in_v.transpose()
    N_el, N_local_nodes_H=H_field.shape
    N_local_nodes_v=N_local_nodes_H+1

    for inde in range(N_el): 
        H_local          = H_field[inde,:]
        global_indices_v = M_Local_to_Global[inde,:]    
        H_in_x_v[global_indices_v] = H_in_x_v[global_indices_v] + in_v_local_values_H @ H_local

    #NB: The exterma have been counted twice
    H_in_x_v[::(N_local_nodes_v-1)]=H_in_x_v[::(N_local_nodes_v-1)]/2
    H_in_x_v[0]=H_in_x_v[0]*2
    H_in_x_v[-1]=H_in_x_v[-1]*2

    return H_in_x_v


def jump_stabilization(v_field,x_v,local_derivatives_v,M_Local_to_Global,M_faces,H_field,DATA): 

    #In local_derivatives_v we have
    #Rows basis functions
    #Columns x_j_v
    #I transpose it to have in each row a specifc x_i_v and in the columns the basis functions
    in_v_local_derivatives_v=local_derivatives_v.transpose()

    N_el, N_local_nodes_v = M_Local_to_Global.shape
    N_f=len(M_faces)

    phi_jump=np.zeros(len(v_field))

    for indf in range(N_f):
        el_L=M_faces[indf][0]
        el_R=M_faces[indf][1]
        if el_L==-1 or el_R==-1:
            #Not in the interior, skip
            continue
        else:

            #Computation of the jump
            global_indices_L = M_Local_to_Global[el_L,:]
            global_indices_R = M_Local_to_Global[el_R,:]

            v_local_L        = v_field[global_indices_L]
            v_local_R        = v_field[global_indices_R]

            x_v_local_L      = x_v[global_indices_L]
            x_v_local_R      = x_v[global_indices_R]

            dxL              = np.sum(x_v_local_L*in_v_local_derivatives_v[-1,:])
            dxR              = np.sum(x_v_local_R*in_v_local_derivatives_v[0,:])

            dxiL             = 1/dxL
            dxiR             = 1/dxR

            dv_L             = np.sum(v_local_L*in_v_local_derivatives_v[-1,:]) * dxiL
            dv_R             = np.sum(v_local_R*in_v_local_derivatives_v[0,:])  * dxiR

            jump_dv=dv_R-dv_L

            #Computation of the spectral radius
            H=(H_field[el_L,-1]+H_field[el_R,0])/2
            sr=np.abs(v_local_R[0])+np.sqrt(DATA.g*H)

            #Computation of dx
            deltaxL=x_v_local_L[-1]-x_v_local_L[0]
            deltaxR=x_v_local_R[-1]-x_v_local_R[0]
            deltax=(deltaxR+deltaxL)/2


            phi_L=np.zeros(N_local_nodes_v)
            phi_R=np.zeros(N_local_nodes_v)

            #CIP coefficient
            alpha=DATA.delta_CIP*sr*deltax**2

            #Element left
            for indi in range(N_local_nodes_v):
                phi_L[indi]=alpha*jump_dv*(-in_v_local_derivatives_v[-1,indi])*dxiL

            #Element right
            for indi in range(N_local_nodes_v):
                phi_R[indi]=alpha*jump_dv*(in_v_local_derivatives_v[0,indi])*dxiR


            phi_jump[global_indices_L]=phi_jump[global_indices_L]+phi_L
            phi_jump[global_indices_R]=phi_jump[global_indices_R]+phi_R


    if DATA.periodic==True:
        phi_jump[0]=phi_jump[0]+phi_jump[-1]
        phi_jump[-1]=phi_jump[0]

    return phi_jump

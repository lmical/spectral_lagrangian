import numpy as np

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
            #print(inde,indi,sum(in_H_local_derivatives_v[indi,:]*x_v_local))
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
            #print(inde,indi,sum(in_H_local_derivatives_v[indi,:]*x_v_local))
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
def Lumped_Mass_Matrix(w_v,x_v,M_Local_to_Global,local_derivatives_v):
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
    
    return M_v

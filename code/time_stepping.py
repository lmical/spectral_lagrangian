import numpy as np
import lagrangian_scheme
import test_dependent
import DeC



#==============================================================
# Euler method
#==============================================================
def Euler_method(dt,H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA):
    #Update
    #x
    x_v=x_v_old+dt*v_field_old
    #v
    v_field=v_field_old+dt*rhs_v_function(H_field_old,v_field_old,x_v_old,B_field_old, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    #H, with strong mass conservation
    H_field=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v,local_derivatives_v_in_H,M_Local_to_Global)
    #B
    x_H=lagrangian_scheme.get_x_H(x_v,local_values_v_in_H,M_Local_to_Global)
    B_field=lagrangian_scheme.get_B(x_H,DATA)

    #---------------------------------------------------------------------------------
    # Test, explicit computation of the rhs
    # M_v=lagrangian_scheme.Lumped_Mass_Matrix(w_v,x_v_old,M_Local_to_Global,local_derivatives_v)
    # phi_v=lagrangian_scheme.Space_Residuals_v(H_field_old, B_field_old, w_v,local_derivatives_H_in_v,M_Local_to_Global)
    # CT_phi_v=lagrangian_scheme.Coupling_Terms_Space_Residuals_v(H_field_old, B_field_old, v_field_old, M_Local_to_Global, M_faces, DATA)
    # ST_i=lagrangian_scheme.Lax_Friedrichs(v_field_old,M_Local_to_Global,H_field_old,x_v_old)
    # v_field=v_field_old+dt*(-DATA.g*phi_v-DATA.g*CT_phi_v-ST_i)/M_v
    # if np.linalg.norm( ((-DATA.g*phi_v-DATA.g*CT_phi_v-ST_i)/M_v)- rhs_v_function(H_field_old,v_field_old,x_v_old,B_field_old, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA) )>1e-14:
    #     print("Problem")
    #     quit()
    #---------------------------------------------------------------------------------

    return H_field, v_field, x_v, B_field
#==============================================================
#
#
#
#==============================================================
# rhs function for the DeC update of v
#==============================================================
def rhs_v_function(H_field,v_field,x_v,B_field, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA):
    M_v=lagrangian_scheme.Lumped_Mass_Matrix(w_v,x_v,M_Local_to_Global,local_derivatives_v)
    phi_v=lagrangian_scheme.Space_Residuals_v(H_field, B_field, w_v,local_derivatives_H_in_v,M_Local_to_Global)
    CT_phi_v=lagrangian_scheme.Coupling_Terms_Space_Residuals_v(H_field, B_field, v_field, M_Local_to_Global, M_faces, x_v, DATA)
    if DATA.LaxFriedrichs==True:
        ST_i=lagrangian_scheme.Lax_Friedrichs(v_field,M_Local_to_Global,H_field,x_v)
    else:
        ST_i=np.zeros(len(phi_v))

    if DATA.periodic==True:
        M_v[0]=M_v[0]+M_v[-1]
        M_v[-1]=M_v[0]

        phi_v[0]=phi_v[0]+phi_v[-1]
        phi_v[-1]=phi_v[0]

        CT_phi_v[0]=CT_phi_v[0]+CT_phi_v[-1]
        CT_phi_v[-1]=CT_phi_v[0]

        ST_i[0]=ST_i[0]+ST_i[-1]
        ST_i[-1]=ST_i[0]


    return (-DATA.g*phi_v-DATA.g*CT_phi_v-ST_i)/M_v
#==============================================================
#
#
#
#==============================================================
# DeC method
#==============================================================
def DeC_method(dt,H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA,dec):

    #Initialization of the structures
    #p = previous iteration
    #a = actual iteration 

    N_el, local_nodes_H = H_field_old.shape
    N_global_nodes_v    = len(v_field_old)
    N_el, N_local_nodes_v = M_Local_to_Global.shape

    H_p     = np.zeros((N_el, local_nodes_H, dec.M_sub+1))
    H_a     = np.zeros((N_el, local_nodes_H, dec.M_sub+1))

    B_p     = np.zeros((N_el, local_nodes_H, dec.M_sub+1))
    B_a     = np.zeros((N_el, local_nodes_H, dec.M_sub+1))


    v_p     = np.zeros((N_global_nodes_v, dec.M_sub+1))
    v_a     = np.zeros((N_global_nodes_v, dec.M_sub+1))

    x_v_p   = np.zeros((N_global_nodes_v, dec.M_sub+1))
    x_v_a   = np.zeros((N_global_nodes_v, dec.M_sub+1))

    rhs_v   = np.zeros((N_global_nodes_v,dec.M_sub+1))
    rhs_x_v = np.zeros((N_global_nodes_v,dec.M_sub+1))


    #Filling p structures with u^{(0)}
    for inds in range(dec.M_sub+1):
        H_p[:,:,inds] = H_field_old[:,:].copy()
        B_p[:,:,inds] = B_field_old[:,:].copy()
        v_p[:,inds]   = v_field_old[:].copy()
        x_v_p[:,inds] = x_v_old[:].copy()


    #Filling a structure in subtimenode 0 (just for the sake of completeness)
    H_a[:,:,0] = H_field_old[:,:].copy()
    B_a[:,:,0] = B_field_old[:,:].copy()
    v_a[:,0]   = v_field_old[:].copy()
    x_v_a[:,0] = x_v_old[:].copy()


    #DeC iteration loop
    rhs_v[:,0]   = rhs_v_function(H_field_old,v_field_old,x_v_old,B_field_old, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    rhs_x_v[:,0] = v_field_old.copy()

    for r in range(1,dec.M_sub+1):
        rhs_v[:,r]   = rhs_v[:,0].copy()
        rhs_x_v[:,r] = rhs_x_v[:,0].copy()
    for k in range(1,dec.n_iter+1):
        H_p=np.copy(H_a)
        B_p=np.copy(B_a)
        v_p=np.copy(v_a)
        x_v_p=np.copy(x_v_a)


        if k>1:
            for r in range(1,dec.M_sub+1):
                rhs_v[:,r]   = rhs_v_function(H_p[:,:,r],v_p[:,r],x_v_p[:,r],B_p[:,:,r], w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
                rhs_x_v[:,r] = v_p[:,r]
        if k < dec.n_iter:
            for m in range(1,dec.M_sub+1):
                v_a[:,m] = v_field_old[:] + dt*sum([dec.theta[r,m]*rhs_v[:,r] for r in range(dec.M_sub+1)])
                x_v_a[:,m] = x_v_old[:]   + dt*sum([dec.theta[r,m]*rhs_x_v[:,r] for r in range(dec.M_sub+1)])
                H_a[:,:,m]=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_a[:,m],local_derivatives_v_in_H,M_Local_to_Global)
                x_H=lagrangian_scheme.get_x_H(x_v_a[:,m],local_values_v_in_H,M_Local_to_Global)
                B_a[:,:,m]=lagrangian_scheme.get_B(x_H,DATA)
        else:
            v_a[:,dec.M_sub]= v_field_old[:]+dt*sum([dec.theta[r,dec.M_sub]*rhs_v[:,r] for r in range(dec.M_sub+1)])
            x_v_a[:,dec.M_sub]= x_v_old[:]+dt*sum([dec.theta[r,dec.M_sub]*rhs_x_v[:,r] for r in range(dec.M_sub+1)])
            H_a[:,:,dec.M_sub]=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_a[:,dec.M_sub],local_derivatives_v_in_H,M_Local_to_Global)
            x_H=lagrangian_scheme.get_x_H(x_v_a[:,dec.M_sub],local_values_v_in_H,M_Local_to_Global)
            B_a[:,:,dec.M_sub]=lagrangian_scheme.get_B(x_H,DATA)

    H_field = H_a[:,:,dec.M_sub]
    v_field = v_a[:,dec.M_sub]
    x_v     = x_v_a[:,dec.M_sub]
    B_field = B_a[:,:,dec.M_sub]

    return H_field, v_field, x_v, B_field




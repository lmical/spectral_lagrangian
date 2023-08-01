import numpy as np
import lagrangian_scheme
import test_dependent
import DeC



#==============================================================
# Euler method
#==============================================================
def Euler_method(H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA):
    #Update
    #x
    x_v=x_v_old+DATA.dt*v_field_old
    #v
    v_field=v_field_old+DATA.dt*rhs_v_function(H_field_old,v_field_old,x_v_old,B_field_old, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
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
    M_v=lagrangian_scheme.Lumped_Mass_Matrix(w_v,x_v,M_Local_to_Global,local_derivatives_v,DATA)
    phi_v=lagrangian_scheme.Space_Residuals_v(H_field, B_field, w_v,local_derivatives_H_in_v,M_Local_to_Global,DATA)
    CT_phi_v=lagrangian_scheme.Coupling_Terms_Space_Residuals_v(H_field, B_field, v_field, M_Local_to_Global, M_faces, x_v, DATA)

    if DATA.LaxFriedrichs==True:
        ST_i=lagrangian_scheme.Lax_Friedrichs(v_field,M_Local_to_Global,H_field,x_v,DATA)
    else:
        ST_i=np.zeros(len(phi_v))

    if DATA.jump=="jc":
        phi_jump=lagrangian_scheme.jump_stabilization(v_field,x_v,local_derivatives_v,M_Local_to_Global,M_faces,H_field,DATA)
    else:
        phi_jump=np.zeros(len(phi_v))



    return (-DATA.g*phi_v-DATA.g*CT_phi_v-ST_i-phi_jump)/M_v
#==============================================================
#
#
#
#==============================================================
# DeC method
#==============================================================
def DeC_method(H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA,dec):

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
                rhs_x_v[:,r] = v_p[:,r].copy()
        if k < dec.n_iter:
            for m in range(1,dec.M_sub+1):
                v_a[:,m] = v_field_old[:] + DATA.dt*sum([dec.theta[r,m]*rhs_v[:,r] for r in range(dec.M_sub+1)])
                x_v_a[:,m] = x_v_old[:]   + DATA.dt*sum([dec.theta[r,m]*rhs_x_v[:,r] for r in range(dec.M_sub+1)])
                H_a[:,:,m]=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_a[:,m],local_derivatives_v_in_H,M_Local_to_Global)
                x_H=lagrangian_scheme.get_x_H(x_v_a[:,m],local_values_v_in_H,M_Local_to_Global)
                B_a[:,:,m]=lagrangian_scheme.get_B(x_H,DATA)
        else:
            v_a[:,dec.M_sub]= v_field_old[:]+DATA.dt*sum([dec.theta[r,dec.M_sub]*rhs_v[:,r] for r in range(dec.M_sub+1)])
            x_v_a[:,dec.M_sub]= x_v_old[:]+DATA.dt*sum([dec.theta[r,dec.M_sub]*rhs_x_v[:,r] for r in range(dec.M_sub+1)])
            H_a[:,:,dec.M_sub]=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_a[:,dec.M_sub],local_derivatives_v_in_H,M_Local_to_Global)
            x_H=lagrangian_scheme.get_x_H(x_v_a[:,dec.M_sub],local_values_v_in_H,M_Local_to_Global)
            B_a[:,:,dec.M_sub]=lagrangian_scheme.get_B(x_H,DATA)

    H_field = H_a[:,:,dec.M_sub]
    v_field = v_a[:,dec.M_sub]
    x_v     = x_v_a[:,dec.M_sub]
    B_field = B_a[:,:,dec.M_sub]

    return H_field, v_field, x_v, B_field
#==============================================================
#
#
#
#==============================================================
# SSPRK4 method
#==============================================================
def SSPRK4_method(H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA):

    #Initialization of the structures

    N_el, local_nodes_H = H_field_old.shape
    N_global_nodes_v    = len(v_field_old)
    N_el, N_local_nodes_v = M_Local_to_Global.shape


    #K0(1:nVar,1:nElemsX,1:nElemsY) = U(1:nVar,1:nElemsX,1:nElemsY)

    H_0   = H_field_old.copy()
    B_0   = B_field_old.copy()
    v_0   = v_field_old.copy()
    x_v_0 = x_v_old.copy()

    #--------------------!
    # First Stage        !
    #--------------------!
    tStage = DATA.time + 0.0*DATA.dt

    # CALL FVTimeDerivative(tStage)

    rhs_v   = rhs_v_function(H_0,v_0,x_v_0,B_0, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    rhs_x_v = v_0.copy()

    # K1(1:nVar,1:nElemsX,1:nElemsY) = &
    #     1.00000000000000*K0(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.39175222700392*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

    v_1   = v_0   + DATA.dt*0.39175222700392*rhs_v
    x_v_1 = x_v_0 + DATA.dt*0.39175222700392*rhs_x_v
    H_1   = lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_1,local_derivatives_v_in_H,M_Local_to_Global)
    x_H   = lagrangian_scheme.get_x_H(x_v_1,local_values_v_in_H,M_Local_to_Global)
    B_1   = lagrangian_scheme.get_B(x_H,DATA)


    # U(1:nVar,1:nElemsX,1:nElemsY)  = K1(1:nVar,1:nElemsX,1:nElemsY)

    #--------------------!
    # Second Stage       !
    #--------------------!
    tStage = DATA.time + 0.39175222700392*DATA.dt

    # CALL FVTimeDerivative(tStage)

    rhs_v   = rhs_v_function(H_1,v_1,x_v_1,B_1, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    rhs_x_v = v_1.copy()

    # K2(1:nVar,1:nElemsX,1:nElemsY) = &
    #     0.44437049406734*K0(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.55562950593266*K1(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.36841059262959*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt


    v_2   = 0.44437049406734*v_0   + 0.55562950593266*v_1   + 0.36841059262959*DATA.dt*rhs_v
    x_v_2 = 0.44437049406734*x_v_0 + 0.55562950593266*x_v_1 + 0.36841059262959*DATA.dt*rhs_x_v
    H_2   = lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_2,local_derivatives_v_in_H,M_Local_to_Global)
    x_H   = lagrangian_scheme.get_x_H(x_v_2,local_values_v_in_H,M_Local_to_Global)
    B_2   = lagrangian_scheme.get_B(x_H,DATA)

    # U(1:nVar,1:nElemsX,1:nElemsY)  = K2(1:nVar,1:nElemsX,1:nElemsY)

    #--------------------!
    # Third Stage        !
    #--------------------!
    tStage = DATA.time + 0.58607968896780*DATA.dt

    # CALL FVTimeDerivative(tStage)

    rhs_v   = rhs_v_function(H_2,v_2,x_v_2,B_2, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    rhs_x_v = v_2.copy()

    # K3(1:nVar,1:nElemsX,1:nElemsY) = &
    #     0.62010185138540*K0(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.37989814861460*K2(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.25189177424738*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

    v_3   = 0.62010185138540*v_0   + 0.37989814861460*v_2   + 0.25189177424738*DATA.dt*rhs_v
    x_v_3 = 0.62010185138540*x_v_0 + 0.37989814861460*x_v_2 + 0.25189177424738*DATA.dt*rhs_x_v
    H_3   = lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_3,local_derivatives_v_in_H,M_Local_to_Global)
    x_H   = lagrangian_scheme.get_x_H(x_v_3,local_values_v_in_H,M_Local_to_Global)
    B_3   = lagrangian_scheme.get_B(x_H,DATA)


    # U(1:nVar,1:nElemsX,1:nElemsY)  = K3(1:nVar,1:nElemsX,1:nElemsY)

    #--------------------!
    # Fourth Stage       !
    #--------------------!
    tStage = DATA.time + 0.474542364687*DATA.dt

    # CALL FVTimeDerivative(tStage)

    rhs_v   = rhs_v_function(H_3,v_3,x_v_3,B_3, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    rhs_x_v = v_3.copy()

    # K4(1:nVar,1:nElemsX,1:nElemsY) = &
    #     0.17807995410773*K0(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.82192004589227*K3(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.54497475021237*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

    v_4   = 0.17807995410773*v_0   + 0.82192004589227*v_3   + 0.54497475021237*DATA.dt*rhs_v
    x_v_4 = 0.17807995410773*x_v_0 + 0.82192004589227*x_v_3 + 0.54497475021237*DATA.dt*rhs_x_v
    H_4   = lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v_4,local_derivatives_v_in_H,M_Local_to_Global)
    x_H   = lagrangian_scheme.get_x_H(x_v_4,local_values_v_in_H,M_Local_to_Global)
    B_4   = lagrangian_scheme.get_B(x_H,DATA)


    # U(1:nVar,1:nElemsX,1:nElemsY)  = K4(1:nVar,1:nElemsX,1:nElemsY)
    # K5(1:nVar,1:nElemsX,1:nElemsY) = Ut(1:nVar,1:nElemsX,1:nElemsY)

    rhs_v_old=rhs_v.copy()
    rhs_x_v_old    =rhs_x_v.copy()

    #--------------------!
    # Fifth Stage        !
    #--------------------!
    tStage  = DATA.time + 0.93501063100924*DATA.dt

    # CALL FVTimeDerivative(tStage)

    rhs_v   = rhs_v_function(H_4,v_4,x_v_4,B_4, w_v, local_derivatives_v, local_derivatives_H_in_v, M_Local_to_Global, M_faces, DATA)
    rhs_x_v = v_4.copy()


    # U(1:nVar,1:nElemsX,1:nElemsY)  = &
    #     0.00683325884039*K0(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.51723167208978*K2(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.12759831133288*K3(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.34833675773694*K4(1:nVar,1:nElemsX,1:nElemsY) &
    # + 0.08460416338212*K5(1:nVar,1:nElemsX,1:nElemsY)*dt &
    # + 0.22600748319395*Ut(1:nVar,1:nElemsX,1:nElemsY)*dt

    v_field = 0.00683325884039*v_0   + 0.51723167208978*v_2   + 0.12759831133288*v_3   + 0.34833675773694*v_4   + 0.08460416338212*DATA.dt*rhs_v_old   + 0.22600748319395*DATA.dt*rhs_v   
    x_v     = 0.00683325884039*x_v_0 + 0.51723167208978*x_v_2 + 0.12759831133288*x_v_3 + 0.34833675773694*x_v_4 + 0.08460416338212*DATA.dt*rhs_x_v_old + 0.22600748319395*DATA.dt*rhs_x_v 
    H_field = lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v,local_derivatives_v_in_H,M_Local_to_Global)
    x_H     = lagrangian_scheme.get_x_H(x_v,local_values_v_in_H,M_Local_to_Global)
    B_field = lagrangian_scheme.get_B(x_H,DATA)


    return H_field, v_field, x_v, B_field







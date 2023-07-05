import numpy as np

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

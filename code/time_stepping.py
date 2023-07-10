import numpy as np
import lagrangian_scheme
import test_dependent




#==============================================================
# Euler method
#==============================================================
def Euler(dt,H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, M_faces, DATA):
    #Update
    #x
    x_v=x_v_old+dt*v_field_old
    #v
    M_v=lagrangian_scheme.Lumped_Mass_Matrix(w_v,x_v_old,M_Local_to_Global,local_derivatives_v)
    phi_v=lagrangian_scheme.Space_Residuals_v(H_field_old, B_field_old, w_v,local_derivatives_H_in_v,M_Local_to_Global)
    CT_phi_v=lagrangian_scheme.Coupling_Terms_Space_Residuals_v(H_field_old, B_field_old, v_field_old, M_Local_to_Global, M_faces, DATA)
    ST_i=lagrangian_scheme.Lax_Friedrichs(v_field_old,M_Local_to_Global,H_field_old,x_v_old)
    v_field=v_field_old+dt*(-DATA.g*phi_v-DATA.g*CT_phi_v-ST_i)/M_v
    #H, with strong mass conservation
    H_field=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v,local_derivatives_v_in_H,M_Local_to_Global)
    #BC
    H_field,v_field,x_v=test_dependent.BC(H_field,v_field,B_field_old,x_v,DATA)

    return H_field, v_field, x_v
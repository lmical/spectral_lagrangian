import numpy as np

#==============================================================
# Riemann solver
#--------------------------------------------------------------
# HLL from: http://www.clawpack.org/riemann_book/html/Shallow_water_approximate.html
# NB: Riemann solver in conserved variables
#==============================================================
def shallow_water_hll(q_l, q_r, xi, grav):
    """
    HLLE approximate solver for the shallow water equations.
    """
    h_l = q_l[0]
    hu_l = q_l[1]
    u_l = hu_l/h_l
    h_r = q_r[0]
    hu_r = q_r[1]
    u_r = hu_r/h_r
        
    # Roe averages
    h_hat = (h_r + h_l)/2.
    u_hat = (np.sqrt(h_r)*u_r + np.sqrt(h_l)*u_l) / \
            (np.sqrt(h_r) + np.sqrt(h_l))
    c_hat = np.sqrt(grav*h_hat)

    lambda_1_l = u_l - np.sqrt(grav*h_l)
    lambda_2_r = u_r + np.sqrt(grav*h_r)
    
    s1 = min(lambda_1_l, u_hat - c_hat)
    s2 = max(lambda_2_r, u_hat + c_hat)
    
    h_m = (hu_r - hu_l - s2*h_r + s1*h_l)/(s1-s2)
    hu_m = (hu_r*u_r - hu_l*u_l + 0.5*grav*(h_r**2 - h_l**2) \
            - s2*hu_r + s1*hu_l)/(s1-s2)
    

    h_out  = (xi<s1)*h_l + (s1<=xi)*(xi<=s2)*h_m + (s2<xi)*h_r
    hu_out = (xi<s1)*hu_l + (s1<=xi)*(xi<=s2)*hu_m + (s2<xi)*hu_r
    return h_out, hu_out



#==============================================================
# WB Riemann solver
#--------------------------------------------------------------
# Based on HLL from: http://www.clawpack.org/riemann_book/html/Shallow_water_approximate.html
# Modification from: A simple well-balanced and positive numerical scheme for the shallow-water system 
# NB: Riemann solver in conserved variables
#==============================================================
def shallow_water_hll_WB(q_l, q_r, b_l, b_r, xi, grav):
    """
    HLLE approximate solver for the shallow water equations.
    With a WB modification
    """
    h_l = q_l[0]
    hu_l = q_l[1]
    u_l = hu_l/h_l
    h_r = q_r[0]
    hu_r = q_r[1]
    u_r = hu_r/h_r
        
    # Roe averages
    h_hat = (h_r + h_l)/2.
    u_hat = (np.sqrt(h_r)*u_r + np.sqrt(h_l)*u_l) / \
            (np.sqrt(h_r) + np.sqrt(h_l))
    c_hat = np.sqrt(grav*h_hat)

    lambda_1_l = u_l - np.sqrt(grav*h_l)
    lambda_2_r = u_r + np.sqrt(grav*h_r)
    
    s1 = min(lambda_1_l, u_hat - c_hat)
    s2 = max(lambda_2_r, u_hat + c_hat)
    
    h_m = (hu_r - hu_l - s2*h_r + s1*h_l)/(s1-s2)
    hu_m = (hu_r*u_r - hu_l*u_l + 0.5*grav*(h_r**2 - h_l**2) \
            - s2*hu_r + s1*hu_l)/(s1-s2)
    

    if xi<s1:
        h_out=h_l
        hu_out=hu_l
        b_out=b_l
    elif (s1<=xi) and (xi<=s2):
        h_out=h_m
        hu_out=hu_m
        #BATHYMETRY AT THE INTERFACE GUARANTEEING WB
        b_out=b_r*(u_hat-s2)/(s1-s2)+b_l*(s1-u_hat)/(s1-s2) 
    else:
        h_out=h_r
        hu_out=hu_r
        b_out=b_r


    return h_out, hu_out, b_out

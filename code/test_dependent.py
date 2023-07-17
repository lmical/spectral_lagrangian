import numpy as np

#==============================================================
# DATA of the simulation
# xL,       left boundary
# xR,       right boundary
# T,        final time
# periodic, periodicity of the mesh
# gravity
#==============================================================
class DATA_CLASS:
    def __init__(self,test,N_el,order_space,time_scheme,order_time,CFL,freq,N_max_iter,scheme,LaxFriedrichs,jump,folder,printing,plotting):

        #Somhow it is better to store also the input paramters in DATA
        self.test               = test
        self.N_el               = N_el
        self.order_space        = order_space
        self.time_scheme        = time_scheme
        self.order_time         = order_time
        self.CFL                = CFL           
        self.freq               = freq              
        self.N_max_iter         = N_max_iter             
        self.scheme             = scheme
        self.LaxFriedrichs      = LaxFriedrichs
        self.jump               = jump
        self.folder             = folder
        self.printing           = printing
        self.plotting           = plotting

        if test=="Sod" or test=="Sod_smooth": #Sod
            # Extrema
            self.xL=0
            self.xR=1
            # Final time
            self.T=0.1
            # Periodicity of the mesh
            self.periodic=False
            # gravity
            self.g=9.81
        elif test=="Smooth_periodic": #Sod
            # Extrema
            self.xL=0
            self.xR=1
            # Final time
            self.T=0.05
            # Periodicity of the mesh
            self.periodic=True
            # gravity
            self.g=9.81
        else:
            print("Error in class DATA_CLASS, in test_dependent, test not available")
            quit()
#==============================================================
#
#
#
#==============================================================
# Some state, e.g., an initial condition in a specific point x at a specific time t
# Used to define the IC
#==============================================================
def Analytical_State(x,t,DATA):
    if DATA.test=="Sod":
        v=0.
        if x<=0.5:
            H=1.
        else:
            H=0.2
    elif DATA.test=="Sod_smooth":
        v=0.
        r0=0.45
        r1=0.55
        if x<=r0:
            H=1.
        elif x>r0 and x<r1:
            H=0.2+0.8*np.exp(1. - 1./(1.-((x-r0)/(r1-r0))**2))    
        else:
            H=0.2
    elif DATA.test=="Smooth_periodic":
        H=2+np.cos(2*np.pi*x)
        v=1.
    else:
        print("Error in function Analytical_State, test not available")
        quit()
    
    return H,v
#==============================================================
#
#
#
#==============================================================
# Bathymetry in a point x
#==============================================================
def Bathymetry(x,DATA):
    if DATA.test=="Sod" or DATA.test=="Sod_smooth" or DATA.test=="Smooth_periodic":
        B=0.
    else:
        print("Error in function Bathymetry, test not available")
        quit()
    return B
#==============================================================
#
#
#
#==============================================================
# Derivative of the bathymetry in a point x
#==============================================================
def Derivative_Bathymetry(x,DATA):
    if DATA.test=="Sod" or DATA.test=="Sod_smooth" or DATA.test=="Smooth_periodic":
        dB=0.
    else:
        print("Error in function Derivative_Bathymetry, test not available")
        quit()
    return dB
#==============================================================
#
#
#
#==============================================================
# Initial condition
#==============================================================
def IC(x_H,x_v,t,DATA):
    """
    Fill the fields H, B and v
    """
    # Matrix H_field[inde,loc_indi_H], rows=elements, columns=loc_indi_H
    H_field = np.zeros(x_H.shape)
    # Matrix B_field[inde,loc_indi_H], with the same convention
    B_field = np.zeros(x_H.shape)
    # v_field[glob_indi_v]
    v_field = np.zeros((len(x_v)))

    N_el, local_nodes_H = x_H.shape
    if DATA.test=="Sod" or DATA.test=="Sod_smooth":
        for inde in range(N_el):
            for indi_l in range(local_nodes_H):
                vec=Analytical_State(x_H[inde,indi_l],0,DATA)
                H_field[inde,indi_l] = vec[0]
                B_field[inde,indi_l] = Bathymetry(x_H[inde,indi_l],DATA)
        v_field[:]=0.
    elif DATA.test=="Smooth_periodic":
        for inde in range(N_el):
            for indi_l in range(local_nodes_H):
                vec=Analytical_State(x_H[inde,indi_l],0,DATA)
                H_field[inde,indi_l] = vec[0]
                B_field[inde,indi_l] = Bathymetry(x_H[inde,indi_l],DATA)
        v_field[:]=1.
    else:
        print("Error in function IC, test not available")
        quit()
    return H_field, B_field, v_field
#==============================================================
#
#
#
#==============================================================
# Boundary condition
#==============================================================
def BC_state(DATA,x,H_other_side, B_other_side, v_other_side,boundary):
    if DATA.periodic==True:
        H=H_other_side
        B=B_other_side
        v=v_other_side
    elif DATA.test=="Sod" or DATA.test=="Sod_smooth":
        if boundary=="L":
            H=1.
            B=0.
            v=0.
        else:
            H=0.2
            B=0.
            v=0.
    else:
        print("Problems in BC_state, test not available")
        quit()
    return H,B,v
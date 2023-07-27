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
    def __init__(self,test,perturbation,N_el,order_space,time_scheme,order_time,CFL,freq,N_max_iter,scheme,LaxFriedrichs,WB,jump,folder,printing,plotting,storing):

        #Somhow it is better to store also the input paramters in DATA
        self.test               = test
        self.perturbation       = perturbation
        self.N_el               = N_el
        self.order_space        = order_space
        self.time_scheme        = time_scheme
        self.order_time         = order_time
        self.CFL                = CFL           
        self.freq               = freq              
        self.N_max_iter         = N_max_iter             
        self.scheme             = scheme
        self.LaxFriedrichs      = LaxFriedrichs
        self.WB                 = WB
        self.jump               = jump
        self.folder             = folder
        self.printing           = printing
        self.plotting           = plotting
        self.storing            = storing

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
            # Analytical solution
            self.analytical_solution=False
        elif test=="Smooth_periodic": #Smooth for convergence analysis with experimental order of accuracy
            # Extrema
            self.xL=0
            self.xR=1
            # Final time
            self.T=0.07
            # Periodicity of the mesh
            self.periodic=True
            # gravity
            self.g=9.81
            # Analytical solution
            self.analytical_solution=False
        elif test=="Lake_At_Rest_Smooth" or test=="Lake_At_Rest_Not_Smooth": #Lake at rest with smooth or non-smooth bathymetry
            # Extrema
            self.xL=0
            self.xR=25
            # Final time
            self.T=10.
            # Periodicity of the mesh
            self.periodic=True
            # gravity
            self.g=9.81
            # Analytical solution
            self.analytical_solution=True
        elif test=="Supercritical_Smooth": #Supercritical with smooth or non-smooth bathymetry
            # Extrema
            self.xL=0
            self.xR=25
            # Final time
            self.T=0.3
            # Periodicity of the mesh
            self.periodic=False
            # gravity
            self.g=9.81
            # Analytical solution
            self.analytical_solution=True
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
    elif DATA.test=="Lake_At_Rest_Smooth" or DATA.test=="Lake_At_Rest_Not_Smooth":
        H=0.5-Bathymetry(x,DATA)
        v=0.
    elif DATA.test=="Supercritical_Smooth":
        q0=24.
        hL=2.
        #Exact
        b=Bathymetry(x,DATA)
        p=[1., (b-q0**2/(2.*DATA.g*hL**2) - hL), 0., q0**2/(2.*DATA.g)]
        hvec=np.roots(p)
        H=hvec[1]
        v=q0/H
    else:
        print("Error in function Analytical_State in test_dependent module, test not available")
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
    elif DATA.test=="Lake_At_Rest_Smooth" or DATA.test=="Supercritical_Smooth":
        x0=10.
        r=5.
        B=0.
        if (x > x0-r) and (x < x0+r):
            B=0.2*np.exp(1. - 1./(1.-((x-x0)/r)**2.))    
        else:
            B=0.
    elif DATA.test=="Lake_At_Rest_Not_Smooth":
        x0=10.
        r=2.
        B=0.
        if (x > x0-r) and (x < x0+r):
            B=0.2-0.05*(x-x0)**2    
        else:
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
    elif DATA.test=="Lake_At_Rest_Smooth":
        x0=10.
        r=5.
        dB=0.
        if (x > x0-r) and (x < x0+r):
            dB = -0.2*np.exp(1.-1./(1.-((x-x0)/r)**2))* 2.*(x-x0)/r**2/((1.-((x-x0)/r)**2)**2)
        else:
            dB=0.
    elif DATA.test=="Lake_At_Rest_Not_Smooth":
        x0=10.
        r=2.
        dB=0.
        if (x > x0-r) and (x < x0+r):
            dB=-2*0.05*(x-x0)    
        else:
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
    N_global_nodes_v=len(x_v)

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
    elif DATA.test=="Lake_At_Rest_Smooth" or DATA.test=="Lake_At_Rest_Not_Smooth":
        for inde in range(N_el):
            for indi_l in range(local_nodes_H):
                vec=Analytical_State(x_H[inde,indi_l],0,DATA)
                H_field[inde,indi_l] = vec[0]
                B_field[inde,indi_l] = Bathymetry(x_H[inde,indi_l],DATA)
        v_field[:]=0.
    elif DATA.test=="Supercritical_Smooth":
        for inde in range(N_el):
            for indi_l in range(local_nodes_H):
                vec=Analytical_State(x_H[inde,indi_l],0,DATA)
                H_field[inde,indi_l] = vec[0]
                B_field[inde,indi_l] = Bathymetry(x_H[inde,indi_l],DATA)
        for indi_g in range(N_global_nodes_v):
            vec=Analytical_State(x_v[indi_g],0,DATA)
            v_field[indi_g]=vec[1]
    else:
        print("Error in test_dependent module, in function function IC, test not available")
        quit()
    return H_field, B_field, v_field
#==============================================================
#
#
#
#==============================================================
# Smooth perturbation
#==============================================================
def smoot_perturbation(x,x0,r,A):
    p=0.
    if(x>x0-r and x<x0+r):
        p=A*np.exp(1-1/(1-((x-x0)/r)**2))

    return p
#==============================================================
#
#
#
#==============================================================
# Insert perturbation
#==============================================================
def insert_perturbation(x_H, x_v, H_field, B_field, v_field, DATA):
    """
    Insert perturbation
    """

    N_el, local_nodes_H = x_H.shape
    if DATA.test=="Sod" or DATA.test=="Sod_smooth" or DATA.test=="Smooth_periodic":
        if DATA.perturbation!= 0:
            print("No perturbation provided for such test",DATA.test)
            quit()
    elif DATA.test=="Lake_At_Rest_Smooth" or DATA.test=="Lake_At_Rest_Not_Smooth" or DATA.test=="Supercritical_Smooth":
        for inde in range(N_el):
            for indi in range(local_nodes_H):
                if DATA.perturbation==0:
                    pass
                else:
                    x0=6
                    r=0.5

                    if DATA.test=="Lake_At_Rest_Smooth" or DATA.test=="Lake_At_Rest_Not_Smooth":
                        DATA.T=1.5
                    elif DATA.test=="Supercritical_Smooth":
                        DATA.T=1.
                    else:
                        print("No final time set in insert_perturbation in test_dependent module, the test was", DATA.test)
                        quit()

                    if DATA.perturbation==1:
                        H_field[inde,indi]=H_field[inde,indi]+smoot_perturbation(x_H[inde,indi],x0,r,5*10**(-1))
                    elif DATA.perturbation==2:
                        H_field[inde,indi]=H_field[inde,indi]+smoot_perturbation(x_H[inde,indi],x0,r,5*10**(-2))
                    elif DATA.perturbation==3:
                        H_field[inde,indi]=H_field[inde,indi]+smoot_perturbation(x_H[inde,indi],x0,r,5*10**(-3))
                    elif DATA.perturbation==4:
                        H_field[inde,indi]=H_field[inde,indi]+smoot_perturbation(x_H[inde,indi],x0,r,5*10**(-4))
                    elif DATA.perturbation==5:
                        H_field[inde,indi]=H_field[inde,indi]+smoot_perturbation(x_H[inde,indi],x0,r,5*10**(-5))
                    else:
                        print("No perturbation",DATA.perturbation,"available for test",DATA.test)
                        quit()
    else:
        print("Error in test_dependent module, in function function insert_perturbation, test not available")
        quit()
    return H_field, v_field, DATA
#==============================================================
#
#
#
#==============================================================
# Boundary condition
#==============================================================
def BC_state(DATA,x,H_inside, B_inside, v_inside, H_other_side, B_other_side, v_other_side,boundary):
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
    elif DATA.test=="Supercritical_Smooth":
        if boundary=="L":
            H=2.      #H on the left
            q=24      #q on the left
            v=q/H
            B=0.
        else:
            H=H_inside
            B=B_inside
            v=v_inside
    else:
        print("Problems in BC_state, test not available")
        quit()
    return H,B,v
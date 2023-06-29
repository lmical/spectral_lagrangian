import numpy as np
import matplotlib.pyplot as plt
import reference_element
import quadr
import DeC
import test_dependent
import mesh

#==============================================================
#INPUT PARAMETERS
#==============================================================
test               = "Sod"            #Test: "Sod"
N_el               = 30              #Number of elements

#Space
order_space        = 3                #Order in space

#--------------------------------------------------------------
#NB: PGLB basis functions are assumed,
#with associated quadrature providing HO mass lumping
#--------------------------------------------------------------
#type_H             = "PGLB"            #Type of basis functions H
#type_v             = "PGLB"            #Type of basis functions v
#Quadratures, for the moment I'm just going to assume to work with "PGLB" and mass lumping
#...this means that the quadratures are fixed once the order in space is fixed
#...however, one may want to change them in the future
#quadrature_type_H  = "PGLB"            #Quadrature for H
#quadrature_type_v  = "PGLB"            #Quadrature for v
#quadrature_type_HO = "PGLB"            #Quadrature for exact integration
#--------------------------------------------------------------

#Time
time_scheme        = "DeC"             #Time scheme
order_time         = order_space       #Order, only important for arbitrary high order approached like DeC

CFL                = 0.4               #CFL
freq               = 10                #Frequency for saving the solution
N_max_iter         = 10000             #Maximal number of iterations


#Space discretization
scheme             = "Galerkin"
LaxFriedrichs      = False
jump               = "jc"
#==============================================================
#
#
#
#==============================================================
# IF SPECIFIED FROM COMMAND LINE REPLACE THE INPUT PARAMETERS
#==============================================================
# ...
#==============================================================
#
#
#
#==============================================================
# PRINT INFORMATION ON THE SIMULATIONS
#==============================================================
print("------------------------------------------")
print("Starting simulation")
print("Test:", test)
print("Number of elements:", N_el)
print("Order space:", order_space) 
#print("...with basis functions for H and v:", type_H, type_v)
print("Time scheme: ", time_scheme)
if time_scheme=="DeC":
    print("...with order:", order_time)
print("CFL: ", CFL)
print("Frequency for storing the data: ", freq)
print("Maximal number of iterations: ", N_max_iter)
print("------------------------------------------")
#==============================================================
#
#
#
#==============================================================
print("------------------------------------------")
print("Getting GLB nodes, weights and derivatives in the nodes")
#H of degree order_space-1   discontinuous psi
#v of degree order_space     continuous    phi

degree_H=order_space-1
degree_v=order_space

N_local_nodes_H=degree_H+1
N_local_nodes_v=degree_v+1

local_nodes_H, w_H = reference_element.get_nodes(N_local_nodes_H,"gaussLobatto")
local_nodes_v, w_v = reference_element.get_nodes(N_local_nodes_v,"gaussLobatto")

#Local derivatives at the nodes
local_derivatives_H=np.zeros((N_local_nodes_H,N_local_nodes_H))
for indi in range(N_local_nodes_H):
    local_derivatives_H[indi,:] = reference_element.lagrange_deriv(local_nodes_H,local_nodes_H,indi)

local_derivatives_v=np.zeros((N_local_nodes_v,N_local_nodes_v))
for indi in range(N_local_nodes_v):
    local_derivatives_v[indi,:] = reference_element.lagrange_deriv(local_nodes_v,local_nodes_v,indi)




#--------------------------------------------------------------
#I don't know if they are needed, I compute them

#Local_values_H_in_v
local_values_H_in_v=np.zeros((N_local_nodes_H,N_local_nodes_v)) #psi_H(x_v)
for indi in range(N_local_nodes_H):
    local_values_H_in_v[indi,:] = reference_element.lagrange_basis(local_nodes_H,local_nodes_v,indi)


#Local_values_v_in_H
local_values_v_in_H=np.zeros((N_local_nodes_v,N_local_nodes_H)) #phi_v(x_H)
for indi in range(N_local_nodes_v):
    local_values_v_in_H[indi,:] = reference_element.lagrange_basis(local_nodes_v,local_nodes_H,indi)



#Local_derivatives_H_in_v
local_derivatives_H_in_v=np.zeros((N_local_nodes_H,N_local_nodes_v)) #d psi_H(x_v)
for indi in range(N_local_nodes_H):
    local_derivatives_H_in_v[indi,:] = reference_element.lagrange_deriv(local_nodes_H,local_nodes_v,indi)


#Local_derivatives_v_in_H
local_derivatives_v_in_H=np.zeros((N_local_nodes_v,N_local_nodes_H)) #d phi_v(x_H)
for indi in range(N_local_nodes_v):
    local_derivatives_v_in_H[indi,:] = reference_element.lagrange_deriv(local_nodes_v,local_nodes_H,indi)




#--------------------------------------------------------------
# print("order space", order_space)

# print("degree H",degree_H, "which should be equal to", order_space-1)
# print("degree v",degree_v, "which should be equal to",degree_H+1)

# print("N local nodes H",N_local_nodes_H, "which should be equal to", degree_H+1)
# print("N local nodes H",N_local_nodes_v, "which should be equal to", degree_v+1)

# print("local nodes H",local_nodes_H)
# print("local weights H",w_H)

# print("local nodes v",local_nodes_v)
# print("local weights v",w_v)

# print("local derivatives H at the nodes")
# print(local_derivatives_H)

# print("local derivatives v at the nodes")
# print(local_derivatives_v)
# quit()
#--------------------------------------------------------------
# psi=reference_element.basis_functions("PGLB",degree_H)
# for indi in range(len(local_nodes_H)):
#     if np.linalg.norm(psi[indi](local_nodes_v)-local_values_H_in_v[indi,:])>1e-14:
#         print("Problem")
#         quit()

# phi=reference_element.basis_functions("PGLB",degree_v)
# for indi in range(len(local_nodes_v)):
#     if np.linalg.norm(phi[indi](local_nodes_H)-local_values_v_in_H[indi,:])>1e-14:
#         print("Problem")
#         quit()


# dpsi=reference_element.derivative_basis_functions("PGLB",degree_H)
# for indi in range(len(local_nodes_H)):
#     if np.linalg.norm(dpsi[indi](local_nodes_v)-local_derivatives_H_in_v[indi,:])>1e-14:
#         print("Problem")
#         quit()

# dphi=reference_element.derivative_basis_functions("PGLB",degree_v)
# for indi in range(len(local_nodes_v)):
#     if np.linalg.norm(dphi[indi](local_nodes_H)-local_derivatives_v_in_H[indi,:])>1e-14:
#         print("Problem")
#         quit()
#--------------------------------------------------------------


#==============================================================
#
#
#
#==============================================================
print("------------------------------------------")
print("Getting DeC structures")


def subtimesteps(subtimenodes_types,order):
    """
    INPUT:
    subtimenodes_types,  distribution of subtimenodes
    order
    OUTPUT:
    M_subtimenodes, number of subtimenodes-1    
    """
    if subtimenodes_types == "equispaced":
        return order-1
    elif subtimenodes_types=="gaussLobatto":
        return int((order+1)//2)

# Getting M_subtimenodes, i.e., number of subtimenodes - 1
M_subtimenodes=subtimesteps("gaussLobatto",order_time)
# Getting DeC structures
# NB: with theta transposed
dec = DeC.DeC(M_sub=M_subtimenodes, n_iter=order_time, nodes_type="gaussLobatto")

#--------------------------------------------------------------
# print("Number of iterations", dec.n_iter,"which should be", order_time)
# print("Total number of subtimenodes", dec.n_subNodes,"which should be", M_subtimenodes+1,"and",dec.M_sub+1)
# print("Beta vector",dec.beta)
# print("Theta matrix")
# print(dec.theta)
# print("NB: The matrix must be transposed")
# quit()
#--------------------------------------------------------------
#==============================================================
print("------------------------------------------")
print("Getting test information")
DATA=test_dependent.DATA_CLASS(test)
#--------------------------------------------------------------
# print("test",DATA.test)
# print("xL",DATA.xL)
# print("xR",DATA.xR)
# print("periodicity",DATA.periodic)
# quit()
#--------------------------------------------------------------
#==============================================================
#
#
#
#==============================================================
print("------------------------------------------")
print("Mesh initialization")
x_H, x_v, M_Local_to_Global, v_Global_to_Local, N_global_nodes_v, M_faces = mesh.build_mesh(DATA,N_el,local_nodes_H,local_nodes_v)
#----------------------------------------------
# print("Local nodes H", local_nodes_H)
# print("Local nodes v", local_nodes_v)
# print("degree H", degree_H)
# print("degree v", degree_v)
# print("Total DoFs v", N_global_nodes_v)
# print("x_H",x_H)
# print("x_v",x_v)
# print("Local to Global",M_Local_to_Global)
# for indi_g in range(N_global_nodes_v):
#     print("DoF",indi_g)
#     print("Contained by",v_Global_to_Local[indi_g].N_el_containing_node, "elements")
#     print("...and these are",v_Global_to_Local[indi_g].vec_el)
#     print("...and the local DoF in these is",v_Global_to_Local[indi_g].vec_indi_l)
#     print()
# print(M_faces)
#----------------------------------------------
#==============================================================
print("------------------------------------------")
print("Variables initialization")

def Analytical_State(DATA,x,t):
    if DATA.test=="Sod":
        v=0.
        if x<=0.5:
            H=1.
        else:
            H=0.5
    else:
        print("Error in function Analytical_State, test not available")
        quit()
    
    return H,v

def Bathymetry(x):
    if DATA.test=="Sod":
        B=0.
    else:
        print("Error in function Bathymetry, test not available")
        quit()
    return B


def Derivative_Bathymetry(x):
    if DATA.test=="Sod":
        dB=0.
    else:
        print("Error in function Derivative_Bathymetry, test not available")
        quit()
    return dB


def IC(x_H,x_v,DATA,t):
    # Matrix H_field[inde,loc_indi_H], rows=elements, columns=loc_indi_H
    H_field = np.zeros(x_H.shape)
    # Matrix B_field[inde,loc_indi_H], with the same convention
    B_field = np.zeros(x_H.shape)
    # v_field[glob_indi_v]
    v_field = np.zeros((len(x_v)))

    N_el, local_nodes_H = x_H.shape
    if DATA.test=="Sod":
        for inde in range(N_el):
            for indi_l in range(local_nodes_H):
                vec=Analytical_State(DATA,x_H[inde,indi_l],0)
                H_field[inde,indi_l] = vec[0]
                B_field[inde,indi_l] = Bathymetry(x_H[inde,indi_l])
        v_field[:]=0.
    else:
        print("Error in function IC, test not available")
        quit()
    return H_field, B_field, v_field

H_field, B_field, v_field = IC(x_H, x_v, DATA,0)




plt.plot(x_v,v_field)
plt.show()



for inde in range(N_el):
    print(x_H[inde,:])
    plt.plot(x_H[inde,:],H_field[inde,:])
plt.show()
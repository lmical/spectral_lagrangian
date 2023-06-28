import numpy as np
import matplotlib.pyplot as plt
import reference_element
import quadr
import DeC

#==============================================================
#INPUT PARAMETERS
#==============================================================
test               = 0                #Test
N_el               = 100              #Number of elements

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
#---------------------------------------------------------------
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
#---------------------------------------------------------------
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

#---------------------------------------------------------------
# print("Number of iterations", dec.n_iter,"which should be", order_time)
# print("Total number of subtimenodes", dec.n_subNodes,"which should be", M_subtimenodes+1,"and",dec.M_sub+1)
# print("Beta vector",dec.beta)
# print("Theta matrix")
# print(dec.theta)
# print("NB: The matrix must be transposed")
#---------------------------------------------------------------
#==============================================================
#
#
#
#==============================================================
print("------------------------------------------")
print("Mesh and fields initialization")

#Thermodynamic field
print("Initializing matrix x_H[inde,loc_indi_H], rows=elements, columns=loc_indi_H")
x_H     = np.zeros((N_el,N_local_nodes_H))
print("Initializing matrix H_field[inde,loc_indi_H], rows=elements, columns=loc_indi_H")
H_field = np.zeros((N_el,N_local_nodes_H))
print("Initializing matrix B_field[inde,loc_indi_H], with the same convention")
B_field = np.zeros((N_el,N_local_nodes_H))

#Kinetic field
N_global_nodes_v=degree_v*N_el+1
print("Initializing vector x_v[glob_indi_v]")
x_v     = np.zeros((N_global_nodes_v))
print("Initializing vector v_field[glob_indi_v]")
v_field = np.zeros((N_global_nodes_v))

#The kinetic field is global, hence, it is useful to have some connectivity structures



print("Initializing matrix M_Local_to_Global[inde,loc_indi_v], rows=elements, columns=loc_indi_H")
print("content=Global index associated to the local node loc_indi_v in the element inde")
M_Local_to_Global=np.zeros((N_el,N_local_nodes_v))
print("Initializing Matrix M_Global_to_Local[glob_indi_v,0:1], rows=global index glob_indi_v, column_0=element containing it, column_1=corresponding local index loc_indi_v in the element")
M_Global_to_Local=np.zeros((N_global_nodes_v,2))

#---------------------------------------------------------------
# print()
# print("Number of elements", N_el)
# print("Local DoFs H", N_local_nodes_H)
# print("Local DoFs v", N_local_nodes_v)
# print()
# print("Order space",order_space)
# print("Degree H"   ,degree_H)
# print("Degree v"   ,degree_v)
# print()
# print("Size of x_H"      ,x_H.shape)
# print("Size of H_field"  ,H_field.shape)
# print("Length of x_v"    ,len(x_v))
# print("Length of v_field",len(v_field))
# print()
# print("Total DoFs v",N_global_nodes_v)
# print("Size of M_Local_to_Global",M_Local_to_Global.shape)
# print("Size of M_Local_to_Global",M_Global_to_Local.shape)
# print()
#---------------------------------------------------------------
#==============================================================
#
#
#
#==============================================================
print("------------------------------------------")
print("Getting test information")

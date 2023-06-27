import numpy as np
import matplotlib.pyplot as plt
import reference_element
import quadr


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
print("order space", order_space)

print("degree H",degree_H, "which should be equal to", order_space-1)
print("degree v",degree_v, "which should be equal to",degree_H+1)

print("N local nodes H",N_local_nodes_H, "which should be equal to", degree_H+1)
print("N local nodes H",N_local_nodes_v, "which should be equal to", degree_v+1)

print("local nodes H",local_nodes_H)
print("local weights H",w_H)

print("local nodes v",local_nodes_v)
print("local weights v",w_v)

print("local derivatives H at the nodes")
print(local_derivatives_H)

print("local derivatives v at the nodes")
print(local_derivatives_v)
#---------------------------------------------------------------


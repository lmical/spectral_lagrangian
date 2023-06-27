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
#H of degree M   discontinuous psi
#v of degree M+1 continuous    phi

M=order_space-1    






xi,w=reference_element.get_nodes(3,"gaussLobatto")
print(reference_element.get_nodes(3,"gaussLobatto"))
print(xi,w)
print(reference_element.lagrange_basis(xi,xi,0))
print(reference_element.lagrange_basis(xi,xi,1))
print(reference_element.lagrange_basis(xi,xi,2))



quit()

degree_H = order_space-1
degree_v = degree_H+1

phiH     = reference_element.basis_functions(type_H,degree_H)
psiv     = reference_element.basis_functions(type_v,degree_v)

#VECTORIZE
x=np.linspace(0,1,100)
y=psiv[3](x)
plt.plot(x,y)
plt.grid()
plt.show()


print("Initialization of the quadrature points")


import numpy as np


#==============================================================
#INPUT PARAMETERS
#==============================================================
test               = 0                #Test
N_el               = 100              #Number of elements

#Space
order_space        = 2                #Order in space
type_H             = "PGL"            #Type of basis functions H
type_v             = "PGL"            #Type of basis functions v
#Quadratures, for the moment I'm just going to assume to work with "PGL" and mass lumping
#...this means that the quadratures are fixed once the order in space is fixed
#...however, one may want to change them in the future
#quadrature_type_H  = "PGL"            #Quadrature for H
#quadrature_type_v  = "PGL"            #Quadrature for v
#quadrature_type_HO = "PGL"            #Quadrature for exact integration


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
#IF SPECIFIED FROM COMMAND LINE REPLACE THE INPUT PARAMETERS
#==============================================================
#...
#==============================================================
#
#
#
#==============================================================
#PRINT INFORMATION ON THE SIMULATIONS
#==============================================================
print("------------------------------------------")
print("Starting simulation")
print("Test:", test)
print("Number of elements:", N_el)
print("Order space:", order_space) 
print("...with basis functions for H and v:", type_H, type_v)
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
print("Initialization of the basis functions")





print("Initialization of the quadrature points")


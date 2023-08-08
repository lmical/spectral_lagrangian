import numpy as np
import matplotlib.pyplot as plt
import reference_element
import quadr
import DeC
import test_dependent
import mesh
import lagrangian_scheme
import output
import time_stepping
import sys

#==============================================================
#INPUT PARAMETERS
#==============================================================
test               = "Smooth_periodic"     #Test: "Sod", "Sod_smooth", "Smooth_periodic", 
                                                #"Lake_At_Rest_Smooth", "Lake_At_Rest_Not_Smooth"
                                                #"Supercritical_Smooth", "Supercritical_Not_Smooth"
                                                #"Constant_Slope_Smooth"


perturbation       = 0                          #Perturbation

N_el               = 50                        #Number of elements

#Space
order_space        = 2                         #Order in space

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
time_scheme        = "DeC"             #Time scheme #"Euler" "DeC" "SSPRK4"
order_time         = order_space       #Order, only important for arbitrary high order approached like DeC

CFL                = 0.5               #CFL
freq               = 800                #Frequency for saving the solution
N_max_iter         = 1000000             #Maximal number of iterations


#Space discretization
scheme             = "Galerkin"
LaxFriedrichs      = True
WB                 = False
jump               = "jc"               #j0,    jc

#Folder where to store
folder             = "Debug"
printing           = True
plotting           = False
storing            = True


#==============================================================
#
#
#
#==============================================================
# IF SPECIFIED FROM COMMAND LINE REPLACE THE INPUT PARAMETERS
#==============================================================
# print(sys.argv,len(sys.argv))
if len(sys.argv)>1:
    test=sys.argv[1]
if len(sys.argv)>2:
    perturbation=int(sys.argv[2])
if len(sys.argv)>3:
    order_space=int(sys.argv[3])
if len(sys.argv)>4:
    time_scheme=sys.argv[4]
    if time_scheme=="DeC":
        order_time=order_space
if len(sys.argv)>5:
    if sys.argv[5]=="True":
        LaxFriedrichs=True
    elif sys.argv[5]=="False":
        LaxFriedrichs=False
    else:
        print("Impossible to get LxF imput from keyboard")
        quit()
if len(sys.argv)>6:
    jump=sys.argv[6]
if len(sys.argv)>7:
    CFL=float(sys.argv[7])
if len(sys.argv)>8:
    N_el=int(sys.argv[8])


#==============================================================
#
#
#
#==============================================================
# PRINT INFORMATION ON THE SIMULATIONS
#==============================================================
print("------------------------------------------")
print("Starting simulation")
print("Test:", test, "with perturbation", perturbation)
print("Number of elements: ",N_el)
print("Order space:", order_space,"->","P"+str(order_space-1)+"P"+str(order_space)) 
#print("...with basis functions for H and v:", type_H, type_v)
print("Time scheme: ", time_scheme)
if time_scheme=="DeC":
    print("...with order:", order_time)
print("Lax-Friedrichs: ", LaxFriedrichs)
print("Jump: ", jump)
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
if time_scheme=="DeC":
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
DATA=test_dependent.DATA_CLASS(test,perturbation,N_el,order_space,time_scheme,order_time,CFL,freq,N_max_iter,scheme,LaxFriedrichs,WB,jump,folder,printing,plotting,storing)
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
#
#
#
#==============================================================
print("------------------------------------------")
print("Variables initialization")
H_field, B_field, v_field = test_dependent.IC(x_H, x_v, 0, DATA)
H_field, v_field, DATA = test_dependent.insert_perturbation(x_H, x_v, H_field, B_field, v_field, DATA)

#Getting the field on the reference element for strong mass conservation
# Hhat_i = H_i(0)*det J(x_i,0)
# J(xi,0) = grad_xi x (xi,0) = sum_j x_j(0) grad_xi phi_j (xi)
Hhat_field=lagrangian_scheme.get_Hhat_on_reference_element(H_field,x_v,local_derivatives_v_in_H,M_Local_to_Global)
#----------------------------------------------
# plt.plot(x_v,v_field)
# plt.show()
# for inde in range(N_el):
#     plt.plot(x_H[inde,:],H_field[inde,:])
# plt.show()
#----------------------------------------------
# H_field_test=lagrangian_scheme.strong_mass_conservation(Hhat_field,x_v,local_derivatives_v_in_H,M_Local_to_Global)
# if (np.linalg.norm(H_field-H_field_test)>1e-14):
#     print("I'm in main. Problem in strong_mass_conservation")
#     quit()
#==============================================================
#
#
#
#==============================================================
# Testing some structures
# M_v=lagrangian_scheme.Lumped_Mass_Matrix(w_v,x_v,M_Local_to_Global,local_derivatives_v)
# print(M_v)
#----------------------------------------------
# phi_v=lagrangian_scheme.Space_Residuals_v(H_field, B_field, w_v,local_derivatives_H_in_v,M_Local_to_Global)
# print(phi_v)
#----------------------------------------------
# ST_i=lagrangian_scheme.Lax_Friedrichs(v_field,M_Local_to_Global)
# print(ST_i)
#----------------------------------------------
# CT_phi_v=lagrangian_scheme.Coupling_Terms_Space_Residuals_v(H_field, B_field, v_field, M_Local_to_Global, M_faces, DATA)
# print(CT_phi_v)
#----------------------------------------------
# x_H = lagrangian_scheme.get_x_H(x_v,local_values_v_in_H,M_Local_to_Global)
# print(x_H)
#==============================================================
#
#
#
#==============================================================
print("------------------------------------------")
print("Timestepping loop")
DATA.time=0     #Time
indt=0  #Counter

if (LaxFriedrichs==True) and (test=="Smooth_periodic" or test=="Lake_At_Rest_Smooth" or test=="Supercritical_Smooth" or test=="Constant_Slope_Smooth"):
    print()
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print("Warning, running a smooth test with first order limiting!")
    print("You may want to really do it, I'm telling you just in case!")
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
    print()


#Printing and plotting IC
if printing==True:
    output.printing_function(indt,DATA.time,H_field,v_field)
if plotting==True:
    x_H=lagrangian_scheme.get_x_H(x_v,local_values_v_in_H,M_Local_to_Global) #Not necessary, but called for coherence
    H_in_x_v=lagrangian_scheme.get_H_in_x_v(H_field,x_v,local_values_H_in_v,M_Local_to_Global)
    output.plotting_function(indt,x_H,H_field,B_field,x_v,v_field,H_in_x_v,DATA,storing_info=False)

while(DATA.time<DATA.T):
    #Computation of the time step
    dt=DATA.T-DATA.time
    dt_max=lagrangian_scheme.Compute_Time_Step(H_field,v_field,x_v,M_Local_to_Global,DATA,degree_v,CFL)
    DATA.dt=min(dt,dt_max)

    #Store solution in the previous time step
    x_v_old=x_v
    H_field_old=H_field
    v_field_old=v_field
    B_field_old=B_field


    if time_scheme=="Euler":
        H_field, v_field, x_v, B_field=time_stepping.Euler_method(H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA)
    elif time_scheme=="DeC":
        H_field, v_field, x_v, B_field=time_stepping.DeC_method(H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA, dec)
    elif time_scheme=="SSPRK4":
        H_field, v_field, x_v, B_field=time_stepping.SSPRK4_method(H_field_old, v_field_old, x_v_old, B_field_old, Hhat_field, w_v, local_derivatives_v, local_derivatives_H_in_v, local_derivatives_v_in_H, M_Local_to_Global, local_values_v_in_H, M_faces, DATA)
    else:
        print("Time scheme not available")
        quit()



    #t
    DATA.time=DATA.time+DATA.dt
    indt=indt+1

    #Only every freq timesteps
    if (indt%freq==0):
        #Printing and plotting IC
        if printing==True:
            output.printing_function(indt,DATA.time,H_field,v_field)
        if plotting==True:
            x_H=lagrangian_scheme.get_x_H(x_v,local_values_v_in_H,M_Local_to_Global)
            H_in_x_v=lagrangian_scheme.get_H_in_x_v(H_field,x_v,local_values_H_in_v,M_Local_to_Global)
            output.plotting_function(indt,x_H,H_field,B_field,x_v,v_field,H_in_x_v,DATA,storing_info=False)


    if indt>=N_max_iter:
        print("Total number of iterations reached", indt, N_max_iter)
        break

#Final print
if printing==True:
    output.printing_function(indt,DATA.time,H_field,v_field)



if storing==True:
    #Final plot to save
    x_H=lagrangian_scheme.get_x_H(x_v,local_values_v_in_H,M_Local_to_Global)
    H_in_x_v=lagrangian_scheme.get_H_in_x_v(H_field,x_v,local_values_H_in_v,M_Local_to_Global)
    output.plotting_function(indt,x_H,H_field,B_field,x_v,v_field,H_in_x_v,DATA,storing_info=True)

    #Storing the final solution
    #indi, x_v, v, H, q, eta
    B_in_x_v=lagrangian_scheme.get_H_in_x_v(B_field,x_v,local_values_H_in_v,M_Local_to_Global)
    output.storing(H_field, v_field, x_v, B_field, H_in_x_v, B_in_x_v, M_Local_to_Global, DATA)

    #Compute error
    if DATA.analytical_solution==True and DATA.perturbation==0:
        output.compute_error(H_field, v_field, x_v, x_H, H_in_x_v, M_Local_to_Global, w_H, w_v, local_derivatives_v_in_H, local_derivatives_v, DATA)
        output.plot_error(H_field, v_field, x_v, x_H, H_in_x_v, DATA)
        

print(test,"N_el",N_el,"order_space",order_space,"CFL",CFL)
print("Maxima",np.max(H_field),np.max(v_field))
print("Minima",np.min(H_field),np.min(v_field))
print("Average",np.average(H_field),np.average(v_field))



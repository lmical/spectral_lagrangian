import numpy as np
import matplotlib.pyplot as plt

#==============================================================
#INPUT PARAMETERS
#==============================================================
test               = 0                #Test
N_el               = 100              #Number of elements

#Space
order_space        = 3                #Order in space
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
print("Starting simulation")
print("Test:", test)
print("Number of elements:", N_el)
print("Order space:", order_space) 
print("...with basis functions for H and v:", type_H, type_v)
print("Time scheme: ", time_scheme)
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
print("Initialization of the basis functions in the reference element")
#M   discontinuous psi
#M+1 continuous    phi

def basis_functions(type_basis,degree):
    #INPUT:
    #type_basis
    #degree
    #OUTPUT
    #array phi of functions 


    phi=np.array([])

    if type_basis=="PGL":

        if degree==0: #---0---#

            def phi0(xi):
            	return 1
            phi0=np.vectorize(phi0)
            phi=np.append(phi,phi0)

        elif degree==1: #0------1#

            def phi0(xi):
            	return 1-xi
            phi0=np.vectorize(phi0)
            phi=np.append(phi,phi0)

            def phi1(xi):
            	return xi
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

        elif degree==2: #0---2---1#

            #x1=xi     barycentric coordinate of the right node
            #x2=1-xi   barycentric coordinate of the left node
 
            def phi0(xi):
            	x1=xi
            	x2=1-xi
            	return x2*(x2-x1)
            phi0=np.vectorize(phi0)
            phi=np.append(phi,phi0)

            def phi1(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-x2)
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

            def phi2(xi):
            	x1=xi
            	x2=1-xi
            	return 4*x1*x2
            phi2=np.vectorize(phi2)
            phi=np.append(phi,phi2)


        elif degree==3: #0--2--3--1#

            #x1=xi     barycentric coordinate of the right node
            #x2=1-xi   barycentric coordinate of the left node
 
            alpha = 0.5-np.sqrt(5)/10    
            beta  = 1-alpha
 
            def phi0(xi):
            	x1=xi
            	x2=1-xi
            	return  (-x1+alpha)*(x1-beta)*(x1-1)/(alpha*beta)
            phi0=np.vectorize(phi0)
            phi=np.append(phi,phi0)

            def phi1(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-alpha)*(x1-beta)/alpha/beta
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

            def phi2(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-beta)*(x1-1)/alpha/(-2*alpha+1)/beta
            phi2=np.vectorize(phi2)
            phi=np.append(phi,phi2)


            def phi3(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-alpha)*(x1-1)/alpha/(2*alpha-1)/beta
            phi3=np.vectorize(phi3)
            phi=np.append(phi,phi3)


        elif degree==4: #0--2-3-4--1#

            #x1=xi     barycentric coordinate of the right node
            #x2=1-xi   barycentric coordinate of the left node
  
            alpha=0.5-np.sqrt(21)/14
            beta = 1 -alpha
  
            def phi0(xi):
            	x1=xi
            	x2=1-xi
            	return (2*(alpha-x1)*(x1-1)*(x1-0.5)*(x1-beta))/(-alpha*beta)
            phi0=np.vectorize(phi0)
            phi=np.append(phi,phi0)

            def phi1(xi):
            	x1=xi
            	x2=1-xi
            	return (2*x1*(-alpha + x1)*(x1 - 0.5)*(x1 - beta))/(alpha*beta)
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

            def phi2(xi):
            	x1=xi
            	x2=1-xi
            	return (x1*(x1-1)*(x1-0.5)*(x1-beta))/(-alpha*(2*alpha-1)*beta*(alpha-0.5))
            phi2=np.vectorize(phi2)
            phi=np.append(phi,phi2)


            def phi3(xi):
            	x1=xi
            	x2=1-xi
            	return  -(4*x1*(alpha -x1)*(x1-1)*(x1-beta))/(alpha-0.5)**2
            phi3=np.vectorize(phi3)
            phi=np.append(phi,phi3)

            def phi4(xi):
            	x1=xi
            	x2=1-xi
            	return -(x1*(alpha-x1)*(x1-1)*(x1-0.5))/(-alpha*(2*alpha-1)*beta*(alpha-0.5))
            phi4=np.vectorize(phi4)
            phi=np.append(phi,phi4)
        else:
            print("Error in basis_functions, degree not defined",type_basis,degree)
            quit()
        
    else:
    	print("Error in basis_functions, basis not defined",type_basis,degree)
    	quit()
    	
    return phi
    
degree_H = order_space-1
degree_v = degree_H+1

phiH     = basis_functions(type_H,degree_H)
psiv     = basis_functions(type_v,degree_v)

#VECTORIZE
x=np.linspace(0,1,100)
y=psiv[3](x)
plt.plot(x,y)
plt.grid()
plt.show()


print("Initialization of the quadrature points")


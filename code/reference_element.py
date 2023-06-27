import numpy as np
import matplotlib.pyplot as plt
import quadr
from mpmath import mp



#==============================================================
# Function to get the (GLB) nodes and weights, w, in the reference interval [0,1]
#==============================================================
def get_nodes(n_points,nodes_type):
    #INPUT
    # n_points
    # nodes_type
    #OUTPUT
    # nodes,  quadrature nodes in [0,1] 
    # w,      weights

    # Added to deal with P0
    if n_points==1:
        nodes=np.zeros(1)
        w=np.zeros(1)
        nodes[0]=1.
        w[0]=1.
        return nodes, w

    # Standard routine
    if nodes_type=="equispaced":
        nodes,w = quadr.equispaced(n_points)
    #elif nodes_type == "gaussLegendre": #Removed because not needed
    #    nodes,w = leggauss(order)
    elif nodes_type == "gaussLobatto":
        nodes, w = quadr.lglnodes(n_points-1,10**-15)
    nodes=nodes*0.5+0.5
    w = w*0.5
    return nodes, w
#==============================================================
#
#
#
#==============================================================
#Function to get the evaluation of Lagrange polynomials associated to certain nodes
#in certain points x
#In, particular the evaluations of l_k are given
#==============================================================
def lagrange_basis(nodes,x,k):
    #INPUT
    # nodes, nodal points
    # x,     vector where to evaluate the Lagrangian function l_k
    # k,     index of the Lagrangian function to be avaluated
    #OUTPUT
    # y,     vector of evaluations of l_k in the abscissae in x


    y=np.zeros(x.size)
    for ix, xi in enumerate(x):
        tmp=[(xi-nodes[j])/(nodes[k]-nodes[j])  for j in range(len(nodes)) if j!=k]
        y[ix]=np.prod(tmp)
    return y
#==============================================================
#
#
#
#==============================================================
#Function to get the derivative of Lagrange polynomials associated to certain nodes
#in certain points x
#In, particular the evaluations of d l_k are given
#==============================================================
# NB: Derivatives compute in multiple precision and then turned into float
#==============================================================
def lagrange_deriv(nodes,x,k):
    f=mp.mpf(0) 
    for j in range(len(nodes)):
        p=mp.mpf(1)
        if k!=j:
            for l in range(len(nodes)):
                if l!=k and l!=j: 
                    p=p*(x-nodes[l])/(nodes[k]-nodes[l])
            f = f + p/(nodes[k]-nodes[j])
    f=convert_vector_mp_2_np(f)
    return f
#==============================================================
#
#
#
#==============================================================
def convert_vector_mp_2_np(v):
    vnp = np.zeros(len(v))
    for i in range(len(v)):
        vnp[i]=float(v[i])
    return vnp
#==============================================================
#
#
#
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
# FUNCTIONS FROM HERE ON USED JUST FOR VERIFICATION
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#==============================================================
#
#
#
#==============================================================
# basis functions
#==============================================================
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
    


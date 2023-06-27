import numpy as np
import matplotlib.pyplot as plt
import quadr
from mpmath import mp



#==============================================================
# Function to get the (GLB) nodes and weights, w, in the reference interval [0,1]
#==============================================================
def get_nodes(n_points,nodes_type):
    """
    INPUT
     n_points
     nodes_type
    OUTPUT
     nodes,  quadrature nodes in [0,1] 
     w,      weights
    """
    
    # Added to deal with P0
    if n_points==1:
        nodes=np.zeros(1)
        w=np.zeros(1)
        nodes[0]=0.5
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
    """
    INPUT
     nodes, nodal points
     x,     vector where to evaluate the Lagrangian function l_k
     k,     index of the Lagrangian function to be evaluated
    OUTPUT
     y,     vector of evaluations of l_k in the abscissae in x
    """

    y=np.zeros(len(x))
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
    """
    INPUT
     nodes, nodal points
     x,     vector where to evaluate the derivative of the Lagrangian function l_k
     k,     index of the Lagrangian function, whose derivative has to be evaluated
    OUTPUT
     y,     vector of evaluations of l_k in the abscissae in x
    """

    # Added to deal with P0
    if len(nodes)==1:
	    return np.zeros(len(x))
    
    y=np.zeros(len(x))
    for indi in range(len(x)):
        f=mp.mpf(0) 
        for j in range(len(nodes)):
            p=mp.mpf(1)
            if k!=j:
                for l in range(len(nodes)):
                    if l!=k and l!=j: 
                        p=p*(x[indi]-nodes[l])/(nodes[k]-nodes[l])
                f = f + p/(nodes[k]-nodes[j])
        f=float(f)
        y[indi]=f
    return y
#==============================================================
#
#
#
#==============================================================
def convert_vector_mp_2_np(v):
    """
    From arbitrary precision to numpy array
    """
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
    """
    INPUT:
    type_basis
    degree
    OUTPUT
    array phi of functions 
    """

    phi=np.array([])

    if type_basis=="PGLB":

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

        elif degree==2: #0---1---2#

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
            	return 4*x1*x2
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

            def phi2(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-x2)
            phi2=np.vectorize(phi2)
            phi=np.append(phi,phi2)


        elif degree==3: #0--1--2--3#

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
            	return x1*(x1-beta)*(x1-1)/alpha/(-2*alpha+1)/beta
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

            def phi2(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-alpha)*(x1-1)/alpha/(2*alpha-1)/beta 
            phi2=np.vectorize(phi2)
            phi=np.append(phi,phi2)


            def phi3(xi):
            	x1=xi
            	x2=1-xi
            	return x1*(x1-alpha)*(x1-beta)/alpha/beta
            phi3=np.vectorize(phi3)
            phi=np.append(phi,phi3)


        elif degree==4: #0--1-2-3--4#

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
            	return (x1*(x1-1)*(x1-0.5)*(x1-beta))/(-alpha*(2*alpha-1)*beta*(alpha-0.5))
            phi1=np.vectorize(phi1)
            phi=np.append(phi,phi1)

            def phi2(xi):
            	x1=xi
            	x2=1-xi
            	return -(4*x1*(alpha -x1)*(x1-1)*(x1-beta))/(alpha-0.5)**2
            phi2=np.vectorize(phi2)
            phi=np.append(phi,phi2)


            def phi3(xi):
            	x1=xi
            	x2=1-xi
            	return  -(x1*(alpha-x1)*(x1-1)*(x1-0.5))/(-alpha*(2*alpha-1)*beta*(alpha-0.5))
            phi3=np.vectorize(phi3)
            phi=np.append(phi,phi3)

            def phi4(xi):
            	x1=xi
            	x2=1-xi
            	return (2*x1*(-alpha + x1)*(x1 - 0.5)*(x1 - beta))/(alpha*beta)
            phi4=np.vectorize(phi4)
            phi=np.append(phi,phi4)
        else:
            print("Error in basis_functions, degree not defined",type_basis,degree)
            quit()
        
    else:
    	print("Error in basis_functions, basis not defined",type_basis,degree)
    	quit()
    	
    return phi
#==============================================================
#
#
#
#==============================================================
# derivative basis functions
#==============================================================    
def derivative_basis_functions(type_basis,degree):
    """
    INPUT:
    type_basis
    degree
    OUTPUT
    array dphi of functions 
    """

    dphi=np.array([])

    if type_basis=="PGLB":

        if degree==0: #---0---#

            def dphi0(xi):
            	return 0
            dphi0=np.vectorize(dphi0)
            dphi=np.append(dphi,dphi0)

        elif degree==1: #0------1#

            def dphi0(xi):
            	return -1
            dphi0=np.vectorize(dphi0)
            dphi=np.append(dphi,dphi0)

            def dphi1(xi):
            	return 1
            dphi1=np.vectorize(dphi1)
            dphi=np.append(dphi,dphi1)

        elif degree==2: #0---1---2#

            #x1=xi     barycentric coordinate of the right node
            #x2=1-xi   barycentric coordinate of the left node
 
            def dphi0(xi):
            	x1=xi
            	x2=1-xi
            	return -3+4*x1
            dphi0=np.vectorize(dphi0)
            dphi=np.append(dphi,dphi0)

            def dphi1(xi):
            	x1=xi
            	x2=1-xi
            	return 4-8*x1
            dphi1=np.vectorize(dphi1)
            dphi=np.append(dphi,dphi1)

            def dphi2(xi):
            	x1=xi
            	x2=1-xi
            	return 4*x1-1
            dphi2=np.vectorize(dphi2)
            dphi=np.append(dphi,dphi2)

        elif degree==3: #0--2--3--1#

            #x1=xi     barycentric coordinate of the right node
            #x2=1-xi   barycentric coordinate of the left node
 
            alpha = 0.5-np.sqrt(5)/10    
            beta  = 1-alpha
 
            def dphi0(xi):
            	x1=xi
            	x2=1-xi
            	return  (- alpha**2 + alpha + 3*x1**2 - 4*x1 + 1)/(alpha*(alpha - 1))
            dphi0=np.vectorize(dphi0)
            dphi=np.append(dphi,dphi0)

            def dphi1(xi):
            	x1=xi
            	x2=1-xi
            	return (2*alpha*x1 - 4*x1 - alpha + 3*x1**2 + 1)/(alpha*(2*alpha**2 - 3*alpha + 1))
            dphi1=np.vectorize(dphi1)
            dphi=np.append(dphi,dphi1)

            def dphi2(xi):
            	x1=xi
            	x2=1-xi
            	return -(alpha - 2*x1 - 2*alpha*x1 + 3*x1**2)/(alpha*(2*alpha**2 - 3*alpha + 1))
            dphi2=np.vectorize(dphi2)
            dphi=np.append(dphi,dphi2)


            def dphi3(xi):
            	x1=xi
            	x2=1-xi
            	return -(- alpha**2 + alpha + 3*x1**2 - 2*x1)/(alpha*(alpha - 1))
            dphi3=np.vectorize(dphi3)
            dphi=np.append(dphi,dphi3)

        elif degree==4: #0--2-3-4--1#

            #x1=xi     barycentric coordinate of the right node
            #x2=1-xi   barycentric coordinate of the left node
  
            alpha=0.5-np.sqrt(21)/14
            beta = 1 -alpha
  
            def dphi0(xi):
            	x1=xi
            	x2=1-xi
            	return (4*alpha**2*x1 - 3*alpha**2 - 4*alpha*x1 + 3*alpha - 8*x1**3 + 15*x1**2 - 8*x1 + 1)/(alpha*(alpha - 1))
            dphi0=np.vectorize(dphi0)
            dphi=np.append(dphi,dphi0)

            def dphi1(xi):
            	x1=xi
            	x2=1-xi
            	return (alpha + 8*x1 - 6*alpha*x1 + 6*alpha*x1**2 - 15*x1**2 + 8*x1**3 - 1)/(alpha*(2*alpha - 1)**2*(alpha - 1))
            dphi1=np.vectorize(dphi1)
            dphi=np.append(dphi,dphi1)

            def dphi2(xi):
            	x1=xi
            	x2=1-xi
            	return (16*(2*x1 - 1)*(- alpha**2 + alpha + 2*x1**2 - 2*x1))/(2*alpha - 1)**2
            dphi2=np.vectorize(dphi2)
            dphi=np.append(dphi,dphi2)


            def dphi3(xi):
            	x1=xi
            	x2=1-xi
            	return -(alpha - 2*x1 - 6*alpha*x1 + 6*alpha*x1**2 + 9*x1**2 - 8*x1**3)/(alpha*(2*alpha - 1)**2*(alpha - 1))
            dphi3=np.vectorize(dphi3)
            dphi=np.append(dphi,dphi3)

            def dphi4(xi):
            	x1=xi
            	x2=1-xi
            	return -(- 4*alpha**2*x1 + alpha**2 + 4*alpha*x1 - alpha + 8*x1**3 - 9*x1**2 + 2*x1)/(alpha*(alpha - 1))
            dphi4=np.vectorize(dphi4)
            dphi=np.append(dphi,dphi4)
        else:
            print("Error in derivative_basis_functions, degree not defined",type_basis,degree)
            quit()

    else:
    	print("Error in derivative_basis_functions, basis not defined",type_basis,degree)
    	quit()
    	
    return dphi
#==============================================================
#
#
#
#==============================================================
# GLB nodes and weights
#==============================================================
def GLB_nodes_weights(n_points):
    """
    INPUT:
    n_points
    OUTPUT:
    nodes
    weights
    """

    nodes = np.zeros(n_points)
    w     = np.zeros(n_points)
    if n_points==1:

        nodes[0] = 0.5
        w[0]     = 1.

    elif n_points==2:

        nodes[0] = 0
        w[0]     = 0.5	
        
        nodes[1] = 1.
        w[1]     =  0.5

    elif n_points==3:

        nodes[0] = 0.
        w[0]     = 1./6.

        nodes[1] = 0.5
        w[1]     = 2./3.
 
        nodes[2] = 1.
        w[2]     = 1./6.

    elif n_points==4:


       s=0.5-np.sqrt(5)/10

       nodes[0] = 0.
       w[0]     = 1./12


       nodes[1] = s
       w[1]     = 5/12

       nodes[2] = 1. -s
       w[2]     = 5./12.

       nodes[3] = 1.
       w[3]     = 1./12.


    elif n_points==5:

       s=0.5-np.sqrt(21.)/14.

       nodes[0] = 0.
       w[0]     = 1./20.


       nodes[1] = s
       w[1]     = 49./180.

       nodes[2] = 0.5
       w[2]     = 16./45.

       nodes[3] = 1. -s
       w[3]     = 49./180.

       nodes[4] = 1.
       w[4]     = 1./20.

    else:
        print("Erroe in GLB_nodes_weights, n_points not available",n_points)
        quit()
    return nodes,w


#==============================================================
# #TEST: compare basis functions and derivative to hard-coded
#==============================================================
# print("Test basis functions")   
# x=np.linspace(0,1,100)
# for degree in range(5):
#     phi=basis_functions("PGLB",degree)
#     nodes,w=get_nodes(degree+1,"gaussLobatto")
#     for indi in range(len(nodes)):
#         y=phi[indi](x)
#         z=lagrange_basis(nodes,x,indi)
#         plt.plot(x,y)
#         plt.plot(x,z)
#         print("degree",degree,"basis",indi,"error in comparison",np.linalg.norm(y-z))
#         if (np.linalg.norm(y-z)>1e-14):    
#             print("Problem",degree,indi)
#             quit()
#         plt.grid()
#         plt.show()

# print("Test derivative basis functions")   
# x=np.linspace(0,1,100)
# for degree in range(5):
#     dphi=derivative_basis_functions("PGLB",degree)
#     nodes,w=get_nodes(degree+1,"gaussLobatto")
#     for indi in range(len(nodes)):
#         y=dphi[indi](x)
#         z=lagrange_deriv(nodes,x,indi)
#         plt.plot(x,y)
#         plt.plot(x,z)
#         print("degree",degree,"basis",indi,"error in comparison",np.linalg.norm(y-z))
#         if (np.linalg.norm(y-z)>1e-13):    
#             print("Problem",degree,indi)
#             quit()
#         plt.grid()
#         plt.show()

print("Test GLB quadrature")   
for degree in range(5):
    x,y=GLB_nodes_weights(degree+1)
    nodes,w=get_nodes(degree+1,"gaussLobatto")
    for indi in range(len(nodes)):
        print("degree",degree,"node",indi,"error in comparison",np.linalg.norm(x-nodes)+np.linalg.norm(y-w))
        if (np.linalg.norm(x-nodes)+np.linalg.norm(y-w)>1e-14):    
            print("Problem",degree,indi)
            print("nodes",x)
            print("nodes",nodes)
            print("weights",y)
            print("weights",w)
            quit()

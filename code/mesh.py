import numpy as np

class node_Global_to_Local_class:
    def __init__(self):
        self.N_el_containing_node = 0            #Number of elements containing the DoF
        self.vec_el               = np.array([],dtype=int) #Indices of the elements containing the DoF
        self.vec_indi_l           = np.array([],dtype=int) #Local indices of the DoF in the elements

def build_mesh(DATA,N_el,local_nodes_H,local_nodes_v):
              #DATA,N_el,local_nodes_H,local_nodes_v
              #DATA,N_el,local_nodes_H,local_nodes_v

    #Reconstruct some informations from the inputs
    N_local_nodes_H  = len(local_nodes_H)
    N_local_nodes_v  = len(local_nodes_v)
 
    degree_H         = N_local_nodes_H-1
    degree_v         = N_local_nodes_v-1



    N_global_nodes_v = degree_v*N_el+1

    print(N_local_nodes_H,N_local_nodes_v)
    print(degree_H,degree_v)
    print(N_global_nodes_v)


    #---------------------------------------
    #Initializing the variables
    #---------------------------------------

    # Thermodynamic field
    # Matrix x_H[inde,loc_indi_H], rows=elements, columns=loc_indi_H
    x_H     = np.zeros((N_el,N_local_nodes_H))



    # Kinetic field
    N_global_nodes_v=degree_v*N_el+1
    # Vector x_v[glob_indi_v]
    x_v     = np.zeros((N_global_nodes_v))


    # The kinetic field is global, hence, it is useful to have some connectivity structures

    # Local         -> Global
    # (inde,l_indi) -> g_indi
    # Matrix M_Local_to_Global[inde,loc_indi_v], rows=elements, columns=loc_indi_v
    # content = Global index associated to the local node loc_indi_v in the element inde
    M_Local_to_Global=np.zeros((N_el,N_local_nodes_v),dtype=int)




    # Global       -> Local
    # g_indi       -> [(inde,l_indi),...,(inde,l_indi)]
    # vector v_Global_to_Local[glob_indi_v]
    # content=vector of vectors of the type [inde,loc_indi_v], inde=element containing the global DoF, loc_indi_v=local index in the element
    v_Global_to_Local=np.array([])



    #NB: I always assume the DoFs orderd by increasing abscissa, locally and globally
    
    #---------------------------------------
    #Filling the variables
    #---------------------------------------

    x_interfaces=np.linspace(DATA.xL,DATA.xR,N_el+1)
    dx=x_interfaces[1]-x_interfaces[0]





    for inde in range(N_el):
        x_H[inde,:]=x_interfaces[inde]+dx*local_nodes_H





    indi_g=0 #counter on the global nodes
    x_v[0]=DATA.xL
    for inde in range(N_el):

        for indi_l in range(1,len(local_nodes_v)): #Loop on the local DoFs excluding the first one
            indi_g=indi_g+1
            x_v[indi_g]=x_interfaces[inde]+dx*local_nodes_v[indi_l]

            print(inde,indi_l)
        print()

    print(x_H.shape)    
    print(len(x_v))
    print(x_interfaces)
    print(dx)
    print(x_H)
    print(x_v)








    indi_g=0 #counter on the global nodes
    for inde in range(N_el):

        for indi_l in range(len(local_nodes_v)): #Loop on the local DoFs

            M_Local_to_Global[inde,indi_l]=indi_g
            indi_g=indi_g+1

        indi_g=indi_g-1 #To start, in the next element, with the last DoF of the current element


    print(M_Local_to_Global)


    #First node
    DoF=node_Global_to_Local_class()
    DoF.N_el_containing_node = 1

    DoF.vec_el        = np.append(DoF.vec_el,0)
    DoF.vec_indi_l    = np.append(DoF.vec_indi_l,0)
    v_Global_to_Local = np.append(v_Global_to_Local,DoF) 



    # print(DoF.N_el_containing_node)
    # print(DoF.vec_el)
    # print(DoF.vec_indi_l)
    print()



    inde=0 #counter on the elements
    indi_l=1 #counter on the local index #NB: starting from 1 because starting from local DoFs

    for indi_g in range(1,N_global_nodes_v-1): #Loop on the internal nodes only

        DoF = node_Global_to_Local_class()

        DoF.N_el_containing_node = 1
        DoF.vec_el               = np.append(DoF.vec_el,inde)
        DoF.vec_indi_l           = np.append(DoF.vec_indi_l,indi_l)

        if indi_l<N_local_nodes_v-1: #Only one element containing the DoF
            v_Global_to_Local        = np.append(v_Global_to_Local,DoF) 
        else:
            DoF.N_el_containing_node = DoF.N_el_containing_node+1
            DoF.vec_el               = np.append(DoF.vec_el,inde+1)
            DoF.vec_indi_l           = np.append(DoF.vec_indi_l,0)
            v_Global_to_Local        = np.append(v_Global_to_Local,DoF) 

        indi_l=indi_l+1

        if indi_l==N_local_nodes_v:
            indi_l=1        
            inde=inde+1        
            


    #Last node
    DoF=node_Global_to_Local_class()
    DoF.N_el_containing_node = 1
    DoF.vec_el        = np.append(DoF.vec_el,N_el-1)
    DoF.vec_indi_l    = np.append(DoF.vec_indi_l,N_local_nodes_v-1)
    v_Global_to_Local = np.append(v_Global_to_Local,DoF) 
                

    if DATA.periodic==True:
        v_Global_to_Local[0].N_el_containing_node=2
        v_Global_to_Local[0].vec_el=np.array([N_el-1,0])
        v_Global_to_Local[0].vec_indi_l=np.array([N_local_nodes_v-1,0])


        v_Global_to_Local[N_global_nodes_v-1].N_el_containing_node=2
        v_Global_to_Local[N_global_nodes_v-1].vec_el=np.array([N_el-1,0])
        v_Global_to_Local[N_global_nodes_v-1].vec_indi_l=np.array([N_local_nodes_v-1,0])


    # for indi_g in range(N_global_nodes_v):
    #     print("DoF",indi_g)
    #     print(v_Global_to_Local[indi_g].N_el_containing_node, "elements contain it")
    #     print("Elements", v_Global_to_Local[indi_g].vec_el)
    #     print("Local indices", v_Global_to_Local[indi_g].vec_indi_l)
    #     print()
    # quit()




    # M_faces[indf,:]
    # M_faces[indf,0] -> left  element
    # M_faces[indf,1] -> right element
    M_faces=np.zeros((N_el+1,2),dtype=int)

    #Content is -1 if no element is present
    M_faces[0,0]    = -1 
    M_faces[0,1]    = 0 
    
    M_faces[N_el,1] = -1
    M_faces[N_el,0] = N_el-1

    for indf in range(1,N_el):
        M_faces[indf,0] = indf-1
        M_faces[indf,1] = indf

    if DATA.periodic==True:
        M_faces[0,0]    = N_el-1 
        M_faces[N_el,1] = 0

        #Elimination of the last face
        M_faces=M_faces[0:N_el,:]


    print(M_faces)


    # #-----------------------------------------------
    # print("Inside build_mesh")
    # print(local_nodes_H)
    # print(local_nodes_v)
    # print(degree_H)
    # print(degree_v)
    # print(N_global_nodes_v)
    # print(x_H)
    # print(x_v)
    # print(M_Local_to_Global)
    # # quit()
    # #-----------------------------------------------
    # for indi_g in range(N_global_nodes_v):
    #     print("DoF",indi_g)
    #     print("Contained by",v_Global_to_Local[indi_g].N_el_containing_node, "elements")
    #     print("...and these are",v_Global_to_Local[indi_g].vec_el)
    #     print("...and the local DoF in these is",v_Global_to_Local[indi_g].vec_indi_l)
    #     print()
    # quit()
    # #-----------------------------------------------
    # # print(M_faces)
    # #-----------------------------------------------

    return x_H, x_v, M_Local_to_Global, v_Global_to_Local, N_global_nodes_v, M_faces
          #x_H, x_v, M_Local_to_Global, v_Global_to_Local, N_global_nodes_v, M_faces
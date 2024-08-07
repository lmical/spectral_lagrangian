import numpy as np
import matplotlib.pyplot as plt
import os

#==============================================================
# INPUT
#==============================================================
test="Transcritical_Smooth"      #Supercritical_Smooth, Subcritical_Smooth, Transcritical_Smooth
perturbation=0
order_space=4
time_scheme="SSPRK4"                           #Time scheme #"Euler" "DeC" "SSPRK4"
LaxFriedrichs="Disabled"       #"Disabled" #"Active" #"ShockDetector_divV" (activated in troubled cells and neighbours) #"ShockDetector_divV_tn" (Same but detection only at time t_n)
K_limiter_divV=0.1                          #0.1 #8  #Important only if ShockDetector_divV_tn or ShockDetector_divV 
N_limited_neighbours=2                      #1   #2  #Important only if ShockDetector_divV_tn or ShockDetector_divV     
jump_CIP_in_v="j0"                          #Keep "j0", not needed, but in case there is jc
jump_eta_in_x=False                         #Keep False, not HO #Stopping term
jump_eta_in_H=False                         #Keep False, not HO #ALE-like
CFL=0.5
#==============================================================
#
#
#
#==============================================================
folder="./New_Results/"+test #Results_Conservative_Formulation #Results_Non_Conservative_Formulation #Results #New_Results
degree_H=order_space-1
degree_v=order_space
local_DoFs_H=degree_H+1
local_DoFs_v=degree_v+1
#==============================================================


if os.path.isdir(folder):  #CONDITION: Is it a folder? If yes go on
    count=0
    errorfiles=[]
    if LaxFriedrichs=="Active" or LaxFriedrichs=="Disabled":
        fileword="errors_pert"+str(perturbation)+"_"+"P"+str(order_space-1)+"P"+str(order_space)+"_"+time_scheme+"_LxF"+str(LaxFriedrichs)+"_"+jump_CIP_in_v+"_jeta"+str(jump_eta_in_x)+"_CFL"+str(CFL)
    elif LaxFriedrichs=="ShockDetector_divV" or LaxFriedrichs=="ShockDetector_divV_tn":
        fileword="errors_pert"+str(perturbation)+"_"+"P"+str(order_space-1)+"P"+str(order_space)+"_"+time_scheme+"_LxF"+str(LaxFriedrichs)+"_K"+str(K_limiter_divV)+"_NLimitedNeighbours"+str(N_limited_neighbours)+"_"+jump_CIP_in_v+"_jeta"+str(jump_eta_in_x)+"_CFL"+str(CFL)
    for file in os.listdir(folder): #CONDITION: Is there more than 1 error files?
        if file.startswith(fileword):
            count=count+1
            errorfiles.append(file)
    if count>1: #If yes go on, MAKE THE CONVERGENCE ANALYSIS
        errorfiles.sort() #Order the files
        # print(errorfiles)
        #Check where you are
        print("You are in the folder"+folder+" where you have "+str(count)+" solution files")
        #Order of convergence assuming half mesh size in subsequent meshes
        n_el_vec = np.zeros(len(errorfiles))
        errors_vec_L1_v  = np.zeros(len(errorfiles))
        errors_vec_L1_H  = np.zeros(len(errorfiles))
        errors_vec_L1_q  = np.zeros(len(errorfiles))
        errors_vec_2_v   = np.zeros(len(errorfiles))
        errors_vec_2_H   = np.zeros(len(errorfiles))


        for indr in range(count): #Index of the refinement

            namefile=folder+"/"+errorfiles[indr]

            lines=[]
            with open(namefile, 'r') as a:
                for idx, line in enumerate(a.readlines()):
                    lines.append(line)


            #Loop over the lines that we have read
            for idx, line in enumerate(lines): #rmk: 0-based numeration
                #Only first line to be read
                if idx!=1: 
                    pass
                else: 
                    data=line.split(",")
                    n_el_vec[indr]        = int(data[0])
                    errors_vec_L1_v[indr] = float(data[1])
                    errors_vec_L1_H[indr] = float(data[2])
                    errors_vec_L1_q[indr] = float(data[3])
                    errors_vec_2_v[indr]  = float(data[4])
                    errors_vec_2_H[indr]  = float(data[5])





        # print("Meshes elements",n_el_vec)
        # print("Errors L^1 v",errors_vec_L1_v)
        # print("Errors L^1 H",errors_vec_L1_H)
        # print("Errors L^1 q",errors_vec_L1_q)
        # print("Errors 2 v",errors_vec_2_v)
        # print("Errors 2 H",errors_vec_2_H)


        errors = np.zeros((len(errorfiles),5))
        errors[:,0]=errors_vec_L1_v.copy()
        errors[:,1]=errors_vec_L1_H.copy()
        errors[:,2]=errors_vec_L1_q.copy()
        errors[:,3]=errors_vec_2_v.copy()
        errors[:,4]=errors_vec_2_H.copy()

        #Opening the file to write the rates of convergence
        fid = open(folder+"/exact_convergence_"+fileword+".tex",'w')
        print("    N          L^1(v)     Order v          L^1(H)     Order H          L^1(q)     Order q         ||v||_2     Order v         ||H||_2     Order H")
        fid.write("  N   & L^1(v)  &  Order v & L^1(H)  &  Order H   & L^1(q)  &  Order q & ||v||_2  &  Order v   & ||H||_2  &  Order H\\ \n")  
        for indi in range(len(n_el_vec)): #Process the files
            if indi>0:
                order_space = np.array([-np.log(errors[indi,j]/errors[indi-1,j])/(np.log(n_el_vec[indi]/n_el_vec[indi-1])) for j in range(5)])
            else:
                order_space = np.zeros(5)
            print(format(int(n_el_vec[indi]), '5d'),"     ",format(errors[indi,0], '.3e'),"     ", format(order_space[0], '.3f'),"     ",format(errors[indi,1], '.3e'),"     ", format(order_space[1], '.3f'),"     ",format(errors[indi,2], '.3e'),"     ", format(order_space[2], '.3f'),"     ",format(errors[indi,3], '.3e'),"     ", format(order_space[3], '.3f'),"     ",format(errors[indi,4], '.3e'),"     ", format(order_space[4], '.3f'))
            fid.write(format(int(n_el_vec[indi]), '5d')+"  &   "+format(errors[indi,0], '.3e')+"  &  "+format(order_space[0], '.3f')+"  &  "+format(errors[indi,1], '.3e')+" & "+format(order_space[1], '.3f')+"  &  "+format(errors[indi,2], '.3e')+" & "+format(order_space[2], '.3f')+"  &  "+format(errors[indi,3], '.3e')+" & "+format(order_space[3], '.3f')+"  &  "+format(errors[indi,4], '.3e')+" & "+format(order_space[4], '.3f')+" \\ \n")  
        fid.close()



        #Plot
        fig=plt.figure()
        plt.loglog(n_el_vec,errors[:,0],"-*",linewidth=1.5,label='L^1(v)')
        plt.grid()
        plt.loglog(n_el_vec,errors[:,1],"-+",linewidth=1.5,label='L^1(H)')
        plt.loglog(n_el_vec,errors[:,2],"-+",linewidth=1.5,label='L^1(q)')
        plt.loglog(n_el_vec,errors[:,3],"-+",linewidth=1.5,label='||v||_2')
        plt.loglog(n_el_vec,errors[:,4],"-+",linewidth=1.5,label='||H||_2')

        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-1)*n_el_vec**(-1),"--",linewidth=1,label="1st order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-2)*n_el_vec**(-2),"--",linewidth=1,label="2nd order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-3)*n_el_vec**(-3),"--",linewidth=1,label="3rd order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-4)*n_el_vec**(-4),"--",linewidth=1,label="4th order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-5)*n_el_vec**(-5),"--",linewidth=1,label="5th order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-6)*n_el_vec**(-6),"--",linewidth=1,label="6th order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-7)*n_el_vec**(-7),"--",linewidth=1,label="7th order")
        plt.loglog(n_el_vec,errors[0,0]/n_el_vec[0]**(-8)*n_el_vec**(-8),"--",linewidth=1,label="8th order")
        plt.legend(loc='lower left')
        params = {'mathtext.default': 'regular' }   
        plt.xlabel("$N_{elements}$")
        plt.ylabel("$L^1$ error")
        plt.savefig(folder+"/exact_convergence_"+fileword+".pdf", format="pdf", bbox_inches="tight")
        # plt.savefig(folder+"/exact_convergence_"+fileword+".png",dpi=600)
        # plt.show()



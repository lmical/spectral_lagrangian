import numpy as np
import matplotlib.pyplot as plt
import os

#==============================================================
# INPUT
#==============================================================
test="Smooth_periodic"      #Smooth_periodic, Supercritical_Smooth, Subcritical_Smooth, Transcritical_Smooth
perturbation=0
order_space=4
time_scheme="DeC"                       #Time scheme #"Euler" "DeC" "SSPRK4"
LaxFriedrichs="Disabled"   #"Disabled" #"Active" #"ShockDetector_divV" (activated in troubled cells and neighbours) #"ShockDetector_divV_tn" (Same but detection only at time t_n)
K_limiter_divV=6.5                      #0.1 #8  #Important only if ShockDetector_divV_tn or ShockDetector_divV 
N_limited_neighbours=2                  #1   #2  #Important only if ShockDetector_divV_tn or ShockDetector_divV     
jump_CIP_in_v="j0"                      #Keep "j0", not needed, but in case there is jc
jump_eta_in_x=False                     #Keep False, not HO #Stopping term
jump_eta_in_H=False                     #Keep False, not HO #ALE-like
CFL=0.5
# if test=="Constant_Slope_Smooth" or test=="No_Slope_Smooth":
#     perturbation=1
#==============================================================
#
#
#
#==============================================================
folder="./New_Results/"+test #Results_Conservative_Formulation #Results_Non_Conservative_Formulation #Results
degree_H=order_space-1
degree_v=order_space
local_DoFs_H=degree_H+1
local_DoFs_v=degree_v+1
#==============================================================


if os.path.isdir(folder):  #CONDITION: Is it a folder? If yes go on
    count=0
    errorfiles=[]
    if LaxFriedrichs=="Active" or LaxFriedrichs=="Disabled":
        fileword="values_pert"+str(perturbation)+"_"+"P"+str(order_space-1)+"P"+str(order_space)+"_"+time_scheme+"_LxF"+str(LaxFriedrichs)+"_"+jump_CIP_in_v+"_jeta"+str(jump_eta_in_x)+"_CFL"+str(CFL)
    elif LaxFriedrichs=="ShockDetector_divV" or LaxFriedrichs=="ShockDetector_divV_tn":
        fileword="values_pert"+str(perturbation)+"_"+"P"+str(order_space-1)+"P"+str(order_space)+"_"+time_scheme+"_LxF"+str(LaxFriedrichs)+"_K"+str(K_limiter_divV)+"_NLimitedNeighbours"+str(N_limited_neighbours)+"_"+jump_CIP_in_v+"_jeta"+str(jump_eta_in_x)+"_CFL"+str(CFL)
    for file in os.listdir(folder): #CONDITION: Is there more than 1 error files?
        if file.startswith(fileword):
            count=count+1
            errorfiles.append(file)
    if count>1: #If yes go on, MAKE THE CONVERGENCE ANALYSIS
        errorfiles.sort() #Order the files
        # print(errorfiles)
        #Check where you are
        print("You are in the folder"+folder+" where you have "+str(count)+" solution files")
        #Experimental order of convergence assuming half mesh size in subsequent meshes
        errors_vec_x  = np.zeros(len(errorfiles)-1)
        errors_vec_v  = np.zeros(len(errorfiles)-1)
        errors_vec_H  = np.zeros(len(errorfiles)-1)
        errors_vec_q  = np.zeros(len(errorfiles)-1)
        errors_vec_eta  = np.zeros(len(errorfiles)-1)

        n_el_vec = np.zeros(len(errorfiles))


        for indr in range(1,count): #Index of the refinement

            namecoarse=folder+"/"+errorfiles[indr-1]

            lines_coarse=[]
            with open(namecoarse, 'r') as a:
                for idx, line in enumerate(a.readlines()):
                    lines_coarse.append(line)


            x_v_coarse=np.array([]); v_coarse=np.array([]); H_coarse=np.array([]); q_coarse=np.array([]); eta_coarse=np.array([]); 
            #Loop over the lines that we have read
            for idx, line in enumerate(lines_coarse): #rmk: 0-based numeration
                if idx==0: #first line we skip, it is the header
                    pass
                else:
                    data=line.split(",")
                    indi=int(data[0])
                    x_v_coarse=np.append(x_v_coarse,float(data[1]))
                    v_coarse=np.append(v_coarse,float(data[2]))
                    H_coarse=np.append(H_coarse,float(data[3]))
                    q_coarse=np.append(q_coarse,float(data[4]))
                    eta_coarse=np.append(eta_coarse,float(data[5]))


            n_el_vec[indr-1]=int(indi/degree_v)


            namerefined=folder+"/"+errorfiles[indr]
            lines_refined=[]
            with open(namerefined, 'r') as a:
                for idx, line in enumerate(a.readlines()):
                    lines_refined.append(line)

            x_v_refined=np.array([]); v_refined=np.array([]); H_refined=np.array([]); q_refined=np.array([]); eta_refined=np.array([]); 
            #Loop over the lines that we have read
            for idx, line in enumerate(lines_refined): #rmk: 0-based numeration
                if idx==0: #first line we skip, it is the header
                    pass
                else:
                    data=line.split(",")
                    indi=int(data[0])
                    x_v_refined=np.append(x_v_refined,float(data[1]))
                    v_refined=np.append(v_refined,float(data[2]))
                    H_refined=np.append(H_refined,float(data[3]))
                    q_refined=np.append(q_refined,float(data[4]))
                    eta_refined=np.append(eta_refined,float(data[5]))

            n_el_vec[indr]=int(indi/degree_v)

            #Isolate correct ones
            x_v_coarse  = x_v_coarse[::degree_v]
            v_coarse    = v_coarse[::degree_v]
            H_coarse    = H_coarse[::degree_v]
            q_coarse    = q_coarse[::degree_v]
            eta_coarse  = eta_coarse[::degree_v]


            x_v_refined = x_v_refined[::degree_v]
            v_refined   = v_refined[::degree_v]
            H_refined   = H_refined[::degree_v]
            q_refined   = q_refined[::degree_v]
            eta_refined = eta_refined[::degree_v]

            x_v_refined = x_v_refined[::2]
            v_refined   = v_refined[::2]
            H_refined   = H_refined[::2]
            q_refined   = q_refined[::2]
            eta_refined = eta_refined[::2]

            errorx=np.linalg.norm(x_v_refined-x_v_coarse)/np.sqrt(len(x_v_coarse))
            errorv=np.linalg.norm(v_refined-v_coarse)/np.sqrt(len(v_coarse))
            errorH=np.linalg.norm(H_refined-H_coarse)/np.sqrt(len(H_coarse))
            errorq=np.linalg.norm(q_refined-q_coarse)/np.sqrt(len(q_coarse))
            erroreta=np.linalg.norm(eta_refined-eta_coarse)/np.sqrt(len(eta_coarse))
            

            errors_vec_x[indr-1]=errorx
            errors_vec_v[indr-1]=errorv
            errors_vec_H[indr-1]=errorH
            errors_vec_q[indr-1]=errorq
            errors_vec_eta[indr-1]=erroreta

        # print("Meshes elements",n_el_vec)
        # print("Errors x",errors_vec_x)
        # print("Errors v",errors_vec_v)
        # print("Errors H",errors_vec_H)
        # print("Errors q",errors_vec_q)
        # print("Errors eta",errors_vec_eta)

        errors = np.zeros((len(errorfiles)-1,5))
        errors[:,0]=errors_vec_x.copy()
        errors[:,1]=errors_vec_v.copy()
        errors[:,2]=errors_vec_H.copy()
        errors[:,3]=errors_vec_q.copy()
        errors[:,4]=errors_vec_eta.copy()

        #Opening the file to write the rates of convergence
        fid = open(folder+"/experimental_convergence_"+fileword+".tex",'w')
        print("    N         Error x     Order x         Error v     Order v         Error H     Order H         Error q     Order q      Error \eta  Order \eta")
        fid.write("  N   & Error x  &  Order x & Error v  &  Order v   & Error H  &  Order H & Error q  &  Order q   & Error \eta  &  Order \eta\\ \n")  
        for indi in range(len(errors)): #Process the files
            if indi>0:
                order_space = np.array([-np.log(errors[indi,j]/errors[indi-1,j])/(np.log(2)) for j in range(5)])
            else:
                order_space = np.zeros(5)
            print(format(int(n_el_vec[indi+1]), '5d'),"     ",format(errors[indi,0], '.3e'),"     ", format(order_space[0], '.3f'),"     ",format(errors[indi,1], '.3e'),"     ", format(order_space[1], '.3f'),"     ",format(errors[indi,2], '.3e'),"     ", format(order_space[2], '.3f'),"     ",format(errors[indi,3], '.3e'),"     ", format(order_space[3], '.3f'),"     ",format(errors[indi,4], '.3e'),"     ", format(order_space[4], '.3f'))
            fid.write(format(int(n_el_vec[indi+1]), '5d')+"  &   "+format(errors[indi,0], '.3e')+"  &  "+format(order_space[0], '.3f')+"  &  "+format(errors[indi,1], '.3e')+" & "+format(order_space[1], '.3f')+"  &  "+format(errors[indi,2], '.3e')+" & "+format(order_space[2], '.3f')+"  &  "+format(errors[indi,3], '.3e')+" & "+format(order_space[3], '.3f')+"  &  "+format(errors[indi,4], '.3e')+" & "+format(order_space[4], '.3f')+" \\ \n")  
        fid.close()

        #Plot
        fig=plt.figure()
        plt.loglog(n_el_vec[1:],errors[:,0],"-*",linewidth=1.5,label='Error x')
        plt.grid()
        plt.loglog(n_el_vec[1:],errors[:,1],"-+",linewidth=1.5,label='Error v')
        plt.loglog(n_el_vec[1:],errors[:,2],"-+",linewidth=1.5,label='Error H')
        plt.loglog(n_el_vec[1:],errors[:,3],"-+",linewidth=1.5,label='Error q')
        plt.loglog(n_el_vec[1:],errors[:,4],"-+",linewidth=1.5,label='Error \eta')

        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-1)*n_el_vec[1:]**(-1),"--",linewidth=1,label="1st order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-2)*n_el_vec[1:]**(-2),"--",linewidth=1,label="2nd order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-3)*n_el_vec[1:]**(-3),"--",linewidth=1,label="3rd order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-4)*n_el_vec[1:]**(-4),"--",linewidth=1,label="4th order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-5)*n_el_vec[1:]**(-5),"--",linewidth=1,label="5th order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-6)*n_el_vec[1:]**(-6),"--",linewidth=1,label="6th order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-7)*n_el_vec[1:]**(-7),"--",linewidth=1,label="7th order")
        plt.loglog(n_el_vec[1:],errors[0,0]/n_el_vec[1]**(-8)*n_el_vec[1:]**(-8),"--",linewidth=1,label="8th order")
        plt.legend(loc='lower left')
        params = {'mathtext.default': 'regular' }   
        plt.xlabel("$N_{elements}$")
        plt.ylabel("$L^1$ error")
        plt.savefig(folder+"/experimental_convergence_"+fileword+".pdf", format="pdf", bbox_inches="tight")
        # plt.savefig(folder+"/experimental_convergence_"+fileword+".png",dpi=600)
        # plt.show()



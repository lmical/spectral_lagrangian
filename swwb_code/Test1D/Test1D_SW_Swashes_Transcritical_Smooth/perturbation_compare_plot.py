nametest="TranscriticalSmoothSwashesDiscreteEquilibrium"
#setting=1
#scheme=[4,14,24]
#if setting==1:
#    jump=[1,2,3,4,5,6,7]
#else:
#    jump=[11,12,13,14,15,16,17]
#bf=["P1","PGL1","B2","P2","P3","B3","PGL2","PGL3","PGL4","B4","P4"]

teststocompare=[] #setting=1, scheme, jump, bf
teststocompare.append([1,14,1,"PGL2"]) 
teststocompare.append([1,14,2,"PGL2"]) 

import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

#test="../AnotherFolder";
test="./" #Actual folder, not used but in case you want to run it from a specific folder you can use it to implement it


params = {'mathtext.default': 'regular' } #parameters plot

#Needed to plot the free surface
def bathymetry(x):
    x0=10.
    r=5.
    b=0.
    if (x > x0-r) and (x < x0+r):
        b=0.2*np.exp(1. - 1./(1.-((x-x0)/r)**2.))    
    else:
        b=0.
    return b




for indt, test in enumerate(teststocompare): #Loop on the schemes
    setting=test[0]
    ischeme=test[1]
    ijump=test[2]
    ibf=test[3]
    schemefolder='setting'+str(setting)+'/scheme'+str(ischeme)
    jumpfolder=schemefolder+'/jump'+str(ijump)
    bffolder=jumpfolder+"/"+ibf #Folder to visit if possible


    if os.path.isdir(bffolder):  #CONDITION: Is it a folder? If yes go on

        folder=bffolder
        namefile=folder+'/'+nametest+'.dat'
        if os.path.isfile(namefile): #If the file exists, read it
            print("You are in the folder "+bffolder+" dealing with "+nametest)
            lines=[]
            with open(namefile, 'r') as a:
                for idx, line in enumerate(a.readlines()):
                    lines.append(line)


            x=np.array([])
            hf=np.array([]); huf=np.array([]); etaf=np.array([]); vf=np.array([]) #final
            h0=np.array([]); hu0=np.array([]); eta0=np.array([]); v0=np.array([]) #initial

            #Loop over the lines that we have read
            for idx, line in enumerate(lines): #rmk: 0-based numeration
                if idx==0: #first line we skip, it is the header
                    pass
                else:
                    data=line.split()
                    indi=int(data[0])
                    x=np.append(x,float(data[1]))
                    hf=np.append(hf,float(data[2]))
                    huf=np.append(huf,float(data[3]))
                    h0=np.append(h0,float(data[4]))
                    hu0=np.append(hu0,float(data[5]))
                    etaf=np.append(etaf, float(data[2]) + bathymetry(float(data[1])) )
                    eta0=np.append(eta0, float(data[4]) + bathymetry(float(data[1])) ) 
                    vf=np.append(vf, float(data[3])/float(data[2]))
                    v0=np.append(v0, float(data[5])/float(data[4]))




            fig, axs = plt.subplots(2,2) #Array of subplots
            fig.suptitle("values_scheme"+str(ischeme)+"_jump"+str(ijump))

            #H
            axs[0,0].plot(x,hf, marker="*")
            axs[0,0].set_title("H")
            axs[0,0].set_xlabel("x")
            axs[0,0].grid()

            #v
            axs[0,1].plot(x,vf)
            axs[0,1].set_title("v")
            axs[0,1].set_xlabel("x")
            axs[0,1].grid()

            #eta
            axs[1,0].plot(x,etaf, marker="*")
            axs[1,0].set_title("\eta")
            axs[1,0].set_xlabel("x")
            axs[1,0].grid()

            #q
            axs[1,1].plot(x,hf*vf)
            axs[1,1].set_title("q")
            axs[1,1].set_xlabel("x")
            axs[1,1].grid()

            fig.tight_layout()

            plt.savefig("values_cons_scheme"+str(ischeme)+"_jump"+str(ijump)+".pdf", format="pdf", bbox_inches="tight")

            plt.show()

            errorH=h0-hf
            errorv=v0-vf
            errorq=hu0-huf

            fig, axs = plt.subplots(1,3) #Array of subplots
            fig.suptitle("errors_scheme"+str(ischeme)+"_jump"+str(ijump))

            #H
            axs[0].plot(x,errorH, marker="*")
            axs[0].set_title("H")
            axs[0].set_xlabel("x")
            axs[0].grid()

            #v
            axs[1].plot(x,errorv)
            axs[1].set_title("v")
            axs[1].set_xlabel("x")
            axs[1].grid()


            #q
            axs[2].plot(x,errorq)
            axs[2].set_title("q")
            axs[2].set_xlabel("x")
            axs[2].grid()


            fig.tight_layout()

            plt.savefig("errors_cons_scheme"+str(ischeme)+"_jump"+str(ijump)+".pdf", format="pdf", bbox_inches="tight")

            plt.show()

import numpy as np
import shutil, os
import re
from joblib import Parallel, delayed
import multiprocessing
import sys



test="Sod_Transcritical_Expansion"       #Sod_Transcritical_Expansion #Sod   #Smooth_periodic, Supercritical_Smooth, Subcritical_Smooth, Transcritical_Smooth
perturbation=0
order_space=3
time_scheme="SSPRK4"                     #Time scheme #"Euler" "DeC" "SSPRK4"
LaxFriedrichs="ShockDetector_divV_tn"    #"Disabled" #"Active" #"ShockDetector_divV" (activated in troubled cells and neighbours) #"ShockDetector_divV_tn" (Same but detection only at time t_n)
K_limiter_divV=0.1                       #0.1 #8   #Important only if ShockDetector_divV_tn or ShockDetector_divV
N_limited_neighbours=2                   #1   #2   #Important only if ShockDetector_divV_tn or ShockDetector_divV
jump_CIP_in_v="j0"                       #Keep "j0", not needed, but in case there is jc
jump_eta_in_x=False                      #Keep False, not HO #Stopping term
jump_eta_in_H=False                      #Keep False, not HO #ALE-like
CFL=0.5


if len(sys.argv)>1:
    test=sys.argv[1]
if len(sys.argv)>2:
    perturbation=int(sys.argv[2])
if len(sys.argv)>3:
    order_space=int(sys.argv[3])
if len(sys.argv)>4:
    time_scheme=sys.argv[4]
if len(sys.argv)>5:
    LaxFriedrichs=sys.argv[5]
if len(sys.argv)>6:
    K_limiter_divV=float(sys.argv[6])
if len(sys.argv)>7:
    N_limited_neighbours=int(sys.argv[7])
if len(sys.argv)>8:
    jump_CIP_in_v=sys.argv[8]
if len(sys.argv)>9:
    if sys.argv[9]=="True":
        jump_eta_in_x=True
    elif sys.argv[9]=="False":
        jump_eta_in_x=False
    else:
        print("Impossible to get jump_eta_in_x imput from keyboard")
        quit()
if len(sys.argv)>10:
    if sys.argv[10]=="True":
        jump_eta_in_H=True
    elif sys.argv[10]=="False":
        jump_eta_in_H=False
    else:
        print("Impossible to get jump_eta_in_H imput from keyboard")
        quit()
if len(sys.argv)>11:
    CFL=float(sys.argv[11])

instruction = ""
instruction +="python3 main_sw.py  "+test\
                                +" "+str(perturbation)\
                                +" "+str(order_space)\
                                +" "+time_scheme\
                                +" "+LaxFriedrichs\
                                +" "+str(K_limiter_divV)\
                                +" "+str(N_limited_neighbours)\
                                +" "+jump_CIP_in_v\
                                +" "+str(jump_eta_in_x)\
                                +" "+str(jump_eta_in_H)\
                                +" "+str(CFL)


if test in ["Sod","Sod_Transcritical_Expansion"]:
    elements=[100,200] #Shock    
elif test in ["Supercritical_Smooth","Subcritical_Smooth","Transcritical_Smooth"]:
    elements=[10,20,40,80,160,320] #Convergence analysis
else:
    print("STOP, careful with elements")
    quit()

for inde in elements:
    instructiontorun=instruction+" "+str(inde)
    print(instructiontorun)
    os.system(instructiontorun)


# python3 main_sw.py Smooth_periodic 0 2 DeC Disabled 0 0 j0 False False 0.5 20
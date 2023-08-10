IC_file="IC_supercritical_smooth_exact.py"

import numpy as np
import shutil, os
import re
from joblib import Parallel, delayed
import multiprocessing
import sys



#PARAMETERS

#setting=1 #my parameters

n_el=[128] #[16, 32, 64, 128, 256] #n_el vectors of the number of elements for the convergence analysis

scheme=14 #4,14,24
jumps=[1,2] #[1,2,3,4,5,6,7]

#parameters for the orders of the polynomials
#rmk: order of the scheme = order of the polynomials+1
first_iorder=2 
first_iiter=2
first_itheta1=0.05
first_itheta2=0.5
first_CFL=0.1

second_iorder=3 
second_iiter=3
second_itheta1=0.3
second_itheta2=0.2
second_CFL=0.1


third_iorder=4 
third_iiter=4
third_itheta1=0.15
third_itheta2=0.2
third_CFL=0.1


fourth_iorder=5 
fourth_iiter=5
fourth_itheta1=0.5 #0.15
fourth_itheta2=0.01 #0.2
fourth_CFL=0.05 #NOT 0.1


##############################################################################
print("Defining class of the specific tests")
class SpecificTest:
    name="" #name of the basis functions
    itype=0 #line 2 #type of the basis functions
    iorder=0 #line 3.1 #order 
    iiter=0 #line 3.2 #iterations
    itheta1=0 #line 5
    itheta2=0 #line 6
    CFL=0 #line 7

print("Defining vector of the specific tests")
TESTS=[] #Vector of tests #P1 #P2 #B2 #PGL2 #P3 #B3 #PGL3 #PGL4 #B4 #and P4 disabled because unstable

##########################################################################33
print("Filling the vector")


#0) FIRST ORDER POLYNOMIALS 
#P1 itype=1
#PGL1 itype=11

#P1
P1=SpecificTest()
P1.name="P1" 
P1.itype=1
P1.iorder=first_iorder
P1.iiter=first_iiter
P1.itheta1=first_itheta1
P1.itheta2=first_itheta2
P1.CFL=first_CFL

TESTS.append(P1)

####
#print(TESTS[0].name)
#print(TESTS[0].itype)
#print(TESTS[0].iorder, TESTS[0].iiter)
#print(TESTS[0].itheta1)
#print(TESTS[0].itheta2)
#print(TESTS[0].CFL)
#sys.exit()
#####

#PGL1
PGL1=SpecificTest()
PGL1.name="PGL1" 
PGL1.itype=11
PGL1.iorder=first_iorder
PGL1.iiter=first_iiter
PGL1.itheta1=first_itheta1
PGL1.itheta2=first_itheta2
PGL1.CFL=first_CFL

TESTS.append(PGL1)

####
#print(TESTS[1].name)
#print(TESTS[1].itype)
#print(TESTS[1].iorder, TESTS[1].iiter)
#print(TESTS[1].itheta1)
#print(TESTS[1].itheta2)
#print(TESTS[1].CFL)
#sys.exit()
#####

#1) SECOND ORDER POLYNOMIALS 
#B2 itype=2
#P2 itype=3
#PGL2 itype=12


#B2
B2=SpecificTest()
B2.name="B2" 
B2.itype=2
B2.iorder=second_iorder
B2.iiter=second_iiter
B2.itheta1=second_itheta1
B2.itheta2=second_itheta2
B2.CFL=second_CFL

TESTS.append(B2)

####
#print(TESTS[2].name)
#print(TESTS[2].itype)
#print(TESTS[2].iorder, TESTS[2].iiter)
#print(TESTS[2].itheta1)
#print(TESTS[2].itheta2)
#print(TESTS[2].CFL)
#sys.exit()
####

#P2
P2=SpecificTest()
P2.name="P2" 
P2.itype=3
P2.iorder=second_iorder
P2.iiter=second_iiter
P2.itheta1=second_itheta1
P2.itheta2=second_itheta2
P2.CFL=second_CFL

TESTS.append(P2)

####
#print(TESTS[3].name)
#print(TESTS[3].itype)
#print(TESTS[3].iorder, TESTS[3].iiter)
#print(TESTS[3].itheta1)
#print(TESTS[3].itheta2)
#print(TESTS[3].CFL)
#sys.exit()
####

#PGL2
PGL2=SpecificTest()
PGL2.name="PGL2" 
PGL2.itype=12
PGL2.iorder=second_iorder
PGL2.iiter=second_iiter
PGL2.itheta1=second_itheta1
PGL2.itheta2=second_itheta2
PGL2.CFL=second_CFL

TESTS.append(PGL2)

####
#print(TESTS[4].name)
#print(TESTS[4].itype)
#print(TESTS[4].iorder, TESTS[4].iiter)
#print(TESTS[4].itheta1)
#print(TESTS[4].itheta2)
#print(TESTS[4].CFL)
#sys.exit()
####

#2) THIRD ORDER POLYNOMIALS 
#P3 itype=4
#B3 itype=5
#PGL3 itype=13


#P3
P3=SpecificTest()
P3.name="P3" 
P3.itype=4
P3.iorder=third_iorder
P3.iiter=third_iiter
P3.itheta1=third_itheta1
P3.itheta2=third_itheta2
P3.CFL=third_CFL

TESTS.append(P3)

####
#print(TESTS[5].name)
#print(TESTS[5].itype)
#print(TESTS[5].iorder, TESTS[5].iiter)
#print(TESTS[5].itheta1)
#print(TESTS[5].itheta2)
#print(TESTS[5].CFL)
#sys.exit()
####

#B3
B3=SpecificTest()
B3.name="B3" 
B3.itype=5
B3.iorder=third_iorder
B3.iiter=third_iiter
B3.itheta1=third_itheta1
B3.itheta2=third_itheta2
B3.CFL=third_CFL

TESTS.append(B3)

####
#print(TESTS[6].name)
#print(TESTS[6].itype)
#print(TESTS[6].iorder, TESTS[6].iiter)
#print(TESTS[6].itheta1)
#print(TESTS[6].itheta2)
#print(TESTS[6].CFL)
#sys.exit()
####

#PGL3
PGL3=SpecificTest()
PGL3.name="PGL3" 
PGL3.itype=13
PGL3.iorder=third_iorder
PGL3.iiter=third_iiter
PGL3.itheta1=third_itheta1
PGL3.itheta2=third_itheta2
PGL3.CFL=third_CFL

TESTS.append(PGL3)

####
#print(TESTS[7].name)
#print(TESTS[7].itype)
#print(TESTS[7].iorder, TESTS[7].iiter)
#print(TESTS[7].itheta1)
#print(TESTS[7].itheta2)
#print(TESTS[7].CFL)
#sys.exit()
####


#3) FOURTH ORDER POLYNOMIALS 
#PGL4 itype=14
#B4 itype=6
#P4 itype=7 #not stable and so not appended

#PGL4
PGL4=SpecificTest()
PGL4.name="PGL4" 
PGL4.itype=14
PGL4.iorder=fourth_iorder
PGL4.iiter=fourth_iiter
PGL4.itheta1=fourth_itheta1
PGL4.itheta2=fourth_itheta2
PGL4.CFL=fourth_CFL

TESTS.append(PGL4)

####
#print(TESTS[8].name)
#print(TESTS[8].itype)
#print(TESTS[8].iorder, TESTS[8].iiter)
#print(TESTS[8].itheta1)
#print(TESTS[8].itheta2)
#print(TESTS[8].CFL)
#sys.exit()
####

#B4
B4=SpecificTest()
B4.name="B4" 
B4.itype=6
B4.iorder=fourth_iorder
B4.iiter=fourth_iiter
B4.itheta1=fourth_itheta1
B4.itheta2=fourth_itheta2
B4.CFL=fourth_CFL

TESTS.append(B4)

####
#print(TESTS[9].name)
#print(TESTS[9].itype)
#print(TESTS[9].iorder, TESTS[9].iiter)
#print(TESTS[9].itheta1)
#print(TESTS[9].itheta2)
#print(TESTS[9].CFL)
#sys.exit()
####

#P4 #<-not stable and so not appended
P4=SpecificTest()
P4.name="P4" 
P4.itype=7
P4.iorder=fourth_iorder
P4.iiter=fourth_iiter
P4.itheta1=fourth_itheta1
P4.itheta2=fourth_itheta2
P4.CFL=fourth_CFL

#TESTS.append(P4) #<-not stable and so not appended

####
#print(TESTS[10].name)
#print(TESTS[10].itype)
#print(TESTS[10].iorder, TESTS[10].iiter)
#print(TESTS[10].itheta1)
#print(TESTS[10].itheta2)
#print(TESTS[10].CFL)
#sys.exit()
####


################################################################
print("Checking the vector")

for indi in range(len(TESTS)):
    print(TESTS[indi].name)
    print(TESTS[indi].itype)
    print(TESTS[indi].iorder, TESTS[indi].iiter)
    print(TESTS[indi].itheta1)        
    print(TESTS[indi].itheta2)  
    print(TESTS[indi].CFL)        
    print()

#sys.exit()    
################################################################

#---------------------------------------------------------------
#run_single_basis_function runs the test for the single basis function
#TEST element in the list TESTS
#n_el vectors of the number of elements for the convergence analysis 
#test_folder NOT USED but in case you need to refer to specific folders you can implement using this
#---------------------------------------------------------------
def run_single_basis_function(TEST, n_el, test_folder):
    for n in n_el: #for each number of elements
        if ((TEST.name in ["B4","P4","PGL4"]) and (n>256)):
            continue
        print(TEST.name, n)              
        for ijump in jumps: 
            #CREATE THE FOLDER OF TEST IF NOT PRESENT AND COPY THE don1d THERE
            #IF IT ALREADY EXISTS ONLY COPY THE don1d THERE
            #COPY ALSO THE FILE TO GENERATE THE INITIAL CONDITION IC_file #<-----
            foldName = "setting1/scheme"+str(scheme)+"/jump"+str(ijump)+"/"+TEST.name
            instruction=""
            instruction +="mkdir -p "+foldName+" \n" #creation of the folder
            instruction +="mkdir -p "+foldName+"/Data \n" #creation of DATA
            instruction +="cp Data/don1d "+foldName+"/Data/don1d \n" #copying DATA/don1d
            instruction +="cp "+IC_file+" "+foldName+"/"+IC_file+"\n" #copying IC_file #<-----
            #print(instruction)
            os.system(instruction)
            #MODIFY THE don1d WITH THE PARAMETERS
            modify_don(foldName+"/Data/don1d", TEST, ijump)
            #MODIFY THE IC_file WITH THE TYPE OF ELEMENTS AND THE number of elements n
            modify_IC_file(foldName+"/"+IC_file, TEST, n)
            #RUN THE IC_file TO GENERATE THE IC
            instruction = ""
            instruction +="cd "+foldName+" \n "  #changing folder into the specific test
            #RUN THE IC_file
            instruction +=  "python3 "+IC_file+ " \n"  #running the IC_file to generate the IC
            os.system(instruction)            
            #FINALLY RUN THE TEST
            instruction = ""
            instruction +="cd "+foldName+" \n "  #changing folder into the specific test
            #RUN FROM THERE
            instruction +=  "../../../../../bin1D/main_dec.out "+str(n)+ " \n"  #running the test with n elements
            os.system(instruction)


#---------------------------------------------------------------    
#modify_don modifies the don1d in DATA with the parameters in TEST
#fname name of the file
#TEST specific test
#---------------------------------------------------------------    
def modify_don(fname, TEST, ijump):
    #first we read the file copying the lines in don
    don=[]
    with open(fname, 'r') as a:
        for idx, line in enumerate(a.readlines()):
            don.append(line)

    #then we write it
    f=open(fname,'w')
    #loop over the lines that we have copied
    for idx, line in enumerate(don): #rmk: 0-based numeration
        if idx==1: #itype
            f.write(str(TEST.itype)+'  itype 1: P1, 2: B2, 3: P2, 4:P3, 5: B3, 6: B4, 7: P4, 11,12,13,14: PGL\n')
        elif idx==2: #iorder #iiter				
            f.write(str(TEST.iorder)+' '+str(TEST.iiter)+' .TRUE.   #order #iterations for DEC #staggered \n')
        elif idx==3: #scheme #jump #GF
            GF=0
            if scheme==24 or ijump in [6,16,7,17]:
                GF=1
            f.write(str(scheme)+' '+str(ijump)+' '+str(GF)+' scheme: 4=GAL, 14=GALWB, 24=GF, jump: 1 (Burman),2 (eta), 3(entropy), 4(residual), 5(residual spectral), 6(GF), 7(GF spectral), GF: 1=WB, 2=STANDARD\n')
        elif idx==4: #theta1				
            f.write(str(TEST.itheta1)+'   theta parameter in Burman stabilization first derivative term  \n')
        elif idx==5: #theta2				
            f.write(str(TEST.itheta2)+'   theta2 parameter in Burman stabilization second derivative term  \n')
        elif idx==6: #CFL				
            f.write(str(TEST.CFL)+'    cfl  \n')		
        elif idx==7: #Maximal number of iterations, increase for security				
            f.write('100000000    ktmax\n')		
        elif idx==9: #ifre increase for security				
            f.write('5000     ifre\n')		
        else:
            f.write(don[idx])
            
    f.close()

    return



#---------------------------------------------------------------
#modify_IC_file modifies the IC_file changing the type of the elements and the number of the elements
#fname name of the file
#TEST specific test
#---------------------------------------------------------------    
def modify_IC_file(fname, TEST, n):
    #first we read the file copying the lines in don
    don=[]
    with open(fname, 'r') as a:
        for idx, line in enumerate(a.readlines()):
            don.append(line)

    #then we write it
    f=open(fname,'w')
    #loop over the lines that we have copied
    for idx, line in enumerate(don): #rmk: 0-based numeration
        if idx==0: #itype
            f.write('itype="'+str(TEST.name)+'" #P1, B2, P2, P3, B3, B4, P4, PGL1, PGL2, PGL3, PGL4\n')
        elif idx==1: #number of elements				
            f.write('n_el='+str(n)+'\n')
        else:
            f.write(don[idx])
            
    f.close()

    return



test_folder=""
#run_single_basis_function(TESTS[1], n_el, test_folder)

num_cores = multiprocessing.cpu_count()-1

#PARALLELIZED
#Parallel(n_jobs=num_cores)(delayed(run_single_basis_function)(TESTS[indi], n_el, test_folder) for indi in range(len(TESTS)))

#NOT PARALLELIZED
#for indi in range(len(TESTS)):
#    run_single_basis_function(TESTS[indi], n_el, test_folder)

run_single_basis_function(TESTS[4], n_el, test_folder)


sys.exit()


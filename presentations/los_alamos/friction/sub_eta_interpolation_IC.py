itype="PGL4" #type of the coarse mesh P1, B2, P2, P3, B3, B4, P4, PGL1, PGL2, PGL3, PGL4
n_el=128
file_refined="SubcriticalNotSmoothFrictionSwashesDiscreteEquilibrium" #to add .dat
itype_ref="P1" #type of the refined mesh
#print("Parameters of the test chosen")
#-----------------------------------------------------------------------------
#Since I want to interpolate eta I need the bathymetry #<--!NEW
def bathymetry(x):
    b=0.
    if((x > 8.) and (x < 12.)): 
        b=0.2-0.05*(x-10.)**2
    return b
#-----------------------------------------------------------------------------
#Length domain
L=25.
#print("Length of the domain", L)
#-----------------------------------------------------------------------------
import numpy as np 
from numpy import linalg as LA
from scipy import optimize as opt
from scipy.interpolate import lagrange
import matplotlib.pyplot as pl
import sys 
import itertools
#-----------------------------------------------------------------------------
#Construction of the coarse mesh where we will interpolate the solution that we have on the refined mesh
print("Defining the elements")

order=0
nDoFs=0
if itype in ["P1","PGL1"]:
    order=1
elif itype in ["B2","P2","PGL2"]:
    order=2
elif itype in ["P3","B3","PGL3"]:
    order=3
elif itype in ["B4","P4","PGL4"]:
    order=4
else:
    print("Element not defined")
    sys.exit()

nDoFs=order+1
print("Element: ",itype)
print("Order of the polynomials: ",order)
print("DoFs per cell: ",nDoFs)

class element:
    left=0.
    right=0.
    coord=[]

print("Defining the mesh")
dx = L/n_el #length of the element
#print(dx)

mesh=[]
#generation of the extrema
for indi in range(n_el):
    el=element()
    el.left=dx*indi
    el.right=dx*(indi+1)
    mesh.append(el)

#for indi in range(n_el):
#    print(indi, mesh[indi].left,mesh[indi].right)    
#    print()

#for indi in range(n_el):
#    print(indi)
#    print(format(mesh[indi].left, '.15f'),format(mesh[indi].right, '.15f')) #me
#    print("{:.15f}".format(mesh[indi].left),"{:.15f}".format(mesh[indi].right)) #francesco
#    print()

#generation of the internal DoFs

#########################################3
#experiment with linspace
#z=np.linspace(0, 1, 10)
#for indi in range(len(z)):
#    print(z[indi])
#    print("{:.15f}".format(z[indi]))
#    print()
#z=np.linspace(0, 1, 3)
#for indi in range(len(z)):
#    print(z[indi])
#    print("{:.15f}".format(z[indi]))
#    print()
#z=np.linspace(0, 1, 1) #only 0.
#for indi in range(len(z)):
#    print(z[indi])
#    print("{:.15f}".format(z[indi]))
#    print()
#z=np.linspace(0, 1, 0) #empty
#for indi in range(len(z)):
#    print(z[indi])
#    print("{:.15f}".format(z[indi]))
#    print()
#########################################3

for indi in range(n_el):
    xmid=(mesh[indi].left+mesh[indi].right)/2.
    if itype in ["P1","PGL1"]: #2 DoFs
        mesh[indi].coor=[mesh[indi].left, mesh[indi].right] 
    elif itype in ["B2","P2","PGL2"]: #3 DoFs
        mesh[indi].coor=np.linspace(mesh[indi].left, mesh[indi].right, 3)
    elif itype in ["P3","B3"]: #4 DoFs
        mesh[indi].coor=np.linspace(mesh[indi].left, mesh[indi].right, 4)
    elif itype in ["P4","B4"]: #5 DoFs
        mesh[indi].coor=np.linspace(mesh[indi].left, mesh[indi].right, 5)        
    elif itype in ["PGL3"]: #4 DoFs but not equispaced
        #e%x(1,3)=0.5_dp-SQRT(5._dp)/10._dp
        #e%x(1,4)=0.5_dp+SQRT(5._dp)/10._dp
        a=np.sqrt(5.)/10.
        mesh[indi].coor=[mesh[indi].left, xmid-dx*a, xmid+dx*a, mesh[indi].right]
        #no precision issues
        #print(np.sqrt(5.)/10.)
        #print(np.sqrt(5)/10)
        #print(np.sqrt(5)/10-0.2236067977499789696409)
        #print(np.sqrt(5)/10-np.sqrt(5.)/10.)
    elif itype in ["PGL4"]: #5 DoFs
        #e%x(1,3)=0.5_dp-SQRT(21._dp)/14._dp
        #e%x(1,4)=0.5_DP
        #e%x(1,5)=0.5_dp+SQRT(21._dp)/14._dp
        b=np.sqrt(21.)/14.
        mesh[indi].coor=[mesh[indi].left, xmid-dx*b, xmid, xmid+dx*b, mesh[indi].right]
    else:
        print("Element not defined")
        sys.exit()
    
#for indi in range(n_el):
#    print(indi, mesh[indi].coor)    
#    print()
#print()
#for indi in range(n_el):
#    print(indi, mesh[indi].coor[0],mesh[indi].coor[1],mesh[indi].coor[2],mesh[indi].coor[3])    
#    print()

#-----------------------------------------------------------------------------
#Discrete solution that later will be interpolated in the coarse mesh that we have built
print("Acquisition of the discrete solution of the refined mesh")


x_ref=[]
h_ref=[]
hu_ref=[]

eta_ref=[] #<--!NEW I need eta to interpolate it

data_refined=[]
with open(file_refined+".dat", 'r') as a:
    for idx, line in enumerate(a.readlines()):
        data_refined.append(line)

#loop over the lines that we have copied
for idx, line in enumerate(data_refined): #rmk: 0-based numeration
        if idx==0: #first line we skip, it is the header
            pass
        else:
            data_line=line.split()
            data_indi=int(data_line[0])
            data_x=float(data_line[1])
            data_h=float(data_line[2])
            data_hu=float(data_line[3])
            #print("line", line)
            #print("Just read: ",data_indi,data_x,data_h,data_hu)
            #print()
            x_ref.append(data_x)
            h_ref.append(data_h)
            hu_ref.append(data_hu)
            

#print(x_ref)
#print(h_ref)
#print(hu_ref)
#print(len(x_ref))

#Filling eta to interpolate it
for indi in range(len(x_ref)):
    etasupp=h_ref[indi]+bathymetry(x_ref[indi])
    eta_ref.append(etasupp)
    #print(indi, x_ref[indi], eta_ref[indi])
    
#-----------------------------------------------------------------------------
#(Re)construction of the refined mesh
print("Generation of the refined mesh")

print("Defining the elements of the mesh refined")

order_ref=0
nDoFs_ref=0
if itype_ref in ["P1","PGL1"]:
    order_ref=1
elif itype_ref in ["B2","P2","PGL2"]:
    order_ref=2
elif itype_ref in ["P3","B3","PGL3"]:
    order_ref=3
elif itype_ref in ["B4","P4","PGL4"]:
    order_ref=4
else:
    print("Element of the refined mesh not defined")
    sys.exit()

nDoFs_ref=order_ref+1
print("Element_ref: ",itype_ref)
print("Order_ref of the polynomials: ",order_ref)
print("DoFs_ref per cell: ",nDoFs_ref)


print("Defining the mesh")
n_el_ref=int((len(x_ref)-1)/order_ref)
dx_ref = L/n_el_ref #length of the element
#print(dx_ref, n_el_ref)


#defining the element of the mesh refined
class element_ref:
    left=0.
    right=0.
    coord=[]
    h=[]
    hu=[]
    eta=[] #<--!NEW, I need eta to interpolate it



mesh_ref=[]
#generation of the extrema
for indi in range(n_el_ref):
    el_ref=element_ref()
    el_ref.left=dx_ref*indi
    el_ref.right=dx_ref*(indi+1)
    mesh_ref.append(el_ref)

#for indi in range(n_el_ref):
#    print(indi, mesh_ref[indi].left,mesh_ref[indi].right)    
#    print()


#generation of the internal DoFs


for indi in range(n_el_ref):
    xmid=(mesh_ref[indi].left+mesh_ref[indi].right)/2.
    if itype_ref in ["P1","PGL1"]: #2 DoFs
        mesh_ref[indi].coor=[mesh_ref[indi].left, mesh_ref[indi].right] 
    elif itype_ref in ["B2","P2","PGL2"]: #3 DoFs
        mesh_ref[indi].coor=np.linspace(mesh_ref[indi].left, mesh_ref[indi].right, 3)
    elif itype_ref in ["P3","B3"]: #4 DoFs
        mesh_ref[indi].coor=np.linspace(mesh_ref[indi].left, mesh_ref[indi].right, 4)
    elif itype_ref in ["P4","B4"]: #5 DoFs
        mesh_ref[indi].coor=np.linspace(mesh_ref[indi].left, mesh_ref[indi].right, 5)        
    elif itype_ref in ["PGL3"]: #4 DoFs but not equispaced
        #e%x(1,3)=0.5_dp-SQRT(5._dp)/10._dp
        #e%x(1,4)=0.5_dp+SQRT(5._dp)/10._dp
        a=np.sqrt(5.)/10.
        mesh_ref[indi].coor=[mesh_ref[indi].left, xmid-dx_ref*a, xmid+dx_ref*a, mesh_ref[indi].right]
        #no precision issues
        #print(np.sqrt(5.)/10.)
        #print(np.sqrt(5)/10)
        #print(np.sqrt(5)/10-0.2236067977499789696409)
        #print(np.sqrt(5)/10-np.sqrt(5.)/10.)
    elif itype_ref in ["PGL4"]: #5 DoFs
        #e%x(1,3)=0.5_dp-SQRT(21._dp)/14._dp
        #e%x(1,4)=0.5_DP
        #e%x(1,5)=0.5_dp+SQRT(21._dp)/14._dp
        b=np.sqrt(21.)/14.
        mesh_ref[indi].coor=[mesh_ref[indi].left, xmid-dx_ref*b, xmid, xmid+dx_ref*b, mesh_ref[indi].right]
    else:
        print("Element not defined")
        sys.exit()
        
#For interpolation purposes
reference_element_coordinates=[]
if itype_ref in ["P1","PGL1"]: #2 DoFs
    reference_element_coordinates=[0., 1.] 
elif itype_ref in ["B2","P2","PGL2"]: #3 DoFs
    reference_element_coordinates=np.linspace(0., 1., 3)
elif itype_ref in ["P3","B3"]: #4 DoFs
    reference_element_coordinates=np.linspace(0., 1., 4)
elif itype_ref in ["P4","B4"]: #5 DoFs
    reference_element_coordinates=np.linspace(0., 1., 5)        
elif itype_ref in ["PGL3"]: #4 DoFs but not equispaced
     #e%x(1,3)=0.5_dp-SQRT(5._dp)/10._dp
     #e%x(1,4)=0.5_dp+SQRT(5._dp)/10._dp
     a=np.sqrt(5.)/10.
     reference_element_coordinates=[0., 0.5-a, 0.5+a, 1.]
elif itype_ref in ["PGL4"]: #5 DoFs
     #e%x(1,3)=0.5_dp-SQRT(21._dp)/14._dp
     #e%x(1,4)=0.5_DP
     #e%x(1,5)=0.5_dp+SQRT(21._dp)/14._dp
     b=np.sqrt(21.)/14.
     reference_element_coordinates=[0., 0.5-b, 0.5, 0.5+b, 1.]
else:
     print("Element not defined")
     sys.exit()
        



        
#Acquisition of the solution in the mesh for interpolation purposes
for indi in range(n_el_ref):
    mesh_ref[indi].h=h_ref[order_ref*indi:order_ref*indi+order_ref+1] #RMK: +1 to include the last one
    mesh_ref[indi].hu=hu_ref[order_ref*indi:order_ref*indi+order_ref+1] #RMK: +1 to include the last one
    mesh_ref[indi].eta=eta_ref[order_ref*indi:order_ref*indi+order_ref+1] #<--!NEW
    
    #for indj in range(nDoFs_ref):
    #    print(mesh_ref[indi].eta[indj]-mesh_ref[indi].h[indj]-bathymetry(mesh_ref[indi].coor[indj]))
    
#for indi in range(n_el_ref):
#    #print(indi) 
#    #print(mesh_ref[indi].coor)    
#    #print(x_ref[order_ref*indi:order_ref*indi+order_ref+1])
#    #print(mesh_ref[indi].h)
#    #print(h_ref[order_ref*indi:order_ref*indi+order_ref+1])
#    #print(mesh_ref[indi].hu)
#    #print(hu_ref[order_ref*indi:order_ref*indi+order_ref+1])
    
#print()
#for indi in range(n_el_ref):
#    print(indi, mesh_ref[indi].coor[0],mesh_ref[indi].coor[1],mesh_ref[indi].coor[2],mesh_ref[indi].coor[3])    
#    print()
#
#print()
#print("x",x_ref)
#print("h",h_ref)
#print("hu",hu_ref)

#print(n_el_ref)
#print(n_el)
#sys.exit()
#-----------------------------------------------------------------------------
#Now let's process the abscissae of the coarse mesh with respect to the refined mesh
print("Interpolation")
xfinal=np.zeros(order*n_el+1)
hfinal=np.zeros(len(xfinal))
hufinal=np.zeros(len(xfinal))
#print(len(xfinal),len(hfinal),len(hufinal))



#We will deal with the right extremum of the domain singularly. We will loop over all the cells of the coarse mesh to store the interpolated solution in all the coor apart from the last one

indi=0

for elem in mesh: #loop over the cells of the coarse mesh
    for indc in range(nDoFs-1): #loop over the coordinates of the coarse mesh apart from the last #rmk:zero based numeration
        x=elem.coor[indc]
        #let's try to understand where the point of the coarse mesh is in the refined mesh
        index_el_ref=int(x//dx_ref)
        e_ref=mesh_ref[index_el_ref]
        #print("Coordinate", x)
        #print("In element", index_el_ref)
        #print("With coordinates", e_ref.coor)
        #print()
        
        #Let's try to understand the normalized coordinate in that element
        norm_coor=(x-e_ref.left)/dx_ref
        #print(norm_coor)
        
        x_lagrange=reference_element_coordinates
        eta_lagrange=e_ref.eta #<--!NEW I interpolate eta rather than h
        hu_lagrange=e_ref.hu
        
        #print(x_lagrange)
        #print(eta_lagrange)
        #print(x)
        #print(index_el_ref)
        
        etafunction=lagrange(x_lagrange,eta_lagrange) #<--!NEW I interpolate eta rather than h
        hufunction=lagrange(x_lagrange,hu_lagrange)
        
        eta_interpolated=etafunction(norm_coor)
        hu_interpolated=hufunction(norm_coor)

        xfinal[indi]=x
        #print(indc, elem.coor[indc], bathymetry(elem.coor[indc]))
        hfinal[indi]=eta_interpolated-bathymetry(elem.coor[indc]) #<--!NEW I interpolate eta rather than h, but then I have to subtract the bathymetry
        hufinal[indi]=hu_interpolated

        
        indi=indi+1

xfinal[len(xfinal)-1]=x_ref[len(x_ref)-1]
hfinal[len(xfinal)-1]=h_ref[len(x_ref)-1] #For the final point no inteprolation problem
hufinal[len(xfinal)-1]=hu_ref[len(x_ref)-1]




def perturbationfunction(x,indp):
    x0=6; r=0.5
    if indp==0:
        A=0
    elif indp==1:
        A=5*10**(-2)
    elif indp==2:
        A=5*10**(-4)
    elif indp==3:
        A=5*10**(-5)
    else:
        print('perturbation not recognised')
        quit()
    p=0.
    if(x>x0-r and x<x0+r):
        p=A*np.exp(1-1/(1-((x-x0)/r)**2))

    return p


bfinal=np.zeros(len(xfinal))
for indi in range(len(bfinal)):
    bfinal[indi]=bathymetry(xfinal[indi])



etafinal=np.zeros(len(xfinal))
for indi in range(len(bfinal)):
    etafinal[indi]=hfinal[indi]+bfinal[indi]+100*perturbationfunction(xfinal[indi],3)




fig=pl.figure()
params = {'mathtext.default': 'regular' }   
pl.plot(xfinal,etafinal,"-",linewidth=2,label='$\eta$')
#pl.plot(xfinal,bfinal,"-",linewidth=2,label='B')
pl.grid()
pl.xlabel("x")
#pl.ylabel("y")
pl.legend(loc='upper right')
#pl.ylim([0, 2.5])
pl.savefig("friction_sub.pdf", format="pdf", bbox_inches="tight")
pl.show()


quit()




#print(xfinal)
#print(hfinal)
#print(hufinal)

#figh = pl.figure()
#pl.plot(x_ref,h_ref)
#pl.plot(xfinal,hfinal)
#pl.show()

#fighu = pl.figure()
#pl.plot(x_ref,hu_ref)
#pl.plot(xfinal,hufinal)
#pl.show()



quit()

#sys.exit()
#-----------------------------------------------------------------------------
#SOLUTION THROUGH THE IMPLICIT SOLVER
print("Saving")
f=open(file_refined+"Interpolated.dat","w+")
f.write("  indi    x    h   hu\n")
for indi in range(len(xfinal)):
    towrite=" "+str(indi+1)+"     "+format(xfinal[indi], '.15f')+"    "+format(hfinal[indi], '.15f')+"    "+format(hufinal[indi], '.15f')+"\n"
    f.write(towrite)
f.close()



sys.exit()



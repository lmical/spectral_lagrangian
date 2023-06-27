itype="PGL4" #P1, B2, P2, P3, B3, B4, P4, PGL2, PGL3, PGL4
n_el=512
#-----------------------------------------------------------------------------
#Precision
from decimal import Decimal
from mpmath import *
import gmpy2
gmpy2.get_context().precision=100
mp.dps = 100
from sympy import *
#getcontext().prec = 50 # 34 or 128 ?
#print("Parameters of the test chosen")
#-----------------------------------------------------------------------------
#Length domain
L=Decimal('25.')
#print("Length of the domain", L)
#-----------------------------------------------------------------------------
#Physical parameters
g=Decimal('9.81')
q0=Decimal('24.')
hL=Decimal('2.')
#print("Physical parameters", g, q0, hL)
#-----------------------------------------------------------------------------
import numpy as np 
from numpy import linalg as LA
from scipy import optimize as opt
import matplotlib.pyplot as pl
import sys 
import itertools
#-----------------------------------------------------------------------------
#Bathymetry
def bathymetry(x):
    x0=Decimal('10.')
    r=Decimal('5.')
    b=Decimal('0.')
    if (x > x0-r) and (x < x0+r):
        b=Decimal('0.2')*np.exp(Decimal('1.') - Decimal('1.')/(Decimal('1.')-((x-x0)/r)**Decimal('2.')))    
    else:
        b=Decimal('0.')
    return b
#xplot = np.arange(0.,L,0.1) #plot abscissa
#bplot = np.zeros(np.size(xplot)) #plot bathymetry
#for indi in range(len(xplot)):
#    bplot[indi]=bathymetry(xplot[indi])
#fig1 = pl.figure()
#pl.plot(xplot,bplot)
#pl.show()
#sys.exit()
#-----------------------------------------------------------------------------
#Construction of the mesh
print("Defining the elements")
order=0
nDoFs=0
if itype in ["P1"]:
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
    left=Decimal('0.')
    right=Decimal('0.')
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
    xmid=(mesh[indi].left+mesh[indi].right)/Decimal('2.')
    if itype in ["P1"]: #2 DoFs
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
        a=np.sqrt(Decimal('5.'))/Decimal('10.')
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
        b=np.sqrt(Decimal('21.'))/Decimal('14.')
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
#SOLUTION THROUGH THE IMPLICIT SOLVER
print("Solution of the implicit problem")
xfinal=np.zeros(order*n_el+1)
hfinal=np.zeros(len(xfinal))
#print(len(xfinal),len(hfinal))

#We will deal with the left extremum of the domain singularly and then we will loop over all the cells to store the solution to the problem in all the coor apart from the first one
xfinal[0]=Decimal('0.')
hfinal[0]=hL

indi=1

for elem in mesh: #loop over the cells
    for indc in range(1,nDoFs): #loop over the coordinates apart from the first
        x=elem.coor[indc]
        xfinal[indi]=x
        b=bathymetry(x)#; print(b)
        #Newton
        #implicitequation = lambda h: h**3 + (b-q0**2/(2.*g*hL**2) - hL)*h**2 + q0**2/(2.*g)
        #h = opt.newton(implicitequation,hL)
        #Exact
        #p=[Decimal('1.'), (b-q0**Decimal('2')/(Decimal('2.')*g*hL**Decimal('2.')) - hL), Decimal('0.'), q0**Decimal('2.')/(Decimal('2.')*g)]
        from sympy.abc import z
        hvec = nroots( z**3+(b-q0**2/(2*g*hL**2) - hL)*z**2+q0**2/(2*g), n=30)
        h=hvec[1] #; print(h)
        #Check the difference    
        #print(indi,hvec[:])   
        hfinal[indi]=h
        indi=indi+1

#print(xfinal)
#print(hfinal)
#figh = pl.figure()
#pl.plot(xfinal,hfinal)
#pl.show()
#sys.exit()
#-----------------------------------------------------------------------------
#SOLUTION THROUGH THE IMPLICIT SOLVER
print("Saving")
f=open("supercriticalsmooth.dat","w+")
f.write("  indi    x    h   hu\n")
for indi in range(len(xfinal)):
    towrite=" "+str(indi+1)+"     "+format(xfinal[indi], '.34f')+"    "+format(hfinal[indi], '.34f')+"    "+format(q0, '.34f')+"\n"
    f.write(towrite)
f.close()


sys.exit()



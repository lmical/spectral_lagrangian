import numpy as np 
from numpy import linalg as LA
from scipy import optimize as opt
import matplotlib.pyplot as pl
import sys 
import itertools


xfinal=np.linspace(0,25,1000)

#Bathymetry
def bathymetry(x):
    b=0.
    if (x > 8) and (x < 12):
        b=0.2 - 0.05*(x - 10)**2    
    else:
        b=0.
    return b

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
    etafinal[indi]=0.5+1000*perturbationfunction(xfinal[indi],3)

fig=pl.figure()
params = {'mathtext.default': 'regular' }   
pl.plot(xfinal,etafinal,"-",linewidth=2,label='$\eta$')
pl.plot(xfinal,bfinal,"-",linewidth=2,label='B')
pl.grid()
pl.xlabel("x")
#pl.ylabel("y")
pl.legend(loc='upper right')
pl.ylim([0, 1.])
pl.savefig("alb_latr_pert_IC.pdf", format="pdf", bbox_inches="tight")
pl.show()



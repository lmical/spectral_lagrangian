import numpy as np 
from numpy import linalg as LA
from scipy import optimize as opt
import matplotlib.pyplot as pl
import sys 
import itertools


xfinal=np.linspace(0,25,1000)

#Bathymetry
def bathymetry(x):
    x0=10.
    r=2.
    b=0.
    if (x > x0-r) and (x < x0+r):
        b=0.2-0.05*(x-10)**2    
    else:
        b=0.
    return b



bfinal=np.zeros(len(xfinal))
for indi in range(len(bfinal)):
    bfinal[indi]=bathymetry(xfinal[indi])

fig=pl.figure()
params = {'mathtext.default': 'regular' }   
pl.plot(xfinal,0.5*np.ones(len(xfinal)),"-",linewidth=2,label='$\eta$')
pl.plot(xfinal,bfinal,"-",linewidth=2,label='B')
pl.grid()
pl.xlabel("x")
#pl.ylabel("y")
pl.legend(loc='upper right')
pl.ylim([0, 1])
pl.savefig("alb_LAKEATREST.pdf", format="pdf", bbox_inches="tight")
pl.show()

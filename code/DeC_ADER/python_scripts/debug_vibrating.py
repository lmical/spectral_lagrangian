from src import DeC
from src import RungeKutta
from src import ODEproblems
from src.DeC import *
from src.ODEproblems import ODEproblem

dec_method=DeC(3,4,"equispaced")
#Test convergence of DeC for several orders
pr=ODEproblem("vibratingDamped")

tt=np.linspace(0,pr.T_fin,1000)   #Plot the evolution for order 8
tt,uu=dec_method.dec(pr.flux, tt, pr.u0)

u_ex = pr.exact_solution_times(pr.u0,tt)
    

plt.plot(tt,uu[0,:],label="approx")
#plt.plot(tt,uu[1,:],label="speed")
plt.plot(tt,u_ex[0,:],label="exact")
plt.legend()
plt.show()

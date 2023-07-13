import numpy as np
import matplotlib
matplotlib.use('tkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
import sys
import os
import re
from optparse import OptionParser

#Command line parser/////////////////////////////////////////
parser = OptionParser()
parser.add_option("-o", "--output", action="store", type="int", dest="o", default ="0")
parser.add_option("-d", "--dir",    action="store", type="string", dest="d", default="outputs")
options, args = parser.parse_args()
output = options.o
folder = options.d

print(folder)

try:
    nt,itype,ischema,test = np.loadtxt(folder+"/Descriptor.dat",usecols=(0,1,2,3),dtype=int,unpack=True)
    print(nt)
    print('nx = %d \n' %nt)
    print('itype = %d \n' %itype)
    print('schema = %d \n' %ischema)
    print('test case = %d \n' %test)
except:
    print("Couldn't open descriptor. To start, type: python viewer1d.py --dir [dir]")
    exit()

fig = plt.figure()
ax1 = plt.subplot2grid((3,1), (0,0))
ax2 = plt.subplot2grid((3,1), (1,0))
ax3 = plt.subplot2grid((3,1), (2,0))

n=0

def visualize(i):
    global n, ordered_list
    n=i
    nr = ordered_list[n]
    name=folder+"/sol_prim_"+str(nr)+".dat"
    try:
        Ma=np.loadtxt(name)
    except:
        print("Couldn't open file: "+name)
        return
    print(Ma)
    rho = Ma[:,1]
    u = Ma[:,2]
    p = Ma[:,3]
    x = Ma[:,0]
    ax1.cla()
    ax1.plot(x,rho,"b")
    ax1.set_title('rho')
    ax2.cla()
    ax2.plot(x,u,"b")
    ax2.set_title('u')
    ax3.cla()
#    ax3.plot(x,p-np.exp(-x),"b")
    ax3.plot(x,p,"b")
    ax3.set_title('p')
    fig.suptitle('Right +1, left -1, down +10, up -10, end NT, home 0, timestep '+str(nr))

def change_file(event):
    sys.stdout.flush()
    global n,cbar,field
    global ordered_list
    if event.key=='right':
        n+=1
    elif event.key=='left':
        if n>0:
            n-=1
        else:
            return
    elif event.key == 'down':
        n=min(n+10, len(list_of_files)-1)
    elif event.key == 'up':
        n=max(n-10,0)
    elif event.key == 'end':
        n = len(list_of_files)-1
    elif event.key == 'home':
        n = 0

    visualize(n)
    fig.canvas.draw_idle()

def animate(i):
    return visualize(i)

if __name__=='__main__':
    list_of_files = [f for f in os.listdir(folder) if 'sol_prim_' in f and 'sol_prim_last' not in f and 'sol_prim_nodal' not in f]
    ordered_list = sorted([int(re.split('_|\.',f)[2]) for f in list_of_files])
    visualize(n)

    fig.canvas.mpl_connect('key_press_event', change_file)
    plt.show()

#    if video_flag == True:
#        range_of_files = np.arange(0,len([f for f in os.listdir(folder) if 'modes' in f])-1)
#        for var in ['rho','vx','vy','p']:
#            field = var
#            min_y, max_y = bounds[field]
#            anim = animation.FuncAnimation(fig,animate,range_of_files,blit=True)
#            anim.save(folder + '_' + field+'.mp4',writer = 'ffmpeg',fps=15)

import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 200                           # Number of points in each direction
xStart,xEnd = -4.0,4.0            # x-direction boundaries
yStart,yEnd = -2.0,2.0            # y-direction boundaries
x = np.linspace(xStart,xEnd,N)    # x 1D-array
y = np.linspace(yStart,yEnd,N)    # y 1D-array
X,Y = np.meshgrid(x,y)            # generation of the mesh grid

#np.shape(X)

uinf=1.0
alphad=0.0
alpha=alphad*pi/180

ufreestream=uinf*cos(alpha)*np.ones((N,N),dtype=float)
vfreestream=uinf*sin(alpha)*np.ones((N,N),dtype=float)

psifreestream= uinf*cos(alpha)*Y-uinf*sin(alpha)*X

def getVelocity(sources,xs,ys,X,Y):
    u = sources/(2*pi)*(X-xs)/((X-xs)**2+(Y-ys)**2)
    v = sources/(2*pi)*(Y-ys)/((X-xs)**2+(Y-ys)**2)
    return u,v

def getStreamFunction(sources,xs,ys,X,Y):
    psi = sources /(2*pi)*np.arctan2((Y-ys),(X-xs))
    return psi
    
sources=5.0
xsource,ysource=-1.0,0.0

usource,vsource=getVelocity(sources,xsource,ysource,X,Y)

psisource= getStreamFunction(sources,xsource,ysource,X,Y)

u=ufreestream+usource
v=vfreestream+vsource
psi=psifreestream+psisource

##stagnation
xstag=xsource-((sources/(2*pi*uinf))*cos(alpha))
ystag=ysource-((sources/(2*pi*uinf))*sin(alpha))
psistag= sources /(2*pi)*np.arctan2((ystag-ysource),(xstag-xsource))
###PLOT
size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)

plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xsource,ysource,c='#CD2305',s=80,marker='o')
plt.scatter(xstag,ystag,c='#006633',s=80,marker='o')
if (alpha==0.0):
    plt.contour(X,Y,psi,\
            levels=[-psistag,+psistag],\
            colors='#CD2305',linewidths=2,linestyles='solid')

sinks = -5.0     # strength of the sink
xsink,ysink = 1.0,0.0   # location of the sink

usink,vsink = getVelocity(sinks,xsink,ysink,X,Y)

psisink = getStreamFunction(sinks,xsink,ysink,X,Y)

u=ufreestream+usource+usink
v=vfreestream+vsource+vsink
psi=psifreestream+psisource+psisink

size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,u,v,density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter([xsource,xsink],[ysource,ysink],c='#CD2305',s=80,marker='o')
if (alpha==0.0):
    plt.contour(X,Y,psi,levels=[0.0],colors='#CD2305',linewidths=2,linestyles='solid')
    
    
Cp=1.0-(u**2+v**2)/uinf**2

size = 10
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label(r'$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.scatter([xsource,xsink],[ysource,ysink],c='#CD2305',s=80,marker='o')
plt.contour(X,Y,psi,\
            levels=[0.0],\
            colors='#CD2305',linewidths=2,linestyles='solid')
        
plt.show()
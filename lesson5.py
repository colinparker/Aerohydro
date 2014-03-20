import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50
xStart, xEnd = -2.0,2.0
yStart, yEnd = -0.5,0.5
x=np.linspace(xStart, xEnd, N)
y=np.linspace(yStart, yEnd, N)
X,Y=np.meshgrid(x,y)

gamma = 5.0

def getVelocityVortex(strength,xv,yv,X,Y):
    u = + strength/(2*pi)*(Y-yv)/((X-xv)**2+(Y-yv)**2)
    v = - strength/(2*pi)*(X-xv)/((X-xv)**2+(Y-yv)**2)
    return u,v
    
def getStreamFunctionVortex(strength,xv,yv,X,Y):
    psi = strength/(4*pi)*np.log((X-xv)**2+(Y-yv)**2)
    return psi
    
##for line at y=0  
nV=10  ##number of vorticies
dxv=0.1 ##spacing
yvortex = 0
ufield=0
vfield=0
psifield=0
for i in range(0,nV):
    xvortex=-(dxv*nV)*0.5+i*dxv
    uvortex,vvortex = getVelocityVortex(gamma,xvortex,yvortex,X,Y)
    psivortex=getStreamFunctionVortex(gamma,xvortex,yvortex,X,Y)
    ufield=ufield+uvortex
    vfield=vfield+vvortex
    psifield=psifield+psivortex
XV=np.linspace(xStart, xEnd, nV)
YV=np.linspace(0,0,nV)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,ufield,vfield,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')

plt.scatter(XV,YV,c='#663399',s=10,marker='o');
plt.show()
print YV
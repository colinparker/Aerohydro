import numpy as np
import matplotlib.pyplot as plt
from math import *

N=50
xStart,xEnd = -2.0,2.0            
yStart,yEnd = -1.0,1.0            
x = np.linspace(xStart,xEnd,N)    
y = np.linspace(yStart,yEnd,N)    
X,Y = np.meshgrid(x,y)  

kappa=1.0
xDoublet,yDoublet=0.0,0.0

def getVelocityDoublet(strength,xd,yd,X,Y):
    u= -strength/(2*pi)*((X-xd)**2-(Y-yd)**2)/((X-xd)**2+(y-yd)**2)**2
    v= -strength/(2*pi)*2*(X-xd)*(Y-yd)/((X-xd)**2+(Y-yd)**2)**2
    return u,v

def getStreamFunctionDoublet(strength,xd,yd,X,Y):
    psi= -strength/(2*pi)*(Y-yd)/((X-xd)**2+(Y-yd)**2)
    return psi

uDoublet,vDoublet=getVelocityDoublet(kappa,xDoublet,yDoublet,X,Y)

psiDoublet=getStreamFunctionDoublet(kappa,xDoublet,yDoublet,X,Y)

size = 10

plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.streamplot(X,Y,uDoublet,vDoublet,\
               density=2.0,linewidth=1,arrowsize=1,arrowstyle='->')
plt.scatter(xDoublet,yDoublet,c='#CD2305',s=80,marker='o')
plt.show()
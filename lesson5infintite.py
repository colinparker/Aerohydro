import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50
xstart, xend = -2.0,2.0
ystart, yend = -1.0,1.0
x=np.linspace(xstart, xend, N)
y=np.linspace(ystart, yend, N)
X,Y=np.meshgrid(x,y)

a=.5
gamma = 5.0

u=gamma/(2*a)*np.sinh(2*pi*Y/a)/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a))
v=gamma/(2*a)*np.sin(2*pi*X/a)/(np.cosh(2*pi*Y/a)-np.cos(2*pi*X/a))

size = 10
plt.figure(figsize=(size,(yend-ystart)/(xend-xstart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart,xend)
plt.ylim(ystart,yend)
plt.streamplot(X,Y,u,v,\
         density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.title('Infinite Row of Vortices', fontsize=16)
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from math import *

N = 50
xstart, xend = -2.0,2.0
ystart, yend = -1.0,1.0
x=np.linspace(xstart, xend, N)
y=np.linspace(ystart, yend, N)
X,Y=np.meshgrid(x,y)

#size = 2
#plt.figure(figsize=((size*(xend-xstart)),(size*(yend-ystart))))
#plt.xlabel('x',fontsize=16)
#plt.ylabel('y',fontsize=16)
#plt.xlim(xstart*1.1,xend*1.1)
#plt.ylim(ystart*1.1,yend*1.1)
#plt.scatter(X,Y,s=8,c='#339900',marker='o',linewidth=0.1)


ss = 5.0
xsource, ysource= -1.0,0.0
usource = np.empty((N,N),dtype=float)
vsource = np.empty((N,N),dtype=float)

for i in range(N):
    for j in range(N):
        usource[i,j] = ss/(2*pi)*(X[i,j]-xsource)/((X[i,j]-xsource)**2+(Y[i,j]-ysource)**2)
        vsource[i,j]=ss/(2*pi)*(Y[i,j]-ysource)/((X[i,j]-xsource)**2+(Y[i,j]-ysource)**2)

#size=3
#plt.figure(figsize=((size*(xend-xstart)),(size*(yend-ystart))))
#plt.xlabel('x',fontsize=16)
#plt.ylabel('y',fontsize=16)
#plt.xlim(xstart*1.1,xend*1.1)
#plt.ylim(ystart*1.1,yend*1.1)
#plt.streamplot(X,Y,usource,vsource,density=2.0,linewidth=1,arrowsize=2,arrowstyle='->')
#plt.scatter(xsource,ysource,c='#111111',s=80,marker='o')

uinf=2.0
sinks=-5.0
xsink, ysink= 1.0,0.0
usink = np.empty((N,N),dtype=float)
vsink = np.empty((N,N),dtype=float)
usink=sinks/(2*pi)*(X-xsink)/((X-xsink)**2+(Y-ysink)**2)
vsink=sinks/(2*pi)*(Y-ysink)/((X-xsink)**2+(Y-ysink)**2)
size=3

plt.figure(figsize=((size*(xend-xstart)),(size*(yend-ystart))))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart*1.1,xend*1.1)
plt.ylim(ystart*1.1,yend*1.1)
plt.streamplot(X,Y,usink,vsink,density=2.0,linewidth=1,arrowsize=2,arrowstyle='->')
plt.scatter(xsink,ysink,c='#111111',s=80,marker='o')
plt.show()


upair=np.empty_like(usource)
vpair=np.empty_like(vsource)
Upair=np.empty_like(vsource)
Z=np.empty_like(vsource)
Z=ss/(2*pi)*np.log(np.sqrt((X-xsource)**2+(Y-ysource)**2))+\
            sinks/(2*pi)*np.log(np.sqrt((X-xsink)**2+(Y-ysink)**2))
upair=usource+usink+uinf
vpair=vsource+vsink
for i in range(N):
    for j in range(N):
        Upair[i][j]=sqrt(upair[i][j]**2+vpair[i][j]**2)
size=3

plt.figure(figsize=((size*(xend-xstart)),(size*(yend-ystart))))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xstart*1.1,xend*1.1)
plt.ylim(ystart*1.1,yend*1.1)
#plt.streamplot(X,Y,upair,vpair)
#plt.scatter(xsink,ysink,c='#111111',s=80,marker='o')
#plt.scatter(xsource,ysource,c='#336699',s=80,marker='o')
plt.contourf(X,Y,Z,1000)
plt.show()

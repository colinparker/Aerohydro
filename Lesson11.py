import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

coords = np.loadtxt(fname='../Aerohydro/resources/naca0012.dat')
xp,yp=coords[:,0],coords[:,1]

valX,valY=.1,.2
xmin,xmax=min(xp),max(xp)
ymin,ymax=min(yp),max(yp)
xStart,xEnd=xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd=ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid=(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k',linewidth=2);

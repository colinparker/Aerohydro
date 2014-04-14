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

class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya=xa,ya
        self.xb,self.yb=xb,yb
        self.xc,self.yc=(xa+xb)/2,(ya+yb)/2
        self.length=sqrt((xb-xa)**2+(yb-ya)**2)
        
        if(xb-xa<=0):self.beta=acos((yb-ya)/self.length)
        elif (xb-xa>0):self.beta=pi+acos(-(yb-ya)/self.length)
        
        if(self.beta<=pi):self.loc='extrados'
        else:self.loc='intrados'
        
        self.sigma=0.
        self.vt=0.
        self.Cp=0.

def definePanels(N,xp,yp):
    
    R=(max(xp)-min(xp))/2
    xCenter =(max(xp)+min(xp))/2
    xCircle=xCenter + R*np.cos(np.linspace(0,2*pi,N+1))
    
    x=np.copy(xCircle)
    y=np.empty_like(x)
    
    xp,yp=np.append(xp,xp[0]),np.append(yp,yp[0])
    
    I=0
    for i in range(N):
        while(I<len(xp)-1):
            if(xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]):break
            else:I+=1
        a=(yp[I+1]-yp[I])/(xp[I+1]-xp[I])
        b=yp[I+1]-a*xp[I+1]
        y[i]=a*x[i]+b
    y[N]=y[0]
    
    panel=np.empty(N,dtype=object)
    for i in range(N):
        panel[i]=Panel(x[i],y[i],x[i+1],y[i+1])
    return panel

N=20
panel=definePanels(N,xp,yp)

valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
#plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel],panel[0].ya),\
        linestyle='-',linewidth=1,\
        marker='o',markersize=6,color='#CD2305');

class Freestream:
    def __init__(self,Uinf,alpha):
        self.Uinf=Uinf
        self.alpha=alpha*pi/180
Uinf=1.0
alpha=1.0
freestream=Freestream(Uinf,alpha)

def I(xci,yci,pj,dxdz,dydz):
    def func(s):
        return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
            +(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
            /((xci-(pj.xa-sin(pj.beta)*s))**2\
            +(yci-(pj.ya+cos(pj.beta)*s))**2)
    return integrate.quad(lambda s:func(s),0.,pj.length)[0]
    
def sourceMatrix(p):
    N=length(p)
    A=np.empty((N,N),dtype=float)
    np.fill_diagonal(A,0.5)
    for i in range(N):
        for j in range(N):
            if(j!=i):
                B[i]-=0.5/pi*I(p[i].xc,p[i].yc,p[j],+sin(p[i].beta),-cos(p[i].beta))
    return B
    

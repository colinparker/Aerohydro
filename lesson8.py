import numpy as np
from scipy import integrate
from math import *
import matplotlib.pyplot as plt

# function to read the coordinates to store into two 1D arrays
def readGeometry():
    inFile = open('C:/Users/Colin/Documents/GitHub/Aerohydro/resources/naca0012.dat','r')

    x,y = [],[]
    for line in inFile:
        data = line.split()
        x.append(float(data[0]))
        y.append(float(data[1]))
    x,y = np.array(x),np.array(y)
    inFile.close()
    return x,y
    
xp,yp = readGeometry()         # read of the geometry from a data file

# plotting the geometry
valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)

class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya                     # 1st end-point
        self.xb,self.yb = xb,yb                     # 2nd end-point
        self.xc,self.yc = (xa+xb)/2,(ya+yb)/2       # control point
        self.length = sqrt((xb-xa)**2+(yb-ya)**2)   # length of the panel
        
        # orientation of the panel
        if (xb-xa<=0.): self.beta = acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = pi+acos(-(yb-ya)/self.length)
        
        # location of the panel
        if (self.beta<=pi): self.loc = 'extrados'
        else: self.loc = 'intrados'
        
        self.sigma = 0.                             # source strength
        self.vt = 0.                                # tangential velocity
        self.Cp = 0.                                # pressure coefficient
        
        # function to descretize the geometry into panels
def definePanels(N,xp,yp):
    Np = len(xp)
    length = sum(np.sqrt((xp[1:Np]-xp[0:Np-1])**2+(yp[1:Np]-yp[0:Np-1])**2))
    length += sqrt((xp[0]-xp[Np-1])**2+(yp[0]-yp[Np-1])**2)
    lc = length/N
    
    x,y = np.empty(N,dtype=float),np.empty(N,dtype=float)
    x[0],y[0] = xp[0],yp[0]
    i = 0
    for j in range(N-1):
        xStart,yStart = x[j],y[j]
        xEnd,yEnd = xp[i+1],yp[i+1]
        d = sqrt((xStart-xEnd)**2+(yStart-yEnd)**2)
        if (d>=lc): d = 0
        elif (d<lc):
            while (i<Np-1):
                i += 1
                xStart,yStart = xEnd,yEnd
                if (i<Np-2): xEnd,yEnd = xp[i+1],yp[i+1]
                else: xEnd,yEnd = xp[0],yp[0]
                dd = sqrt((xStart-xEnd)**2+(yStart-yEnd)**2)
                if (d+dd<lc): d += dd
                else: break
        x[j+1] = xStart+(lc-d)/sqrt((xEnd-xStart)**2+(yEnd-yStart)**2)*(xEnd-xStart)
        y[j+1] = yStart+(lc-d)/sqrt((xEnd-xStart)**2+(yEnd-yStart)**2)*(yEnd-yStart)
                
    panel = np.empty(N,dtype=object)
    for i in range(N-1):
        panel[i] = Panel(x[i],y[i],x[i+1],y[i+1])
        panel[N-1] = Panel(x[N-1],y[N-1],x[0],y[0])
    
    return panel
    
N = 500                        # number of panels
panel = definePanels(N,xp,yp)  # discretization of the geometry into panels

# plotting the geometry with the panels
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
#plt.plot(append([p.xa for p in panel],panel[0].xa),\
#        append([p.ya for p in panel],panel[0].ya),\
#        linestyle='-',linewidth=1,\
#        marker='o',markersize=6,color='#CD2305')




class Freestream:
	def __init__(self,Uinf,alpha):
		self.Uinf = Uinf                   # velocity magnitude
		self.alpha = alpha*pi/180          # angle of attack (degrees --> radians)
		
# definition of the object freestream
Uinf = 1.0                               # freestream speed
alpha = 5.0                              # angle of attack (in degrees)
freestream = Freestream(Uinf,alpha)      # instantiation of the object freestream

# function to evaluate the integral Iij(zi)
def I(xci,yci,pj,dxdz,dydz):
	def func(s):
		return (+(xci-(pj.xa-sin(pj.beta)*s))*dxdz\
				+(yci-(pj.ya+cos(pj.beta)*s))*dydz)\
			   /((xci-(pj.xa-sin(pj.beta)*s))**2\
			   + (yci-(pj.ya+cos(pj.beta)*s))**2)
	return integrate.quad(lambda s:func(s),0.,pj.length)[0]
	
	# function to build the source matrix
def sourceMatrix(p):
	N = len(p)
	A = np.empty((N,N),dtype=float)
	np.fill_diagonal(A,0.5)
	for i in range(N):
		for j in range(N):
			if (i!=j):
				A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],+cos(p[i].beta),+sin(p[i].beta))
	return A

# function to build the vortex array
def vortexArray(p):
	N = len(p)
	B = np.zeros(N,dtype=float)
	for i in range(N):
		for j in range(N):
			if (j!=i):
				B[i] -= 0.5/pi*I(p[i].xc,p[i].yc,p[j],+sin(p[i].beta),-cos(p[i].beta))
	return B

# function to build the Kutta condition array
def kuttaArray(p):
    N = len(p)
    B = np.zeros(N+1,dtype=float)
    for j in range(N):
        if (j==0):
            B[j] = 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],sin(p[N-1].beta),-cos(p[N-1].beta))
        elif (j==N-1):
            B[j] = 0.5/pi*I(p[0].xc,p[0].yc,p[j],sin(p[0].beta),-cos(p[0].beta))
        else:
            B[j] = 0.5/pi*I(p[0].xc,p[0].yc,p[j],-sin(p[0].beta),+cos(p[0].beta))\
                + 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],-sin(p[N-1].beta),+cos(p[N-1].beta))
            B[N] -= 0.5/pi*I(p[0].xc,p[0].yc,p[j],+cos(p[0].beta),+sin(p[0].beta))\
                + 0.5/pi*I(p[N-1].xc,p[N-1].yc,p[j],+cos(p[N-1].beta),+sin(p[N-1].beta))
    return B

# function to assemble the global matrix
def buildMatrix(panel):
	N = len(panel)
	A = np.empty((N+1,N+1),dtype=float)
	AS = sourceMatrix(panel)
	BV = vortexArray(panel)
	BK = kuttaArray(panel)
	A[0:N,0:N],A[0:N,N],A[N,:] = AS[:,:],BV[:],BK[:]
	return A

# function to build the right hand-side of the linear system
def buildRHS(p,fs):
	N = len(p)
	B = np.zeros(N+1,dtype=float)
	for i in range(N):
		B[i] = -fs.Uinf*cos(fs.alpha-p[i].beta)
	B[N] = -fs.Uinf*(sin(fs.alpha-p[0].beta)+sin(fs.alpha-p[N-1].beta))
	return B
	
A = buildMatrix(panel)					# calculate the singularity matrix
B = buildRHS(panel,freestream)			# calculate the freestream RHS

# solve the linear system
var = np.linalg.solve(A,B)
for i in range(len(panel)):
	panel[i].sigma = var[i]
gamma = var[-1]

# function to calculate the tangential velocity at each control point
def getTangentVelocity(p,fs,gamma):
	N = len(p)
	A = np.zeros((N,N+1),dtype=float)
	for i in range(N):
		for j in range(N):
			if (i!=j):
				A[i,j] = 0.5/pi*I(p[i].xc,p[i].yc,p[j],-sin(p[i].beta),+cos(p[i].beta))
				A[i,N] -= 0.5/pi*I(p[i].xc,p[i].yc,p[j],+cos(p[i].beta),-sin(p[i].beta))
	B = fs.Uinf*np.sin([fs.alpha-pp.beta for pp in p])
	var = np.empty(N+1,dtype=float)
	var = np.append([pp.sigma for pp in p],gamma)
	vt = np.dot(A,var)+B
	for i in range(N):
		p[i].vt = vt[i]
		
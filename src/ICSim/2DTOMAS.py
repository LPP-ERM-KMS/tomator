import pickle
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
from ngsolve import *
from netgen.occ import *
from ngsolve.webgui import Draw
import netgen.geom2d as geom2d
from netgen.geom2d import CSG2d, Circle, Rectangle

fig, ax1 = plt.subplots()

freq = 11.44e6 #driving frequency
prefix = '11.44Mhz'

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

my_data = np.genfromtxt('data/'+prefix+'.csv', delimiter=',')
x_ = my_data[1:,0]
temp_ = my_data[1:,1]
cstemp = CubicSpline(x_,temp_)
x_range = np.arange(x_[0], x_[-1], 0.01)
tempsmooth = movingaverage(cstemp(x_range),2000)
dens_ = my_data[1:,2]
csdens = CubicSpline(x_,dens_)
x_range = np.arange(x_[0], x_[-1], 0.01)
denssmooth = movingaverage(csdens(x_range),2000)

#--------------------#
# Machine definition #
#--------------------#
ne = np.flip(denssmooth)#1e17 #density
gasn = 2e19
Te = np.flip(tempsmooth)
Ti = np.flip(tempsmooth)
I = 1600 #A
Power = 6000 #6kW IC
R0 = 0.780 #major radius
Ra = 0.260 #minor radius

MAXH=0.02
order_mesh = 3

try: 
    with open('pickled/mesh.pkl', 'rb') as file:
        mesh = pickle.load(file)
except:
    #--------------------#
    #        Mesh        #
    #--------------------#
    geo = CSG2d()
    wall = Circle( center=(0,0), radius=0.26+0.78, mat="Plasma", bc="Wall" )
    hole = Circle( center=(0,0), radius=0.78-0.26, mat="copper", bc="Wall" )
    antenna = Rectangle( pmin=(0.21+0.78,-0.044), pmax=(0.21+0.78+0.0035,0.044), mat="copper", bc="Antenna" )
    #NoLeftSide = Rectangle( pmin=(-0.27-0.78,-0.27), pmax=(0,0.27), mat="copper", bc="Wall" )
    System = wall-hole-antenna#-NoLeftSide
    geo.Add(System)

    with TaskManager():
        mesh = Mesh(geo.GenerateMesh (maxh=MAXH))
    mesh.Curve(5)
    with open('pickled/mesh.pkl', 'wb') as file:
        pickle.dump(mesh, file)

from pyRFplasma.system import System
from pyRFplasma.solve import Solve
from pyRFplasma.constants import Constants

TOMAS = System({"e":1,"H2":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMAS.Epsilon2D(MAXH,temperature="CXD")

solution = Solve(TOMAS)
solution.GetSolution()

# Find scaling
P = solution.PowerDeposition2D()

resolution = 301
angle = np.pi*10/180
TP = []
R = np.linspace(R0-Ra+0.01,R0+Ra-0.01,resolution)
for r in R:
    TP.append(P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real)
PowerScalingFactor = 6000/sum(TP)
TP = [TP[i]*PowerScalingFactor for i,j in enumerate(TP)]
plt.plot(R,TP)
plt.show()

# Compute ion heating
TOMASH2 = System({"H2":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMASH2.Epsilon2D(MAXH,temperature="CXD")
H2diel = TOMASH2.eps
P = solution.PowerDeposition2D(H2diel)
resolution = 301
angle = np.pi*10/180
H2P = []
R = np.linspace(R0-Ra+0.01,R0+Ra-0.01,resolution)
for r in R:
    H2P.append(P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real)
H2P = [H2P[i]*PowerScalingFactor for i,j in enumerate(H2P)]
plt.plot(R,H2P)
plt.show()
# Compute electron heating
TOMASe = System({"e":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMASe.Epsilon2D(MAXH,temperature="CXD")
ediel = TOMASe.eps
P = solution.PowerDeposition2D(ediel)
resolution = 301
angle = np.pi*10/180
eP = []
R = np.linspace(R0-Ra+0.01,R0+Ra-0.01,resolution)
for r in R:
    eP.append(P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real)
TP = [eP[i]*PowerScalingFactor for i,j in enumerate(eP)]
plt.plot(R,eP)
plt.show()

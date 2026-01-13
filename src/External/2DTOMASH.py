import csv
import pickle
import numpy as np
from numpy import genfromtxt
from scipy.interpolate import CubicSpline
from ngsolve import *
from netgen.occ import *
from ngsolve.webgui import Draw
import netgen.geom2d as geom2d
from netgen.geom2d import CSG2d, Circle, Rectangle

with open('/tmp/DensAndTemp.csv') as f:
    info = str(f.readline().strip('\n'))
freqstr = info.split(",")[-1]
freq = float(info.split(",")[-1])*1e6

prefix = freqstr+'Mhz'

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

Profiles = np.genfromtxt('/tmp/DensAndTemp.csv', delimiter=',',skip_header=1)

#--------------------#
# Machine definition #
#--------------------#
gasn = {"H":1e18,"H2":1e19}
r_ = Profiles[1:,0]*1e-2 #in cm
ne = Profiles[1:,1]*1e6 #in cm^-3
Te = Profiles[1:,2]
ni = Profiles[1:,3]*1e6 #in cm^-3
Ti = Profiles[1:,4]

I = 1600 #A
Power = 6000 #6kW IC
R0 = 0.780 #major radius
Ra = 0.260 #minor radius

MAXH=0.01
order_mesh = 3

try: 
    with open('/tmp/mesh.pkl', 'rb') as file:
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
    with open('/tmp/mesh.pkl', 'wb') as file:
        pickle.dump(mesh, file)

#pyRFplasma is my custom library
from pyRFplasma.system import System 
from pyRFplasma.solve import Solve
from pyRFplasma.constants import Constants

TOMAS = System({"e":1,"H":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn,orientation="horizontal",ni=ni)
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

# Compute ion heating
TOMASH2 = System({"H":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMASH2.Epsilon2D(MAXH,temperature="CXD")
H2diel = TOMASH2.eps
P = solution.PowerDeposition2D(H2diel)
angle = np.pi*10/180
HP = []
R = np.linspace(R0-Ra+0.01,R0+Ra-0.01,resolution)
for r in R:
    HP.append(P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real*PowerScalingFactor)

# Compute electron heating
TOMASe = System({"e":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMASe.Epsilon2D(MAXH,temperature="CXD")
ediel = TOMASe.eps
P = solution.PowerDeposition2D(ediel)
angle = np.pi*10/180
eP = []
R = np.linspace(R0-Ra+0.01,R0+Ra-0.01,resolution)
for r in R:
    eP.append(P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real*PowerScalingFactor)

with open('/tmp/PowerDeposition.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    spamwriter.writerow(['eP','HP'])
    for i in range(len(TP)):
        spamwriter.writerow([eP[i],HP[i]])

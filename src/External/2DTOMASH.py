import csv
import time
import pickle
import numpy as np
from numpy import genfromtxt
from scipy.interpolate import CubicSpline
from ngsolve import *
from netgen.occ import *
from ngsolve.webgui import Draw
import netgen.geom2d as geom2d
from netgen.geom2d import CSG2d, Circle, Rectangle
import logging

with open('/tmp/DensAndTemp.csv') as f:
    info = str(f.readline().strip('\n'))
freqstr = info.split(",")[-1]
infostring = info.split(",")
freq = float(info.split(",")[-1])

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')

Profiles = np.genfromtxt('/tmp/DensAndTemp.csv', delimiter=',',skip_header=1)

logger = logging.getLogger(__name__)
logging.basicConfig(filename='/tmp/pyrfplasma.log', level=logging.INFO)
logger.info('Started')

#--------------------#
# Machine definition #
#--------------------#
r_ = Profiles[1:,infostring.index("Ra")]*1e-2 #in cm
ne = Profiles[1:,infostring.index("Ne")]*1e6 #in cm^-3
Te = Profiles[1:,infostring.index("Te")]
ni = Profiles[1:,infostring.index("nHi")]*1e6 #in cm^-3
Ti = Profiles[1:,infostring.index("THi")]
gasn = {"H":Profiles[1:,infostring.index("nH")]*1e6,"H2":Profiles[1:,infostring.index("nH2")]*1e6} #moet nog gelezen worden

I = 1600 #A
Power = 5000 #5kW IC
R0 = 0.780 #major radius
Ra = 0.260 #minor radius

MAXH=0.02
meshscalefactor=20 #determined through multiple simulations to be minimal with high accuracy
order_mesh = 2

tic = time.time()
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
    logger.info(f'refining mesh around antenna')
    for el in mesh.Elements():
        for v in el.vertices:
            if (mesh[v].point[0]>0) and (abs(mesh[v].point[1])<0.78):
                mesh.SetRefinementFlag(el, True)
            else:
                mesh.SetRefinementFlag(el, False)
    mesh.Refine()
    logger.info(f'refining mesh around antenna')
    for el in mesh.Elements():
        for v in el.vertices:
            if (mesh[v].point[0]>0) and (abs(mesh[v].point[1])<0.3):
                mesh.SetRefinementFlag(el, True)
            else:
                mesh.SetRefinementFlag(el, False)
    mesh.Refine()
toc = time.time()
logger.info(f'loading/creating mesh took {toc-tic}s')

#pyRFplasma is my custom library
from pyRFplasma.system import System 
from pyRFplasma.solve import Solve
from pyRFplasma.constants import Constants

tic = time.time()
logger.info(f'creating dielectric')
TOMAS = System({"e":1,"H":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn,orientation="horizontal",ni=ni)
TOMAS.Epsilon2D(MAXH/meshscalefactor,temperature="CXD")
toc = time.time()
logger.info(f'creating dielectric took {toc-tic}s')

solution = Solve(TOMAS)
solution.GetSolution(order_mesh=order_mesh)

# Find scaling
P = solution.PowerDeposition2D()

resolution = 301
nAngles = 10
angles = np.pi*np.linspace(19/180,21/180,nAngles) #around LP
TP = np.zeros(resolution)
R = np.linspace(R0-Ra,R0+Ra,resolution)
for angle in angles:
    for i,r in enumerate(R):
        TP[i] += (P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real)/nAngles
PowerScalingFactor = (Power*26*2)/(sum(TP)*resolution)
TP = [TP[i]*PowerScalingFactor for i,j in enumerate(TP)]

# Compute ion heating
tic = time.time()
TOMASH2 = System({"H":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMASH2.Epsilon2D(MAXH/meshscalefactor,temperature="CXD")
H2diel = TOMASH2.eps
P = solution.PowerDeposition2D(H2diel)
HP = np.zeros(resolution)
for angle in angles:
    for i,r in enumerate(R):
        HP[i] += (P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real)/nAngles
toc = time.time()
logger.info(f'ion power deposition calculation total took {toc-tic}s')

# Compute electron heating
tic = time.time()
TOMASe = System({"e":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn)
TOMASe.Epsilon2D(MAXH/meshscalefactor,temperature="CXD")
ediel = TOMASe.eps
P = solution.PowerDeposition2D(ediel)
angle = np.pi*10/180
eP = np.zeros(resolution)
for angle in angles:
    for i,r in enumerate(R):
        eP[i] += (P(mesh(r*np.cos(angle),r*np.sin(angle)))[0].real)/nAngles
toc = time.time()
logger.info(f'electron power deposition calculation total took {toc-tic}s')

with open('/tmp/PowerDeposition.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    spamwriter.writerow(['eP','HP'])
    for i in range(len(TP)):
        spamwriter.writerow([eP[i],HP[i]])

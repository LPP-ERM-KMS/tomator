import csv
import time
import pickle
import numpy as np
from scipy.interpolate import CubicSpline
import ngsolve as ngs
from netgen.occ import *
from netgen.meshing import *
from pyRFplasma.helpeqs import Gaussian
import logging

with open('/tmp/DensAndTemp.csv') as f:
    info = str(f.readline().strip('\n'))
infostring = info.split(",")
freq = float(info.split(",")[-2])
Power = float(info.split(",")[-1])*1e3

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

TEST=False
TEMP="CXD"
logger = logging.getLogger(__name__)
logging.basicConfig(filename='/tmp/pyrfplasma.log', level=logging.INFO)
logger.info('Started')

tic = time.time()

I = 1600 #A
R0 = 0.780 #major radius
Ra = 0.260 #minor radius
Rant = 0.210 #antenna
MAXH=0.0008
resolution = 1001 #for tomator
meshscale = 2
order_mesh = 3
maxn = 40
neutral_temperature = 55
mode_numbers = [i for i in range(-maxn,maxn)]
R = np.linspace(R0-Ra,R0+Rant,int(2*Ra/MAXH))
R_ = np.linspace(R0-Ra,R0+Ra,resolution)
# split at antenna
#R1 = R_[R_ <= R0+Rant]
#R2 = R_[R_ >= R0+Rant]

try: 
    with open('/tmp/1Dmesh.pickle', 'rb') as file:
        mesh = pickle.load(file)
except:
    nnodes = R.size
    part = [R[0], R[-1]]

    n = []
    for i in range(len(part) - 1):
        n.append(nnodes)

    unit_cell = Mesh(dim = 1)

    diff = 0
    nodes = np.empty(0)
    for i in range(len(n)):
        a = np.linspace(part[i] + diff, part[i + 1], n[i])
        nodes = np.append(nodes, a)
        diff = np.abs(nodes[len(nodes) - 1] - nodes[len(nodes) - 2])
    # nodes = np.delete(nodes, 0)

    pnums = []
    for i in range(sum(n)):
        pnums.append(unit_cell.Add(MeshPoint(Pnt(nodes[i], 0, 0))))

    idx1 = unit_cell.AddRegion("plasma", dim=1)

    for i in range(len(nodes)-1):
        if nodes[i] <= part[1]:
            unit_cell.Add (Element1D([pnums[i], pnums[i+1]], index = idx1))

    id5 = unit_cell.AddRegion("l", dim=0)
    id6 = unit_cell.AddRegion("r", dim=0)
    id7 = unit_cell.AddRegion("ant", dim=0)

    unit_cell.Add (Element0D(pnums[0] , index=id5))
    unit_cell.Add (Element0D(pnums[-1], index=id6))
    unit_cell.Add (Element0D(pnums[len(R)-10], index=id7))

    #%% Meshing
    mesh = ngs.Mesh(unit_cell)
    meshw = open('/tmp/1Dmesh.pickle', 'wb')
    pickle.dump(mesh,meshw)



toc = time.time()
logger.info(f'loading/creating mesh took {toc-tic}s')

tic = time.time()

from pyRFplasma.system import System
from pyRFplasma.solve import Solve
from pyRFplasma.constants import Constants


TOMAS = System({"e":1,"H":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn,neutral_temperature=neutral_temperature,mode_numbers=maxn)
TOMAS.Epsilon1D(MAXH/meshscale,temperature=TEMP)

toc = time.time()
logger.info(f'Creating dielectric took {toc-tic}s')

tic = time.time()

for i,mode_number in enumerate(mode_numbers):
    solution = Solve(TOMAS)
    solution.GetSolution1D(order_mesh=order_mesh,mode_number=mode_number)
    gfx,gfy,gfz = solution.result

    if i == 0:
        Ex = gfx
        Ey = gfy
        Ez = gfz
    else:
        Ex += gfx
        Ey += gfy
        Ez += gfz

toc = time.time()
logger.info(f'Computing electric fields took {toc-tic}s')

tic = time.time()

Ptot = solution.PowerDeposition1D(E=[Ex,Ey,Ez])
TP = np.zeros(len(R_))
for i,r in enumerate(R_):
    if r < R0+Rant:
        TP[i] = Ptot(mesh(r))[0].real
    else:
        TP[i] = 1e-12 #not zero for stability reasons
TP = [tp if tp>= 0 else 0 for tp in TP]

toc = time.time()
logger.info(f'Computing power deposition took {toc-tic}s')

PowerScalingFactor = (Power)/(sum(TP)) #Watt per meshpoint
TOMASe = System({"e":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn,neutral_temperature=neutral_temperature)
TOMASe.Epsilon1D(MAXH/meshscale,temperature=TEMP)
Pe = solution.PowerDeposition1D(customdielectric=TOMASe.eps,E=[Ex,Ey,Ez])
eP = np.zeros(len(R_))
for i,r in enumerate(R_):
    if r < R0+Rant:
        eP[i] = Pe(mesh(r))[0].real*PowerScalingFactor
    else:
        eP[i] = 1e-12
eP = [ep if ep>= 0 else 0 for ep in eP]

TOMASi = System({"H":1},I,freq,Power,ne,Ti,Te,R0,Ra,mesh,gasnd=gasn,neutral_temperature=neutral_temperature)
TOMASi.Epsilon1D(MAXH/meshscale,temperature=TEMP)
Pi = solution.PowerDeposition1D(customdielectric=TOMASi.eps,E=[Ex,Ey,Ez])
HP = np.zeros(len(R_))
for i,r in enumerate(R_):
    if r < R0+Rant:
        HP[i] = Pi(mesh(r))[0].real*PowerScalingFactor
    else:
        HP[i] = 1e-12
HP = [hp if hp>= 0 else 0 for hp in HP]

with open('/tmp/PowerDeposition.csv', 'w', newline='') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',')
    spamwriter.writerow(['eP','HP'])
    for i in range(len(TP)):
        spamwriter.writerow([eP[i],HP[i]])

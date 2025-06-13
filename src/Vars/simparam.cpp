#include "simparam.h"

int nmeshp = 160;    // index of the center grid point for output to logfile
int cc = 80;    // index of the center grid point for output to logfile

// MAIN DISCHARGE PARAMETERS:
// magnetic field
double Bt = 1.5;                                                 // T, toroidal magnetic field on axis
double Bv = 0.0075;                                              // T, vertical magnetic field
double Bh = 0.0001;                                              // T, radial magnetic field
                                                                 // toroidal machine geometry
double R = 88.0;                                                 // cm, Vessel major radius
double a = 28.0;                                                 // c, Vessel minor radius from center to wall
double b = 75.0;                                                 // c, Vessel vertical radius from center to wall
double lHFS = R - 0.914 * a;                                     // cm, radial location of limiter at HFS ;
double lLFS = R + 0.914 * a;                                     // c, radial location of limiter at LFS ;
double nlimiters = 6.0;                                          // number of poloidal limiters ;
double Vpl = 2.0 * b * pi * (pow(R + a, 2.0) - pow(R - a, 2.0)); // simulation volume (here calculated, but can be specified...)
                                                                 // neutral pressure
double pHe = 1.68E-04;                                           // He pressure at outer grid points
double pH2 = 4.32E-05;                                           // H2 pressure at outer grid points

// COUPLED (RF) POWER
double Prf = 4e2;            // kW, launched or coupled power level
bool bselfcol = false;        // include self collisions in case of sollisional absorption (should always be false, this bool can be deleted later)
double freq = 82.7e3;        // MHz, RF frequency
double dtpramp = 0.0;        // ramp up power
                             //// ECH
double Rdep = 89.4;          // the radial location where the power is coupled (resonance layer)
bool bgray = false;          // either use gray lookup table for power coupling: dependent on ne and Te
bool bram = false;           // or use Tulchhi's lookup table for power coupling: dependent on ne, Te and nue
bool bfixpowerfrac = false;  // couple percentage of Prf
double pecabs0 = 0.1;        // fraction of Prf that is coupled, used in 'bfixpowerfrac' and 'bnefix'
bool bnefix = true;          // tune coupled EC power to obtain the central powerdensity
int ic = 86;                 // grid point where the ne should be fixed
double necfix = 1.0584e13;   // at 90.1m, fitted manual, 86
double Pini = 3.0;           // initial P coefficient
double tauP = 0.00005;       // P tuning parameter
double widthech = a / 15.0;  // width of resonance : R * exp(-pow((R - Rdep) / (widthech, 2.0)); AUG: a / 15.0, TCV: a / 30.0
double echbackground = 1e-7; //                          + echbackground * R * exp(-pow((R - Rdep) / a, 2.0))
bool bTOMAS = false;         // Johans bunch
double harmonic = 2.0;       // number of harmonic frequency
double muw = 0.01;
double Rdep1 = 86.1875;
double PRdep1 = 0.0; // % O-wave
double Rdep2 = 97.65;
double PRdep2 = 1.0; // % X-wave
double Rdep3 = 72.4325;
double PRdep3 = 0.0; // % B-wave coupling efficiency
double Rdep4 = 72.4325;
double PRdep4 = 0.0;     // % UHR coupling percentage

bool bmanuel = false;        // 
std::string smanuel = "";  //

bool bICWC = false; // use ICWC coupling

                         //// ICH
bool bkipt = false;      // use kipt RF module
bool bantlr = false;  // calculate coupling using antenna resistance, works only for kipt module
double avlr = 0.23;      // if bantlr == true, then this is the antenna vacuum resistance
double alphalaw = 0.002; // else use alphalaw
double HtoHD = 1.0;      // not used now (only for the KIPT RF module)
        //// ICH
bool blhr = false;      // use kipt RF module
double Rant = 220.0;       // radial position of the antenna
double widthlhr = a / 15.0;  // width of resonance : R * exp(-pow((R - Rdep) / (widthech, 2.0)); AUG: a / 15.0, TCV: a / 30.0
double lhrbackground = 1e-7; //                          + echbackground * R * exp(-pow((R - Rdep) / a, 2.0))
double fracpne = 1.0;     // fraction coupled to electrons
double fraclhr = 0.25;     // fraction coupled per lhr (forward and backward)
                         //// OTHER
bool bproptone = false;
bool fixedTe = false;
bool bnopower = false;

bool bDfix = false;        // use transport scaling of TCV paper
bool bDfix_neutr = false;
double Dfix = 10000.0;         // cm2/s
double Dfix_neut = 1000000.0;
bool bDbohm = false;       // use transport scaling of TCV paper
bool bDscaling = true;    // use transport scaling of TCV paper
bool bVscaling = true;     // use transport scaling of TCV paper

// TRANSPORT PARAMETERS
bool bVfix = false;        // use transport scaling of TCV paper
double Vfix = 100.0;       // cm/s
// bool bscaling = false;     // use transport scaling of TCV paper
int veq = 8;               // which equation to use for convection scaling (8 or 9)
double Dfact = 1.0;        // diffusion prefactor
double Vfact = 5.0;        // convection prefactor
bool btunedv = true;       // then tune D and V
bool btunevleft = true;    // true for Bv scan (high density), false for P scan (low density)
int il = 73;               // left side grid point for which D and V is tuned to match density
double nelfix = 3.1898e12; // left side density for which D and V is tuned
int ir = 114;              // right side grid point for which D and V is tuned to match density
double nerfix = 6.4301e12; // right side density for which D and V is tuned
double Vini = -2.0;        // initial V coefficient
double tauV = 0.01;        // V tuning parameter
double Dini = 1.0;         // initial D coefficient
double tauD = 0.0005;      // D tuning parameter

// PHYSICS TO Include
bool bH = true;
bool bH2 = true;
bool bHe = true;
bool bADAS = false;
bool bion = true;
bool bcx = true;
bool belas = true;
bool bcoulomb = true;
bool bimpur = false;
bool btranspions = true;
bool btranspneut = true;
bool bedge = true;
bool bpol = true;
bool vdrift = true;
bool bcoll = true;

// INITIAL CONDITIONS
double rmaxini = 89.4;            // initial density distribution: location of maximum
double widthini = 0.1 * a;        // initial density distribution: width of gaussian
double nebackgroundl = 1e-5;      // initial hfs background density
double nebackgroundr = 1e-3;      // initial lfs background density
double Ta0 = 0.02587;             // initial atom and molecule temperature
double Te0 = 3.0;                 // initial charged particle temperature
double nevac = 1.0e0;             // vacuum density (if any density goes lower, it is no longer accurately traced)
double nH0 = 2.0e6;               // initial H atom density
double nHi0 = 1.0e6;              // initial Hi ion density
double nH2i0 = 1.0e6;             // initial H2i ion density
double nH3i0 = 2.0e6;             // initial H3i ion density
double nHeII0 = 1.0584e13 * 0.95; // initial HeII ion density
double nHeIII0 = 2.0e6;           // initial HeIII ion density
double nCII0 = 0.1e10;            // initial CII ion density
double nCIII0 = 0.1e1;            // initial CIII ion density
double nCIV0 = 0.1e1;             // initial CIV ion density
double nCV0 = 0.1e1;              // initial CV ion density
                                  // the initial electron density is calculated from the above ion densities

// EDGE CONDITIONS
double RH = 0.5;  // H atom reflection coefficient
double REH = 0.9; // H atom energy reflection coefficient

double gEd = 2.5;  // ;3! ( 3/2 + d(log(Te))/d(log(Te)) ) E=3/2*n*T
double gEv = 1.5;  // ;1! ( 3/2 to have equal n and 3/2nT transport...)
double gEdn = 2.0; // ;3! ( 3/2 + d(log(sqrt(Te)))/d(log(Te)) )
double gEe = 2.5;  // ;3!
bool fixBCs = false;
double lam_dec_length_ions = 2.0;
double lam_dec_length_energy = 1.0;


// SIMULATION CONTROL SETTINGS
// input file
bool bfinput = false; // start calculation from saved "input.txt" // leave no empty spaces at start or ending of file.
std::string sinputfile = "noinput";
// time step (note machine accuracy is 2.22045e-16)
double t0 = 0.0;           // simulation start time
double tmainend = 1.0;     // simulation end time
double accur = 0.05;       // state cannot change more than accur each time step
double dtmax = 2.0e-6;     // maximum time step
double dtmin = 0.25e-7;    // minimum time step (works only after initial transition phase starting from dtinit)
double dtinit = 1.0e-11;   // minimum time step (works only after initial transition phase starting from dtinit)
                           // time step for RF coupling (usefull is subroutine is slow)
bool bupdateRFstep = true; // update RF every N time steps of every dtRF
int updateRF = 1;          // Update RF deposition profile every N interations
double dtRF = 1.0e-7;      // Update RF deposition profile every dtRF times
bool dtRFvar = false;      // change dtRF based step differences
double dtRFmax = 5e-6;
double dtRFmin = 0.1e-7;
bool dtRFconv = true;
// advanced time step settings
bool dtsmooth = true;           // reduce time step to 10% if discontinuity in time is seen  (for ICRF)
double shokparam = 2;           // reduce time step to 10% if shok (?) is seen  (for ICRF)
double maxtstepincrement = 1.5; // value > 1: max increase of dt from one tstep to another.
bool btendvar = false;           // end simulation if converged
double convcrit = 100.0;
double convsavetime = 0.01;
bool baccurvar = true; // increase requisted accuracy of certain convergence is reached
double accurcrit = 500;
double minaccur = 0.05;
// output parameters
int Nlog = 500; // output to logfile every Nlog iterations
std::string sOutputfolder = "";
int Nloopsave = 5000;   // save to file every N iterations
bool bOutdt = true;     // save additionally at dtsave (if dtsave < N*<dt>)
double dtsave = 1.0e-3; // dtsave needs to be larger than dtmax!
                        // solver parameters
double solvertolerance = 1e-10;

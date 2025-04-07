// simulation parameters
// AUG4
#ifndef SIMPARAM_H
#define SIMPARAM_H

#define NMESHP 161 // number of radial grid points, recompile code to change value
//#define ISCAN 1    // number of parallel simulations in scan

#include <cmath>
#include <cstring>
#include <string>

#include "constants.h"

extern int nmeshp; // number of radial grid points

extern int cc; // index of the center grid point ()

extern double Bt; // TEXTOR, magnetic field on axis
extern double Bv; // vertical magnetic field, losses along magnetic field line (idealy Bv determines convection)
extern double Bh; // vertical magnetic field, losses along magnetic field line

extern double freq;     // MHz, RF frequency
extern double harmonic; // number of harmonic frequency

extern double a;         // AUG, vessel minor radius from center to wall
extern double b;         // AUG, vessel minor radius from center to wall
extern double R;         // AUG, vessel major radius
extern double lHFS;      // radial location of limiter at HFS ;
extern double lLFS;      // radial location of limiter at LFS ;
extern double nlimiters; // number of poloidal limiters ;
extern double Vpl;       // 2.0*pow(pi,2.0)*pow(a,2.0)*R;	// Plasma volume [cm3]

extern double Prf; // paper v6
extern double dtpramp;

extern bool bkipt;      // use kipt RF module
extern bool bantlr;  // calculate coupling using antenna resistance, works only for kipt module
extern double avlr;     // if bantlr = , then vacuum resistance
extern double alphalaw; // else use alphalaw
extern bool bselfcol;    // else use alphalaw
extern bool bgray;
extern double Rdep;     // 301 grid points
extern double widthech; // width of resonance : R * exp(-pow((R - Rdep) / (widthech, 2.0)); AUG: a / 15.0, TCV: a / 30.0
extern double echbackground;

extern bool bmanuel;
extern std::string smanuel;

extern bool bICWC; // use ICWC coupling

extern bool bram;

extern bool bTOMAS; // extern    double pecabs0[] = {0.0020, 0.0020, 0.0020, 0.0020, 0.0020, 0.0020, 0.0020}; // percentage of power coupled
extern double muw;
extern double Rdep1;
extern double PRdep1; // % O-wave
extern double Rdep2;
extern double PRdep2; // % X-wave
extern double Rdep3;
extern double PRdep3; // % B-wave coupling efficiency
extern double Rdep4;
extern double PRdep4; // % UHR coupling percentage

extern bool blhr;      // use kipt RF module
extern double Rant;       // radial position of the antenna
extern double fracpne;     // else use alphalaw
extern double fraclhr;     // else use alphalaw
extern double widthlhr; // width of resonance : R * exp(-pow((R - Rdep) / (widthech, 2.0)); AUG: a / 15.0, TCV: a / 30.0
extern double lhrbackground;

extern bool bfixpowerfrac;
extern bool bnefix;
extern double necfix; // at 90.1m, fitted manual, 86
extern double nelfix; // at 0.856m 73
extern double nerfix; // at 1.039m

extern int ic; // not used here
extern int il; // not used here
extern int ir; // not used here
//

extern double Vini; // mbar
extern double Dini; // mbar
extern double Pini; // mbar

extern double tauV; // mbar
extern double tauD; // mbar
extern double tauP; // mbar

extern double Dfact; // 0.3Bv%+1
// extern    double Vfact5}; // eq9
extern double Vfact;   // eq8
extern double pecabs0; // 5Bv%^2+2

extern bool bDfix;   // use scalings instead of fitting
extern double Dfix; // 0.3Bv%+1
extern bool bDbohm;  // 
extern bool bDscaling;   


extern bool bVfix;   // use scalings instead of fitting
extern double Vfix; // 0.3Bv%+1
extern bool bVscaling;
// extern bool bscaling;   // use scalings instead of fitting
extern int veq;         // which equation to use for convection scaling (8 or 9)
extern bool btunedv;    // then tune D and V
extern bool btunevleft; // true for Bv scan (high density), false for P scan (low density)

extern bool bcentral;
extern int centerval;

extern bool bproptone;
extern bool bnopower;

extern double pHe;   // mbar
extern double pH2;   // mbar
extern double HtoHD; // only for the KIPT RF module, not functional anymore at the moment

extern double rmaxini;       // initial density distribution: location of maximum
extern double widthini;      // initial density distribution: width of gaussian
extern double nebackgroundl; // initial hfs background density
extern double nebackgroundr; // initial lfs background density

extern double Ta0; //
extern double Te0;
// extern    double K = 2.1; // kappa parameter to rescale temperature dependency in transport  Tmaxw = Tkappa (K-3/2)/K

extern double nevac;

extern double nH0;
extern double nHi0;
extern double nH2i0;
extern double nH3i0;
extern double nHeII0; // 4.03e5
// extern    double nHeII0[] = {2036473068884.78*0.5,	9.9574E+12*0.5,    9.7360E+12*0.5,    5.0220E+12*0.5,    2.7836E+12*0.5,    2.6233E+12*0.5} ;
extern double nHeIII0;
extern double nCII0;
extern double nCIII0;
extern double nCIV0;
extern double nCV0;

extern double RH; // reflection coefficient
// extern    double RHi = 0.998; // ion recycle coefficient
// extern    double RH2 = 1.0; // reflection coefficient
// extern    double RHe = 1.0; // reflection coefficient
// extern    double RHei = 0.998; // ion recycle coefficient
extern double REH; // energy reflection coefficient
// extern    double REH2 = 0.9;  // energy reflection coefficient
// extern    double REHe = 0.9; //check

extern double gEd;  // ;3! ( 3/2 + d(log(Te))/d(log(Te)) ) E=3/2*n*T
extern double gEv;  // ;1! ( 3/2 to have equal n and 3/2nT transport...)
extern double gEdn; // ;3! ( 3/2 + d(log(sqrt(Te)))/d(log(Te)) )
extern double gEe;  // ;3!

// extern    extern    bool pcst   ; // if true, ncst should be false
// extern    extern    bool ncst  ; // if true, pcst should be false
// extern    extern    bool quick   ; // if true, others should be false
// extern    extern    bool bpump   ;
// extern    extern    bool binj   ;
// double PSH2 = Vpl/0.01; // [cm3/s]
// double PSHe = Vpl/0.01; // [cm3/s]
// extern    double taupumpH2 = Vpl/PSH2; // [s]
// extern    double taupumpHe = Vpl/PSHe; // [s]

extern bool bupdateRFstep;
extern int updateRF; // INT_MAX; // Update RF deposition profile every 300 interations
extern double dtRF;  // Update RF deposition profile every dtRF times

extern double accur; // must be smaller than << 1
extern double dtmax;
// extern    double dtmin = 1.1*2.22045e-16; // Machine accuracy!
extern double dtmin; // 0.3e-7
extern double dtinit;
extern double t0;
extern double tmainend;
extern double tloopendinit;

extern std::string sOutputfolder;

extern int Nlog;
extern int Nloopsave;
extern bool bOutdt;
extern double dtsave; // needs to be larger than dtmax!

extern bool bH;
extern bool bH2;
extern bool bHe;
extern bool bion;
extern bool bcx;
extern bool belas;
extern bool bcoulomb;
extern bool bimpur;
extern bool btranspions;
extern bool btranspneut;
extern bool bedge;
extern bool bpol;
extern bool vdrift;
extern bool bcoll;
extern bool fixedTe;

extern bool bfinput; // start calculation from saved "input.txt" // leave no empty spaces at start or ending of file.
extern std::string sinputfile;

extern bool btendvar;
extern double convcrit;
extern double convsavetime;

extern bool baccurvar;
extern double accurcrit;
extern double minaccur;

extern double solvertolerance;

extern bool dtRFvar;
extern double dtRFmax;
extern double dtRFmin;
extern bool dtRFconv;

extern bool dtsmooth;
extern double maxtstepincrement; // Value > 1 to smooth increment from one tstep to another.
extern double shokparam;

#endif

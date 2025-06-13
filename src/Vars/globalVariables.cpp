#include "globalVariables.h"

DENS nr, nr_m1, nr_p1;
TEMP Tr;
ENER Er, Er_m1, Er_p1;
dDENS dnr;
dENER dEr;
dDENS dnr_cn;
dENER dEr_cn;
CRATE colrate;
CRATE colrateRF;
POWER Power;

double Te_RF[NMESHP];
double ne_RF[NMESHP];
double nHi_RF[NMESHP];
double nH2i_RF[NMESHP];
double nH3i_RF[NMESHP];
double nHeII_RF[NMESHP];
double nHeIII_RF[NMESHP];
double dne_RF[NMESHP];
double dnHi_RF[NMESHP];
double dnH2i_RF[NMESHP];
double dnH3i_RF[NMESHP];
double dnHeII_RF[NMESHP];
double dnHeIII_RF[NMESHP];
double PRFe_id[NMESHP];
double PRFHi_id[NMESHP];
double PRFH2i_id[NMESHP];
double PRFH3i_id[NMESHP];
double PRFHeII_id[NMESHP];
double PRFHeIII_id[NMESHP];
double PRFe_array_stat[NMESHP];
double PRFHi_array_stat[NMESHP];
double PRFH2i_array_stat[NMESHP];
double PRFH3i_array_stat[NMESHP];
double PRFHeII_array_stat[NMESHP];
double PRFHeIII_array_stat[NMESHP];

double tmain = t0;
double tRF = t0;
double dtcalc = 0.0;

double dt = 0;
double dt2 = 0;
double dt3 = 0;
double dtnew = dtinit;
double drval = 0.0;
double drval2 = 0.0;
double drval3 = 0.0;
double drval4 = 0.0;
double drmax = 0;
double drmax2 = 0;
double drmax3 = 0;
double drmax4 = 0;

int cnt_save = 0;
int cnt_loop = 1;
double dtout = 0;
int Nit = 0;
int NitOut = Nit + Nloopsave;

double nefact;

double pecabs;
double nHeI0; // constant op de uiterste gridpunten
double nH20;  // constant op de uiterste gridpunten
double nCI0;

double Pabs1[NMESHP];
double Pabs2[NMESHP];
double Pabs3[NMESHP];
double Pabs4[NMESHP];

double Dionh[NMESHP];
double Dion[NMESHP];
double Vionh[NMESHP];
double Vion[NMESHP];
double mfpi = 0.0;
double gri = 0.0;

double nui = 0.0;

double coef1 = 0, coef2 = 0, coef3 = 0;

Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver2;
Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver3;
Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver4;

int mytid = 0;

double PIerrorD = Dini;
double Dfsave = 10.0;
double PIerrorV = Vini;
double Vfsave = 100.0;
double PIerrorP = Pini;

double alphaval;

double PRFe_array[NMESHP];
double PRFHi_array[NMESHP];
double PRFH2i_array[NMESHP];
double PRFH3i_array[NMESHP];
double PRFHeII_array[NMESHP];
double PRFHeIII_array[NMESHP];

double AverageTimer = 0.0;
double AvTimeStep = 0.0;
double StepTimer[9] = {0.0};
double ColTimer = 0.0;
double ParTimer = 0.0;
double RFTimer = 0.0;
double TimeStepTimer = 0.0;

double oldtstep = 0.0;
double olddaccur = 0.0;

DENS n_save = nr;
ENER E_save = Er;
int convsave = (int)(convsavetime / dtsave);
int convsavemin = (int)(convsavetime * 0.1 / dtsave);
int convsavemax = convsave;
int convsavecount = 0;
double tdifmax = 0.0;
double tdifne = 0.0;
double tdifTe = 0.0;
double minstopcrit = 5e6;
double tsum = 0.0;

DENS nRF_save;
ENER ERF_save;
double Rfdif = 0;
double Rfdifmax = 0;
double dtRFtemp = 0.0;
double dtRFnunc = 0.0;

DENS nr1 = nr;
DENS nr1bis;
ENER Er1 = Er;
ENER Er1bis;

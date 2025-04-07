#ifndef GLOBALVARIABLES_H
#define GLOBALVARIABLES_H

#include "../Eigen/Sparse"
#include "datastructures.h"
#include "simparam.h"

const double ENERGY_FACTOR = 3.0 / 2.0;

extern DENS nr, nr_m1, nr_p1;
extern TEMP Tr;
extern ENER Er, Er_m1, Er_p1;
extern dDENS dnr;
extern dENER dEr;
extern dDENS dnr_cn;
extern dENER dEr_cn;
extern CRATE colrate;
extern CRATE colrateRF;
extern POWER Power;

extern double Te_RF[NMESHP];
extern double ne_RF[NMESHP];
extern double nHi_RF[NMESHP];
extern double nH2i_RF[NMESHP];
extern double nH3i_RF[NMESHP];
extern double nHeII_RF[NMESHP];
extern double nHeIII_RF[NMESHP];
extern double dne_RF[NMESHP];
extern double dnHi_RF[NMESHP];
extern double dnH2i_RF[NMESHP];
extern double dnH3i_RF[NMESHP];
extern double dnHeII_RF[NMESHP];
extern double dnHeIII_RF[NMESHP];
extern double PRFe_id[NMESHP];
extern double PRFHi_id[NMESHP];
extern double PRFH2i_id[NMESHP];
extern double PRFH3i_id[NMESHP];
extern double PRFHeII_id[NMESHP];
extern double PRFHeIII_id[NMESHP];
extern double PRFe_array_stat[NMESHP];
extern double PRFHi_array_stat[NMESHP];
extern double PRFH2i_array_stat[NMESHP];
extern double PRFH3i_array_stat[NMESHP];
extern double PRFHeII_array_stat[NMESHP];
extern double PRFHeIII_array_stat[NMESHP];

extern double tmain;
extern double tRF;
extern double dtcalc;

extern double dt;
extern double dt2;
extern double dt3;
extern double dtnew;
extern double drval;
extern double drval2;
extern double drval3;
extern double drval4;
extern double drmax;
extern double drmax2;
extern double drmax3;
extern double drmax4;

extern int cnt_save;
extern int cnt_loop;
extern double dtout;
extern int Nit;
extern int NitOut;

extern double nefact;

extern double pecabs;
extern double nHeI0; // constant op de uiterste gridpunten
extern double nH20;  // constant op de uiterste gridpunten
extern double nCI0;

extern double Pabs1[NMESHP];
extern double Pabs2[NMESHP];
extern double Pabs3[NMESHP];
extern double Pabs4[NMESHP];

extern double Dionh[NMESHP];
extern double Dion[NMESHP];
extern double Vionh[NMESHP];
extern double Vion[NMESHP];
extern double mfpi;
extern double gri;

extern double nui;

extern double coef1, coef2, coef3;

extern Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver;
extern Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver2;
extern Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver3;
extern Eigen::BiCGSTAB<Eigen::SparseMatrix<double>, Eigen::IncompleteLUT<double>> solver4;

extern int mytid;

extern double PIerrorD;
extern double Dfsave;
extern double PIerrorV;
extern double Vfsave;
extern double PIerrorP;

extern double alphaval;

extern double PRFe_array[NMESHP];
extern double PRFHi_array[NMESHP];
extern double PRFH2i_array[NMESHP];
extern double PRFH3i_array[NMESHP];
extern double PRFHeII_array[NMESHP];
extern double PRFHeIII_array[NMESHP];

extern double AverageTimer;
extern double AvTimeStep;
extern double StepTimer[9];
extern double ColTimer;
extern double ParTimer;
extern double RFTimer;
extern double TimeStepTimer;

extern double oldtstep;
extern double olddaccur;

extern DENS n_save;
extern ENER E_save;
extern int convsave;
extern int convsavemin;
extern int convsavemax;
extern int convsavecount;
extern double tdifmax;
extern double tdifne;
extern double tdifTe;
extern double minstopcrit;
extern double tsum;

extern DENS nRF_save;
extern ENER ERF_save;
extern double Rfdif;
extern double Rfdifmax;
extern double dtRFtemp;
extern double dtRFnunc;

extern DENS nr1;
extern DENS nr1bis;
extern ENER Er1;
extern ENER Er1bis;

#endif // GLOBALVARIABLES_H

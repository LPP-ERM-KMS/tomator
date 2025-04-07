#ifndef SOLVER_PRECOMP_H
#define SOLVER_PRECOMP_H

#include "datastructures.h"
#include "simparam.h"
#include "../Eigen/Sparse"

extern DENS nsend;
extern ENER Esend;
extern double vth; //, vcs, mfp, gr, deltat;
extern double DionHFS, DionLFS, VionHFS, VionLFS;
extern double lambda;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> Ds;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> Vs;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> A;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> B;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> C;
extern Eigen::SparseMatrix<double, Eigen::RowMajor> D;
extern double nuion[NMESHP];
extern double mfp;

extern double Z;  // charge state
extern double mu; // mass unit

extern Eigen::VectorXd x;
extern Eigen::VectorXd x2;
extern Eigen::VectorXd x3;
extern Eigen::VectorXd x4;

extern Eigen::VectorXd bnion;
extern Eigen::VectorXd bEion;
extern Eigen::VectorXd bEelec;

extern Eigen::VectorXd bnion1;
extern Eigen::VectorXd bEion1;
extern Eigen::VectorXd bEelec1;

extern Eigen::VectorXd Snneut;
extern Eigen::VectorXd SEneut;
extern Eigen::VectorXd Snion;
extern Eigen::VectorXd SEion;
extern Eigen::VectorXd SEelec;

extern double TH_array[NMESHP];
extern double TH2_array[NMESHP];

extern double nuH_array[NMESHP];
extern double nuH2_array[NMESHP];
extern double nuHi_array[NMESHP];
extern double nuH2i_array[NMESHP];
extern double nuH3i_array[NMESHP];
extern Eigen::VectorXd bnH, bnH2, bnHi, bnH2i, bnH3i;
extern Eigen::VectorXd bEH, bEH2;
extern Eigen::VectorXd bn1H, bn1H2;
extern Eigen::VectorXd bE1H, bE1H2;
extern Eigen::VectorXd SnH, SnH2;
extern Eigen::VectorXd SEH, SEH2;
extern Eigen::VectorXd SsnH2;
extern Eigen::VectorXd SsEH2;

extern double THeI_array[NMESHP];

extern double nuHeI_array[NMESHP];
extern double nuHeII_array[NMESHP];
extern double nuHeIII_array[NMESHP];

extern Eigen::VectorXd bnHeI;
extern Eigen::VectorXd bnHeII;
extern Eigen::VectorXd bnHeIII;
extern Eigen::VectorXd bEHeI;

extern Eigen::VectorXd bn1HeI;
extern Eigen::VectorXd bE1HeI;

extern Eigen::VectorXd SnHeI;
extern Eigen::VectorXd SEHeI;

extern Eigen::VectorXd SsnHeI;
extern Eigen::VectorXd SsEHeI;

extern double nuCI_array[NMESHP];

extern Eigen::VectorXd bnCI;

extern Eigen::VectorXd bn1CI;

extern Eigen::VectorXd bnCII;
extern Eigen::VectorXd bnCIII;
extern Eigen::VectorXd bnCIV;
extern Eigen::VectorXd bnCV;

extern Eigen::VectorXd SnCI;

extern Eigen::VectorXd SsnCI; //

extern double DH[NMESHP];
extern double DH2[NMESHP];

extern Eigen::SparseMatrix<double, Eigen::RowMajor> DHs;

extern Eigen::SparseMatrix<double, Eigen::RowMajor> DH2s;

#endif // SOLVER_PRECOMP_H

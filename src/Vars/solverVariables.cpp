#include "solverVariables.h"

DENS nsend;
ENER Esend;
double vth; //, vcs, mfp, gr, deltat;
double DionHFS, DionLFS, VionHFS, VionLFS;
double lambda;
Eigen::SparseMatrix<double, Eigen::RowMajor> Ds(2 * NMESHP, 2 * NMESHP);
Eigen::SparseMatrix<double, Eigen::RowMajor> Vs(2 * NMESHP, 2 * NMESHP);
Eigen::SparseMatrix<double, Eigen::RowMajor> A;
Eigen::SparseMatrix<double, Eigen::RowMajor> B;
Eigen::SparseMatrix<double, Eigen::RowMajor> C;
Eigen::SparseMatrix<double, Eigen::RowMajor> D;
double nuion[NMESHP];
double mfp;

double Z;  // charge state
double mu; // mass unit

Eigen::VectorXd x(2 * NMESHP);
Eigen::VectorXd x2(2 * NMESHP);
Eigen::VectorXd x3(2 * NMESHP);
Eigen::VectorXd x4(2 * NMESHP);

Eigen::VectorXd bnion(2*NMESHP);
Eigen::VectorXd bEion(2*NMESHP);
Eigen::VectorXd bEelec(2*NMESHP);

Eigen::VectorXd bnion1(2*NMESHP);
Eigen::VectorXd bEion1(2*NMESHP);
Eigen::VectorXd bEelec1(2*NMESHP);

Eigen::VectorXd Snneut(2*NMESHP);
Eigen::VectorXd SEneut(2*NMESHP);
Eigen::VectorXd Snion(2*NMESHP);
Eigen::VectorXd SEion(2*NMESHP);
Eigen::VectorXd SEelec(2*NMESHP);

double TH_array[NMESHP];
double TH2_array[NMESHP];

double nuH_array[NMESHP];
double nuH2_array[NMESHP];
double nuHi_array[NMESHP];
double nuH2i_array[NMESHP];
double nuH3i_array[NMESHP];
Eigen::VectorXd bnH(2*NMESHP), bnH2(2*NMESHP), bnHi(2*NMESHP), bnH2i(2*NMESHP), bnH3i(2*NMESHP);
Eigen::VectorXd bEH(2*NMESHP), bEH2(2*NMESHP);
Eigen::VectorXd bn1H(2*NMESHP), bn1H2(2*NMESHP);
Eigen::VectorXd bE1H(2*NMESHP), bE1H2(2*NMESHP);
Eigen::VectorXd SnH(2*NMESHP), SnH2(2*NMESHP);
Eigen::VectorXd SEH(2*NMESHP), SEH2(2*NMESHP);
Eigen::VectorXd SsnH2(2*NMESHP);
Eigen::VectorXd SsEH2(2*NMESHP);

double THeI_array[NMESHP];

double nuHeI_array[NMESHP];
double nuHeII_array[NMESHP];
double nuHeIII_array[NMESHP];

Eigen::VectorXd bnHeI(2*NMESHP);
Eigen::VectorXd bnHeII(2*NMESHP);
Eigen::VectorXd bnHeIII(2*NMESHP);
Eigen::VectorXd bEHeI(2*NMESHP);

Eigen::VectorXd bn1HeI(2*NMESHP);
Eigen::VectorXd bE1HeI(2*NMESHP);

Eigen::VectorXd SnHeI(2*NMESHP);
Eigen::VectorXd SEHeI(2*NMESHP);

Eigen::VectorXd SsnHeI(2*NMESHP);
Eigen::VectorXd SsEHeI(2*NMESHP);

double nuCI_array[NMESHP];

Eigen::VectorXd bnCI(2*NMESHP);

Eigen::VectorXd bn1CI(2*NMESHP);

Eigen::VectorXd bnCII(2*NMESHP);
Eigen::VectorXd bnCIII(2*NMESHP);
Eigen::VectorXd bnCIV(2*NMESHP);
Eigen::VectorXd bnCV(2*NMESHP);
Eigen::VectorXd SnCI(2*NMESHP);

Eigen::VectorXd SsnCI(2*NMESHP); //

double DH[NMESHP];
double DH2[NMESHP];

Eigen::SparseMatrix<double, Eigen::RowMajor> DHs(2*NMESHP,2*NMESHP);

Eigen::SparseMatrix<double, Eigen::RowMajor> DH2s(2*NMESHP,2*NMESHP);

#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#include "simparam.h"

// Data structures
struct DENS {
    double ne[NMESHP], nH[NMESHP], nH2[NMESHP], nHi[NMESHP], nH2i[NMESHP], nH3i[NMESHP];
    double nHeI[NMESHP], nHeII[NMESHP], nHeIII[NMESHP];
    double nCI[NMESHP], nCII[NMESHP], nCIII[NMESHP], nCIV[NMESHP], nCV[NMESHP];
    double xne[NMESHP], xnH[NMESHP], xnH2[NMESHP], xnHi[NMESHP], xnH2i[NMESHP], xnH3i[NMESHP];
    double xnHeI[NMESHP], xnHeII[NMESHP], xnHeIII[NMESHP];
    double xnCI[NMESHP], xnCII[NMESHP], xnCIII[NMESHP], xnCIV[NMESHP], xnCV[NMESHP];
};
struct ENER {
    double Ee[NMESHP], EH[NMESHP], EH2[NMESHP], EHi[NMESHP], EH2i[NMESHP], EH3i[NMESHP];
    double EHeI[NMESHP], EHeII[NMESHP], EHeIII[NMESHP];
    double xEe[NMESHP], xEH[NMESHP], xEH2[NMESHP], xEHi[NMESHP], xEH2i[NMESHP], xEH3i[NMESHP];
    double xEHeI[NMESHP], xEHeII[NMESHP], xEHeIII[NMESHP];
};
struct TEMP {
    double Te[NMESHP], TH[NMESHP], TH2[NMESHP], THi[NMESHP], TH2i[NMESHP], TH3i[NMESHP];
    double THeI[NMESHP], THeII[NMESHP], THeIII[NMESHP];
};
struct dDENS {
    double dne[NMESHP], dnH[NMESHP], dnH2[NMESHP], dnHi[NMESHP], dnH2i[NMESHP], dnH3i[NMESHP];
    double dnHeI[NMESHP], dnHeII[NMESHP], dnHeIII[NMESHP];
    double dnCI[NMESHP], dnCII[NMESHP], dnCIII[NMESHP], dnCIV[NMESHP], dnCV[NMESHP];
    double dxne[NMESHP], dxnH[NMESHP], dxnH2[NMESHP], dxnHi[NMESHP], dxnH2i[NMESHP], dxnH3i[NMESHP];
    double dxnHeI[NMESHP], dxnHeII[NMESHP], dxnHeIII[NMESHP];
    double dxnCI[NMESHP], dxnCII[NMESHP], dxnCIII[NMESHP], dxnCIV[NMESHP], dxnCV[NMESHP];
};
struct dENER {
    double dEe[NMESHP], dEH[NMESHP], dEH2[NMESHP], dEHi[NMESHP], dEH2i[NMESHP], dEH3i[NMESHP];
    double dEHeI[NMESHP], dEHeII[NMESHP], dEHeIII[NMESHP];
    double dxEe[NMESHP], dxEH[NMESHP], dxEH2[NMESHP], dxEHi[NMESHP], dxEH2i[NMESHP], dxEH3i[NMESHP];
    double dxEHeI[NMESHP], dxEHeII[NMESHP], dxEHeIII[NMESHP];
};
struct CRATE {
    double nue[NMESHP], nuH[NMESHP], nuH2[NMESHP], nuHi[NMESHP], nuH2i[NMESHP], nuH3i[NMESHP];
    double nuHeI[NMESHP], nuHeII[NMESHP], nuHeIII[NMESHP];
};
struct TAUP {
    double taupH[NMESHP], taupH2[NMESHP], taupHi[NMESHP], taupH2i[NMESHP], taupH3i[NMESHP];
    double taupHeI[NMESHP], taupHeII[NMESHP], taupHeIII[NMESHP];
};
struct POWER {
    double Pabs1[NMESHP];
    double Pabs2[NMESHP];
    double Pabs3[NMESHP];
    double Pabs4[NMESHP];
    double PTEST[NMESHP];
};

#endif // DATASTRUCTURES_H

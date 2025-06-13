#include "solver.h"

void solverInit() {
    solver.setTolerance(solvertolerance);
    solver2.setTolerance(solvertolerance);
    solver3.setTolerance(solvertolerance);
    solver4.setTolerance(solvertolerance);

    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                Ds.insert(2 * i, 2 * j) = 0.0;
                Ds.insert(2 * i + 1, 2 * j) = 0.0;
                Ds.insert(2 * i, 2 * j + 1) = 0.0;
                Ds.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                Vs.insert(2 * i, 2 * j) = 0.0;
                Vs.insert(2 * i + 1, 2 * j) = 0.0;
                Vs.insert(2 * i, 2 * j + 1) = 0.0;
                Vs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                DHs.insert(2 * i, 2 * j) = 0.0;
                DHs.insert(2 * i + 1, 2 * j) = 0.0;
                DHs.insert(2 * i, 2 * j + 1) = 0.0;
                DHs.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i + 1) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            } else if (j == i - 1) {
                DH2s.insert(2 * i, 2 * j) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j) = 0.0;
                DH2s.insert(2 * i, 2 * j + 1) = 0.0;
                DH2s.insert(2 * i + 1, 2 * j + 1) = 0.0;
            }
        }
    }
}

void solverAb_x() {
    //////////////////////////////////////////////////////
    //// Transport
    //////////////////////////////////////////////////////

    double tstep = dtnew;

    // Initialization moved to the precompiler
    // double tsolvestart;
    ////////////////////////////////////////////////////
    //// CALC transport of e, process 1
    ////////////////////////////////////////////////////

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        nsend.ne[im] = nsend.nH[im] = nsend.nH2[im] = nsend.nHi[im] = nsend.nH2i[im] = nsend.nH3i[im] = 0.0;
        nsend.nHeI[im] = nsend.nHeII[im] = nsend.nHeIII[im] = 0.0;
        nsend.nCI[im] = nsend.nCII[im] = nsend.nCIII[im] = nsend.nCIV[im] = nsend.nCV[im] = 0.0;
        nsend.xne[im] = nsend.xnH[im] = nsend.xnH2[im] = nsend.xnHi[im] = nsend.xnH2i[im] = nsend.xnH3i[im] = 0.0;
        nsend.xnHeI[im] = nsend.xnHeII[im] = nsend.xnHeIII[im] = 0.0;
        nsend.xnCI[im] = nsend.xnCII[im] = nsend.xnCIII[im] = nsend.xnCIV[im] = nsend.xnCV[im] = 0.0;

        Esend.Ee[im] = Esend.EH[im] = Esend.EH2[im] = Esend.EHi[im] = Esend.EH2i[im] = Esend.EH3i[im] = 0.0;
        Esend.EHeI[im] = Esend.EHeII[im] = Esend.EHeIII[im] = 0.0;
        Esend.xEe[im] = Esend.xEH[im] = Esend.xEH2[im] = Esend.xEHi[im] = Esend.xEH2i[im] = Esend.xEH3i[im] = 0.0;
        Esend.xEHeI[im] = Esend.xEHeII[im] = Esend.xEHeIII[im] = 0.0;
        // cout << nsend.ne[im] << "   " << endl;
    }
    DionLFS = 0.0; // reset
    DionHFS = 0.0;
    int nil = 0;
    int nih = 0;
    VionLFS = 0.0;
    VionHFS = 0.0;
    int im = 1;
    while (aR[im] < lHFS) {
        DionHFS += Dion[im];
        VionHFS += Vionh[im];
        ++nih;
        ++im;
    }

    im = NMESHP - 2;
    while (aR[im] > lLFS) {
        DionLFS += Dion[im];
        VionLFS += Vionh[im];
        ++nil;
        --im;
    }

    im = 0;
    while (aR[im] < lHFS) {
        Dion[im] = DionHFS / nih;
        Vionh[im] = VionHFS / nih;
        ++im;
    }

    im = NMESHP - 1;
    while (aR[im] > lLFS) {
        Dion[im] = DionLFS / nil;
        Vionh[im] = VionLFS / nil;
        --im;
    }

    DionHFS = Dion[0];           // min 1.0e4;
    DionLFS = Dion[NMESHP - 1];  // min 1.0e4;
    VionHFS = Vionh[0];          // min 1.0e4;
    VionLFS = Vionh[NMESHP - 1]; // min 1.0e4;

    // double Tion[NMESHP];
    // Initialization moved to precompiler

    // Following section moved from transportions to the following column (because it does not need te be recalculated for every situation)

    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DLfbba[i] * Dion[i - 1] + DLfbbb[i] * Dion[i] + DRfbbc[i] * Dion[i + 1] + DRfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * Dion[i - 1] + dDLfbbb[i] * Dion[i] + dDRfbbc[i] * Dion[i + 1] + dDRfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * Dion[i - 1] + DLdfbbb[i] * Dion[i] + DRdfbbc[i] * Dion[i + 1] + DRdfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * Dion[i - 1] + dDLdfbbb[i] * Dion[i] + dDRdfbbc[i] * Dion[i + 1] + dDRdfbbb[i] * Dion[i];
            } else if (j == i - 1) {
                Ds.coeffRef(2 * i, 2 * j) = DLfaba[i] * Dion[i - 1] + DLfabb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * Dion[i - 1] + dDLfabb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * Dion[i - 1] + DLdfabb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * Dion[i - 1] + dDLdfabb[i] * Dion[i];
            } else if (j == i + 1) {
                Ds.coeffRef(2 * i, 2 * j) = DRfcbc[i] * Dion[i + 1] + DRfcbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * Dion[i + 1] + dDRfcbb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * Dion[i + 1] + DRdfcbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * Dion[i + 1] + dDRdfcbb[i] * Dion[i];
            }
        }
    }
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DRfbbc[i] * Dion[i + 1] + DRfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfbbc[i] * Dion[i + 1] + dDRfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfbbc[i] * Dion[i + 1] + DRdfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfbbc[i] * Dion[i + 1] + dDRdfbbb[i] * Dion[i];
            } else if (j == i + 1) {
                Ds.coeffRef(2 * i, 2 * j) = DRfcbc[i] * Dion[i + 1] + DRfcbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * Dion[i + 1] + dDRfcbb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * Dion[i + 1] + DRdfcbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * Dion[i + 1] + dDRdfcbb[i] * Dion[i];
            }
        }
    }
    #pragma omp parallel for collapse(2)
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DLfbba[i] * Dion[i - 1] + DLfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * Dion[i - 1] + dDLfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * Dion[i - 1] + DLdfbbb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * Dion[i - 1] + dDLdfbbb[i] * Dion[i];
            } else if (j == i - 1) {
                Ds.coeffRef(2 * i, 2 * j) = DLfaba[i] * Dion[i - 1] + DLfabb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * Dion[i - 1] + dDLfabb[i] * Dion[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * Dion[i - 1] + DLdfabb[i] * Dion[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * Dion[i - 1] + dDLdfabb[i] * Dion[i];
            }
        }
    }

    ////////////////////////////////////////////////////
    //// CALC transport of Hi, process 2
    ////////////////////////////////////////////////////

    Z = 1.0;  // charge state
    mu = 1.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  Tr.THi[im];
        nuion[im] = colrate.nuHi[im];
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }
        bnion(2 * im) = nr.nHi[im];
        bnion(2 * im + 1) = nr.xnHi[im];
        bEion(2 * im) = Er.EHi[im];
        bEion(2 * im + 1) = Er.xEHi[im];
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nHi[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nHi[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nHi[im];
        bnion1(2 * im + 1) = nr1.xnHi[im];
        bEion1(2 * im) = Er1.EHi[im];
        bEion1(2 * im + 1) = Er1.xEHi[im];
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nHi[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nHi[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        SEneut(2 * im) = 0.0;
        SEneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnHi[im];
        Snion(2 * im + 1) = 0.0;
        SEion(2 * im) = dEr.dEHi[im];
        SEion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nHi[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // cout << "1.... Er1.Ee[cc] = "<< Er1.xEe[im] << endl;
    //  for (int id=1; id<NMESHP-1; ++id)
    //  {
    //      Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //      SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //      Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //      SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //      SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //  }
    //  Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    //  SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    //  Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    //  SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    //  SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    //  Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    //  SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    //  Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    //  SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    //  SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    // #include "transportions_v3.03.h"
    

        
    transportions(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}

        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nHi[im] += bnion(2 * im);
        Esend.EHi[im] += bEion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnHi[im] += bnion(2 * im + 1);
        Esend.xEHi[im] += bEion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of H2i, process 4
    ////////////////////////////////////////////////////

    Z = 1.0;  // charge state
    mu = 2.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  Tr.TH2i[im];
        nuion[im] = colrate.nuH2i[im];
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nH2i[im];
        bnion(2 * im + 1) = nr.xnH2i[im];
        bEion(2 * im) = Er.EH2i[im];
        bEion(2 * im + 1) = Er.xEH2i[im];
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nH2i[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nH2i[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nH2i[im];
        bnion1(2 * im + 1) = nr1.xnH2i[im];
        bEion1(2 * im) = Er1.EH2i[im];
        bEion1(2 * im + 1) = Er1.xEH2i[im];
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nH2i[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nH2i[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        SEneut(2 * im) = 0.0;
        SEneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnH2i[im];
        Snion(2 * im + 1) = 0.0;
        SEion(2 * im) = dEr.dEH2i[im];
        SEion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nH2i[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    // #include "transportions_v3.03.h"
    transportions(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nH2i[im] += bnion(2 * im);
        Esend.EH2i[im] += bEion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnH2i[im] += bnion(2 * im + 1);
        Esend.xEH2i[im] += bEion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of H3i, process 5
    ////////////////////////////////////////////////////

    Z = 1.0;  // charge state
    mu = 3.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  Tr.TH3i[im];
        nuion[im] = colrate.nuH3i[im];
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nH3i[im];
        bnion(2 * im + 1) = nr.xnH3i[im];
        bEion(2 * im) = Er.EH3i[im];
        bEion(2 * im + 1) = Er.xEH3i[im];
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nH3i[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nH3i[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nH3i[im];
        bnion1(2 * im + 1) = nr1.xnH3i[im];
        bEion1(2 * im) = Er1.EH3i[im];
        bEion1(2 * im + 1) = Er1.xEH3i[im];
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nH3i[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nH3i[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        SEneut(2 * im) = 0.0;
        SEneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnH3i[im];
        Snion(2 * im + 1) = 0.0;
        SEion(2 * im) = dEr.dEH3i[im];
        SEion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nH3i[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    transportions(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nH3i[im] += bnion(2 * im);
        Esend.EH3i[im] += bEion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnH3i[im] += bnion(2 * im + 1);
        Esend.xEH3i[im] += bEion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of HeII, process 7
    ////////////////////////////////////////////////////

    Z = 1.0;  // charge state
    mu = 4.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  Tr.THeII[im];
        nuion[im] = colrate.nuHeII[im];
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nHeII[im];
        bnion(2 * im + 1) = nr.xnHeII[im];
        bEion(2 * im) = Er.EHeII[im];
        bEion(2 * im + 1) = Er.xEHeII[im];
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nHeII[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nHeII[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nHeII[im];
        bnion1(2 * im + 1) = nr1.xnHeII[im];
        bEion1(2 * im) = Er1.EHeII[im];
        bEion1(2 * im + 1) = Er1.xEHeII[im];
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nHeII[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nHeII[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        SEneut(2 * im) = 0.0;
        SEneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnHeII[im];
        Snion(2 * im + 1) = 0.0;
        SEion(2 * im) = dEr.dEHeII[im];
        SEion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nHeII[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    transportions(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nHeII[im] += bnion(2 * im);
        Esend.EHeII[im] += bEion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnHeII[im] += bnion(2 * im + 1);
        Esend.xEHeII[im] += bEion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of HeIII, process 8
    ////////////////////////////////////////////////////

    Z = 2.0;  // charge state
    mu = 4.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  Tr.THeIII[im];
        nuion[im] = colrate.nuHeIII[im];
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nHeIII[im];
        bnion(2 * im + 1) = nr.xnHeIII[im];
        bEion(2 * im) = Er.EHeIII[im];
        bEion(2 * im + 1) = Er.xEHeIII[im];
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nHeIII[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nHeIII[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nHeIII[im];
        bnion1(2 * im + 1) = nr1.xnHeIII[im];
        bEion1(2 * im) = Er1.EHeIII[im];
        bEion1(2 * im + 1) = Er1.xEHeIII[im];
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nHeIII[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nHeIII[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        SEneut(2 * im) = 0.0;
        SEneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnHeIII[im];
        Snion(2 * im + 1) = 0.0;
        SEion(2 * im) = dEr.dEHeIII[im];
        SEion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nHeIII[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    transportions(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nHeIII[im] += bnion(2 * im);
        Esend.EHeIII[im] += bEion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnHeIII[im] += bnion(2 * im + 1);
        Esend.xEHeIII[im] += bEion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of CII, process 10
    ////////////////////////////////////////////////////

    Z = 1.0;   // charge state
    mu = 12.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  3.0;
        nuion[im] = 0.0;
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nCII[im];
        bnion(2 * im + 1) = nr.xnCII[im];
        bEion(2 * im) = 3.0 / 2.0 * nr.nCII[im] * 3.0;
        bEion(2 * im + 1) = 3.0 / 2.0 * nr.xnCII[im] * 3.0;
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nCII[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nCII[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nCII[im];
        bnion1(2 * im + 1) = nr1.xnCII[im];
        bEion1(2 * im) = 3.0 / 2.0 * nr1.nCII[im] * 3.0;
        bEion1(2 * im + 1) = 3.0 / 2.0 * nr1.xnCII[im] * 3.0;
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nCII[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nCII[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnCII[im];
        Snion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nCII[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    if (bimpur) {
        transportions(tstep);
    }

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nCII[im] += bnion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnCII[im] += bnion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of CIII, process 11
    ////////////////////////////////////////////////////

    Z = 2.0;   // charge state
    mu = 12.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  3.0;
        nuion[im] = 0.0;
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nCIII[im];
        bnion(2 * im + 1) = nr.xnCIII[im];
        bEion(2 * im) = 3.0 / 2.0 * nr.nCIII[im] * 3.0;
        bEion(2 * im + 1) = 3.0 / 2.0 * nr.xnCIII[im] * 3.0;
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nCIII[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nCIII[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nCIII[im];
        bnion1(2 * im + 1) = nr1.xnCIII[im];
        bEion1(2 * im) = 3.0 / 2.0 * nr1.nCIII[im] * 3.0;
        bEion1(2 * im + 1) = 3.0 / 2.0 * nr1.xnCIII[im] * 3.0;
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nCIII[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nCIII[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnCIII[im];
        Snion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nCIII[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    if (bimpur) {
        transportions(tstep);
    }

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        if (nsend.nCIII[im] != 0.0) {
            cout << "(nsend.nCIII[im] ~= 0.0) before" << endl;
        }
        nsend.nCIII[im] += bnion(2 * im);
        if (nsend.nCIII[im] < 0.0) {
            cout << "(nsend.nCIII[im] < 0.0) after" << endl;
        }
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnCIII[im] += bnion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of CIV, process 12
    ////////////////////////////////////////////////////

    Z = 3.0;   // charge state
    mu = 12.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  3.0;
        nuion[im] = 0.0;
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nCIV[im];
        bnion(2 * im + 1) = nr.xnCIV[im];
        bEion(2 * im) = 3.0 / 2.0 * nr.nCIV[im] * 3.0;
        bEion(2 * im + 1) = 3.0 / 2.0 * nr.xnCIV[im] * 3.0;
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nCIV[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nCIV[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nCIV[im];
        bnion1(2 * im + 1) = nr1.xnCIV[im];
        bEion1(2 * im) = 3.0 / 2.0 * nr1.nCIV[im] * 3.0;
        bEion1(2 * im + 1) = 3.0 / 2.0 * nr1.xnCIV[im] * 3.0;
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nCIV[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nCIV[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnCIV[im];
        Snion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nCIV[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    if (bimpur) {
        transportions(tstep);
    }

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nCIV[im] += bnion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnCIV[im] += bnion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of CV, process 13
    ////////////////////////////////////////////////////

    Z = 4.0;   // charge state
    mu = 12.0; // mass unit

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Tion[im]        =  3.0;
        nuion[im] = 0.0;
        if (nuion[im] < 1.0e4) {
            nuion[im] = 1.0e4;
        }

        bnion(2 * im) = nr.nCV[im];
        bnion(2 * im + 1) = nr.xnCV[im];
        bEion(2 * im) = 3.0 / 2.0 * nr.nCV[im] * 3.0;
        bEion(2 * im + 1) = 3.0 / 2.0 * nr.xnCV[im] * 3.0;
        bEelec(2 * im) = Er.Ee[im] * Z * nr.nCV[im] / nr.ne[im];
        bEelec(2 * im + 1) = Er.xEe[im] * Z * nr.nCV[im] / nr.ne[im];

        bnion1(2 * im) = nr1.nCV[im];
        bnion1(2 * im + 1) = nr1.xnCV[im];
        bEion1(2 * im) = 3.0 / 2.0 * nr1.nCV[im] * 3.0;
        bEion1(2 * im + 1) = 3.0 / 2.0 * nr1.xnCV[im] * 3.0;
        bEelec1(2 * im) = Er1.Ee[im] * Z * nr1.nCV[im] / nr1.ne[im];
        bEelec1(2 * im + 1) = Er1.xEe[im] * Z * nr1.nCV[im] / nr1.ne[im];

        Snneut(2 * im) = 0.0;
        Snneut(2 * im + 1) = 0.0;
        Snion(2 * im) = dnr.dnCV[im];
        Snion(2 * im + 1) = 0.0;
        SEelec(2 * im) = dEr.dEe[im] * Z * nr.nCV[im] / nr.ne[im];
        SEelec(2 * im + 1) = 0.0;

        // Ssnion(2*im)    =  0.0;                 SsEion(2*im)     =  0.0;                 SsEelec(2*im)     =  0.0;
        // Ssnion(2*im+1)  =  0.0;                 SsEion(2*im+1)   =  0.0;                 SsEelec(2*im+1)   =  0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }

    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     Snneut(2*id+1) = (Snneut(2*(id+1)) - Snneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEneut(2*id+1) = (SEneut(2*(id+1)) - SEneut(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     Snion(2*id+1)  = (Snion(2*(id+1))  - Snion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEion(2*id+1)  = (SEion(2*(id+1))  - SEion(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEelec(2*id+1) = (SEelec(2*(id+1)) - SEelec(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    // }
    // Snneut(1) = (Snneut(2*(1)) - Snneut(2*(0))) / (aR[1] - aR[0]) ;
    // SEneut(1) = (SEneut(2*(1)) - SEneut(2*(0))) / (aR[1] - aR[0]) ;
    // Snion(1)  = (Snion(2*(1))  - Snion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEion(1)  = (SEion(2*(1))  - SEion(2*(0)))  / (aR[1] - aR[0]) ;
    // SEelec(1) = (SEelec(2*(1)) - SEelec(2*(0))) / (aR[1] - aR[0]) ;
    // Snneut(2*(NMESHP-1)+1) = (Snneut(2*(NMESHP-1)) - Snneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEneut(2*(NMESHP-1)+1) = (SEneut(2*(NMESHP-1)) - SEneut(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // Snion(2*(NMESHP-1)+1)  = (Snion(2*(NMESHP-1))  - Snion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEion(2*(NMESHP-1)+1)  = (SEion(2*(NMESHP-1))  - SEion(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEelec(2*(NMESHP-1)+1) = (SEelec(2*(NMESHP-1)) - SEelec(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    if (bimpur) {
        transportions(tstep);
    }

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // if (im==0) {cout << nsend.ne[im] << endl;}
        nsend.ne[im] += bnion(2 * im) * Z;
        Esend.Ee[im] += bEelec(2 * im);
        nsend.nCV[im] += bnion(2 * im);
        nsend.xne[im] += bnion(2 * im + 1) * Z;
        Esend.xEe[im] += bEelec(2 * im + 1);
        nsend.xnCV[im] += bnion(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of H + H2
    ////////////////////////////////////////////////////

    // Moved to precompiler

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // nH_neut[im]   =  0.0;                  EH_neut[im]    =  0.0;
        // nH2_neut[im]  =  0.0;                  EH2_neut[im]   =  0.0;
        TH_array[im] = Tr.TH[im];
        TH2_array[im] = Tr.TH2[im];
        // THi_array[im]  =  Tr.THi[im] ;
        // TH2i_array[im]  =  Tr.TH2i[im] ;
        // TH3i_array[im]  =  Tr.TH3i[im] ;
        nuH_array[im] = colrate.nuH[im];
        if (nuH_array[im] < 1.0e3) {
            nuH_array[im] = 1.0e3;
        }
        nuH2_array[im] = colrate.nuH2[im];
        if (nuH2_array[im] < 1.0e3) {
            nuH2_array[im] = 1.0e3;
        }
        nuHi_array[im] = colrate.nuHi[im];
        if (nuHi_array[im] < 1.0e3) {
            nuHi_array[im] = 1.0e3;
        }
        nuH2i_array[im] = colrate.nuH2i[im];
        if (nuH2i_array[im] < 1.0e3) {
            nuH2i_array[im] = 1.0e3;
        }
        nuH3i_array[im] = colrate.nuH3i[im];
        if (nuH3i_array[im] < 1.0e3) {
            nuH3i_array[im] = 1.0e3;
        }

        bnH(2 * im) = nr.nH[im];
        bEH(2 * im) = Er.EH[im];
        bnH(2 * im + 1) = nr.xnH[im];
        bEH(2 * im + 1) = Er.xEH[im];
        bnH2(2 * im) = nr.nH2[im];
        bEH2(2 * im) = Er.EH2[im];
        bnH2(2 * im + 1) = nr.xnH2[im];
        bEH2(2 * im + 1) = Er.xEH2[im];

        bn1H(2 * im) = nr1.nH[im];
        bE1H(2 * im) = Er1.EH[im];
        bn1H(2 * im + 1) = nr1.xnH[im];
        bE1H(2 * im + 1) = Er1.xEH[im];
        bn1H2(2 * im) = nr1.nH2[im];
        bE1H2(2 * im) = Er1.EH2[im];
        bn1H2(2 * im + 1) = nr1.xnH2[im];
        bE1H2(2 * im + 1) = Er1.xEH2[im];

        bnHi(2 * im) = nr.nHi[im];
        bnHi(2 * im + 1) = nr.xnHi[im];
        bnH2i(2 * im) = nr.nH2i[im];
        bnH2i(2 * im + 1) = nr.xnH2i[im];
        bnH3i(2 * im) = nr.nH3i[im];
        bnH3i(2 * im + 1) = nr.xnH3i[im];

        SnH(2 * im) = dnr.dnH[im];
        SEH(2 * im) = dEr.dEH[im];
        SnH(2 * im + 1) = 0.0;
        SEH(2 * im + 1) = 0.0;
        SnH2(2 * im) = dnr.dnH2[im];
        SEH2(2 * im) = dEr.dEH2[im];
        SnH2(2 * im + 1) = 0.0;
        SEH2(2 * im + 1) = 0.0;

        SsnH2(2 * im) = 0.0;
        SsEH2(2 * im) = 0.0;
        SsnH2(2 * im + 1) = 0.0;
        SsEH2(2 * im + 1) = 0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
        // cout << "    " << SnH(2*im) << endl;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     SnH(2*id+1) = (SnH(2*(id+1)) - SnH(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SEH(2*id+1) = (SEH(2*(id+1)) - SEH(2*(id-1))) / (aR[id+1] - aR[id-1]) ;
    //     SnH2(2*id+1)  = (SnH2(2*(id+1))  - SnH2(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    //     SEH2(2*id+1)  = (SEH2(2*(id+1))  - SEH2(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    // }
    // SnH(1) = (SnH(2*(1)) - SnH(2*(0))) / (aR[1] - aR[0]) ;
    // SEH(1) = (SEH(2*(1)) - SEH(2*(0))) / (aR[1] - aR[0]) ;
    // SnH2(1)  = (SnH2(2*(1))  - SnH2(2*(0)))  / (aR[1] - aR[0]) ;
    // SEH2(1)  = (SEH2(2*(1))  - SEH2(2*(0)))  / (aR[1] - aR[0]) ;
    // SnH(2*(NMESHP-1)+1) = (SnH(2*(NMESHP-1)) - SnH(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEH(2*(NMESHP-1)+1) = (SEH(2*(NMESHP-1)) - SEH(2*(NMESHP-2))) / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SnH2(2*(NMESHP-1)+1)  = (SnH2(2*(NMESHP-1))  - SnH2(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    // SEH2(2*(NMESHP-1)+1)  = (SEH2(2*(NMESHP-1))  - SEH2(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;
    transpH(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        nsend.nH[im] += bnH(2 * im);
        Esend.EH[im] += bEH(2 * im);
        nsend.nH2[im] += bnH2(2 * im);
        Esend.EH2[im] += bEH2(2 * im);
        nsend.xnH[im] += bnH(2 * im + 1);
        Esend.xEH[im] += bEH(2 * im + 1);
        nsend.xnH2[im] += bnH2(2 * im + 1);
        Esend.xEH2[im] += bEH2(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of HeI, process 6
    ////////////////////////////////////////////////////

    // Moved to precompiler

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        bnHeI(2 * im) = nr.nHeI[im];
        bEHeI(2 * im) = Er.EHeI[im];
        bn1HeI(2 * im) = nr1.nHeI[im];
        bE1HeI(2 * im) = Er1.EHeI[im];
        THeI_array[im] = Tr.THeI[im];
        // THeII_array[im]   =  Tr.THeII[im];
        // THeIII_array[im]   =  Tr.THeIII[im];
        nuHeI_array[im] = colrate.nuHeI[im];
        if (nuHeI_array[im] < 1.0e3) {
            nuHeI_array[im] = 1.0e3;
        }
        nuHeII_array[im] = colrate.nuHeII[im];
        if (nuHeII_array[im] < 1.0e3) {
            nuHeII_array[im] = 1.0e3;
        }
        nuHeIII_array[im] = colrate.nuHeIII[im];
        if (nuHeIII_array[im] < 1.0e3) {
            nuHeIII_array[im] = 1.0e3;
        }

        bnHeII(2 * im) = nr.nHeII[im];
        bnHeIII(2 * im) = nr.nHeIII[im];
        bnHeI(2 * im + 1) = nr.xnHeI[im];
        bEHeI(2 * im + 1) = Er.xEHeI[im];

        bn1HeI(2 * im + 1) = nr1.xnHeI[im];
        bE1HeI(2 * im + 1) = Er1.xEHeI[im];

        bnHeII(2 * im + 1) = nr.xnHeII[im];
        bnHeIII(2 * im + 1) = nr.xnHeIII[im];

        SnHeI(2 * im) = dnr.dnHeI[im];
        SEHeI(2 * im) = dEr.dEHeI[im];
        SnHeI(2 * im + 1) = 0.0;
        SEHeI(2 * im + 1) = 0.0;

        SsnHeI(2 * im) = 0.0;
        SsEHeI(2 * im) = 0.0;
        SsnHeI(2 * im + 1) = 0.0;
        SsEHeI(2 * im + 1) = 0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    transpHe(tstep);

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        nsend.nHeI[im] += bnHeI(2 * im);
        Esend.EHeI[im] += bEHeI(2 * im);
        nsend.xnHeI[im] += bnHeI(2 * im + 1);
        Esend.xEHeI[im] += bEHeI(2 * im + 1);
    }

    ////////////////////////////////////////////////////
    //// CALC transport of Carbon NEUTRALS
    ////////////////////////////////////////////////////

    // Moved to precompiler

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        nuCI_array[im] = 0.33 * (colrate.nuHeI[im] + colrate.nuH[im] + colrate.nuH2[im]);
        if (nuCI_array[im] < 1.0e3) {
            nuCI_array[im] = 1.0e3;
        }
        // nuCII_array[im]  =  0.25*(colrate.nuHeII[im]+colrate.nuHi[im]+colrate.nuH2i[im]+colrate.nuH3i[im]);    if (nuCII_array[im]<1.0e3) {nuCII_array[im]=1.0e3;}
        // nuCIII_array[im]  =  colrate.nuHeIII[im];                                                           if (nuCIII_array[im]<1.0e3) {nuCIII_array[im]=1.0e3;}

        bnCI(2 * im) = nr.nCI[im];
        bnCI(2 * im + 1) = nr.xnCI[im];

        bn1CI(2 * im) = nr1.nCI[im];
        bn1CI(2 * im + 1) = nr1.xnCI[im];

        bnCII(2 * im) = nr.nCII[im];
        bnCII(2 * im + 1) = nr.xnCII[im];
        bnCIII(2 * im) = nr.nCIII[im];
        bnCIII(2 * im + 1) = nr.xnCIII[im];
        bnCIV(2 * im) = nr.nCIV[im];
        bnCIV(2 * im + 1) = nr.xnCIV[im];
        bnCV(2 * im) = nr.nCV[im];
        bnCV(2 * im + 1) = nr.xnCV[im];
        SnCI(2 * im) = dnr.dnCI[im];
        SnCI(2 * im + 1) = 0.0;
        SsnCI(2 * im) = 0.0;
        SsnCI(2 * im + 1) = 0.0;

        x(2 * im) = 0.0;
        x(2 * im + 1) = 0.0;
    }
    // for (int id=1; id<NMESHP-1; ++id)
    // {
    //     SnCI(2*id+1)  = (SnCI(2*(id+1))  - SnCI(2*(id-1)))  / (aR[id+1] - aR[id-1]) ;
    // }
    // SnCI(1)  = (SnCI(2*(1))  - SnCI(2*(0)))  / (aR[1] - aR[0]) ;
    // SnCI(2*(NMESHP-1)+1)  = (SnCI(2*(NMESHP-1))  - SnCI(2*(NMESHP-2)))  / (aR[NMESHP-1] - aR[NMESHP-2]) ;

    if (bimpur) {
        // #include "transportC_v2.41.h"
        transpC(tstep);
    }

//#pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        nsend.nCI[im] += bnCI(2 * im);
        nsend.xnCI[im] += bnCI(2 * im + 1);
    }

////////////////////////////////////////////////////
//// END TRANSPORT,
////////////////////////////////////////////////////
//////////////////////////////////////////////////////
    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        dnr_cn.dnH[im] = 0.0;
        dnr_cn.dnH2[im] = 0.0;
        dnr_cn.dnHi[im] = 0.0;
        dnr_cn.dnH2i[im] = 0.0;
        dnr_cn.dnH3i[im] = 0.0;
        dnr_cn.dnHeI[im] = 0.0;
        dnr_cn.dnHeII[im] = 0.0;
        dnr_cn.dnHeIII[im] = 0.0;
        dnr_cn.dnCI[im] = 0.0;
        dnr_cn.dnCII[im] = 0.0;
        dnr_cn.dnCIII[im] = 0.0;
        dnr_cn.dnCIV[im] = 0.0;
        dnr_cn.dnCV[im] = 0.0;
        dEr_cn.dEH[im] = 0.0;
        dEr_cn.dEH2[im] = 0.0;
        dEr_cn.dEHi[im] = 0.0;
        dEr_cn.dEH2i[im] = 0.0;
        dEr_cn.dEH3i[im] = 0.0;
        dEr_cn.dEHeI[im] = 0.0;
        dEr_cn.dEHeII[im] = 0.0;
        dEr_cn.dEHeIII[im] = 0.0;
        dnr_cn.dne[im] = 0.0;
        dEr_cn.dEe[im] = 0.0;

        dnr_cn.dxnH[im] = 0.0;
        dnr_cn.dxnH2[im] = 0.0;
        dnr_cn.dxnHi[im] = 0.0;
        dnr_cn.dxnH2i[im] = 0.0;
        dnr_cn.dxnH3i[im] = 0.0;
        dnr_cn.dxnHeI[im] = 0.0;
        dnr_cn.dxnHeII[im] = 0.0;
        dnr_cn.dxnHeIII[im] = 0.0;
        dnr_cn.dxnCI[im] = 0.0;
        dnr_cn.dxnCII[im] = 0.0;
        dnr_cn.dxnCIII[im] = 0.0;
        dnr_cn.dxnCIV[im] = 0.0;
        dnr_cn.dxnCV[im] = 0.0;
        dEr_cn.dxEH[im] = 0.0;
        dEr_cn.dxEH2[im] = 0.0;
        dEr_cn.dxEHi[im] = 0.0;
        dEr_cn.dxEH2i[im] = 0.0;
        dEr_cn.dxEH3i[im] = 0.0;
        dEr_cn.dxEHeI[im] = 0.0;
        dEr_cn.dxEHeII[im] = 0.0;
        dEr_cn.dxEHeIII[im] = 0.0;
        dnr_cn.dxne[im] = 0.0;
        dEr_cn.dxEe[im] = 0.0;
    }
    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        dnr_cn.dnH[im] = (nsend.nH[im] - nr.nH[im]) / dtnew;             // cout << "im = " << im << ", dnH " << dnr_cn.dnH[im] << endl;
        dnr_cn.dnH2[im] = (nsend.nH2[im] - nr.nH2[im]) / dtnew;          // cout << "im = " << im << ", dnH2 " << dnr_cn.dnH2[im] << endl;
        dnr_cn.dnHi[im] = (nsend.nHi[im] - nr.nHi[im]) / dtnew;          // cout << "im = " << im << ", nHi " << nr.nHi[im] << ", nHi send " << nsend.nHi[im] << ", dnHi " << dnr_cn.dnHi[im] << endl;
        dnr_cn.dnH2i[im] = (nsend.nH2i[im] - nr.nH2i[im]) / dtnew;       // cout << "im = " << im << ", dnH2i " << dnr_cn.dnH2i[im] << endl;
        dnr_cn.dnH3i[im] = (nsend.nH3i[im] - nr.nH3i[im]) / dtnew;       // cout << "im = " << im << ", dnH3i " << dnr_cn.dnH3i[im] << endl;
        dnr_cn.dnHeI[im] = (nsend.nHeI[im] - nr.nHeI[im]) / dtnew;       // cout << "im = " << im << ", dnHeI " << dnr_cn.dnHeI[im] << endl;
        dnr_cn.dnHeII[im] = (nsend.nHeII[im] - nr.nHeII[im]) / dtnew;    // cout << "im = " << im << ", dnHeII " << dnr_cn.dnHeII[im] << endl;
        dnr_cn.dnHeIII[im] = (nsend.nHeIII[im] - nr.nHeIII[im]) / dtnew; // cout << "im = " << im << ", dnHeIII " << dnr_cn.dnHeIII[im] << endl;
        dnr_cn.dnCI[im] = (nsend.nCI[im] - nr.nCI[im]) / dtnew;          // cout << "im = " << im << ", dnCI " << dnr_cn.dnCI[im] << endl;
        dnr_cn.dnCII[im] = (nsend.nCII[im] - nr.nCII[im]) / dtnew;       // cout << "im = " << im << ", dnCII " << dnr_cn.dnCII[im] << endl;
        dnr_cn.dnCIII[im] = (nsend.nCIII[im] - nr.nCIII[im]) / dtnew;    // cout << "im = " << im << ", dnCIII " << dnr_cn.dnCIII[im] << endl;
        dnr_cn.dnCIV[im] = (nsend.nCIV[im] - nr.nCIV[im]) / dtnew;       // cout << "im = " << im << ", dnCIII " << dnr_cn.dnCIII[im] << endl;
        dnr_cn.dnCV[im] = (nsend.nCV[im] - nr.nCV[im]) / dtnew;          // cout << "im = " << im << ", dnCIII " << dnr_cn.dnCIII[im] << endl;
        dEr_cn.dEH[im] = (Esend.EH[im] - Er.EH[im]) / dtnew;             // cout << "im = " << im << ", nEH " << Er.EH[im] << ", nEH send " << Esend.EH[im] << ", dEH " << dEr_cn.dEH[im] << endl;
        dEr_cn.dEH2[im] = (Esend.EH2[im] - Er.EH2[im]) / dtnew;          // cout << "im = " << im << ", dEH2 " << dEr_cn.dEH2[im] << endl;
        dEr_cn.dEHi[im] = (Esend.EHi[im] - Er.EHi[im]) / dtnew;          // cout << "im = " << im << ", dEHi " << dEr_cn.dEHi[im] << endl;
        dEr_cn.dEH2i[im] = (Esend.EH2i[im] - Er.EH2i[im]) / dtnew;       // cout << "im = " << im << ", dEH2i " << dEr_cn.dEH2i[im] << endl;
        dEr_cn.dEH3i[im] = (Esend.EH3i[im] - Er.EH3i[im]) / dtnew;       // cout << "im = " << im << ", dEH3i " << dEr_cn.dEH3i[im] << endl;
        dEr_cn.dEHeI[im] = (Esend.EHeI[im] - Er.EHeI[im]) / dtnew;       // cout << "im = " << im << ", dEHeI " << dEr_cn.dEHeI[im] << endl;
        dEr_cn.dEHeII[im] = (Esend.EHeII[im] - Er.EHeII[im]) / dtnew;    // cout << "im = " << im << ", dEHeII " << dEr_cn.dEHeII[im] << endl;
        dEr_cn.dEHeIII[im] = (Esend.EHeIII[im] - Er.EHeIII[im]) / dtnew; // cout << "im = " << im << ", dEHeIII " << dEr_cn.dEHeIII[im] << endl;


        dnr_cn.dne[im] = (nsend.ne[im] - nr.ne[im]) / dtnew;
        // if (im==0) {cout << "im = " << im << ", ne " << nr.ne[im] << ", ne send " << nsend.ne[im] << ", dne " << dnr_cn.dne[im] << endl;}
        dEr_cn.dEe[im] = (Esend.Ee[im] - Er.Ee[im]) / dtnew; // cout << "im = " << im << ", dEe " << dEr_cn.dEe[im] << endl;

        dnr_cn.dxnH[im] = (nsend.xnH[im] - nr.xnH[im]) / dtnew;             // cout << "im = " << im << ", dnH " << dnr_cn.dnH[im] << endl;
        dnr_cn.dxnH2[im] = (nsend.xnH2[im] - nr.xnH2[im]) / dtnew;          // cout << "im = " << im << ", dnH2 " << dnr_cn.dnH2[im] << endl;
        dnr_cn.dxnHi[im] = (nsend.xnHi[im] - nr.xnHi[im]) / dtnew;          // cout << "im = " << im << ", dnHi " << dnr_cn.dnHi[im] << endl;
        dnr_cn.dxnH2i[im] = (nsend.xnH2i[im] - nr.xnH2i[im]) / dtnew;       // cout << "im = " << im << ", dnH2i " << dnr_cn.dnH2i[im] << endl;
        dnr_cn.dxnH3i[im] = (nsend.xnH3i[im] - nr.xnH3i[im]) / dtnew;       // cout << "im = " << im << ", dnH3i " << dnr_cn.dnH3i[im] << endl;
        dnr_cn.dxnHeI[im] = (nsend.xnHeI[im] - nr.xnHeI[im]) / dtnew;       // cout << "im = " << im << ", dnHeI " << dnr_cn.dnHeI[im] << endl;
        dnr_cn.dxnHeII[im] = (nsend.xnHeII[im] - nr.xnHeII[im]) / dtnew;    // cout << "im = " << im << ", dnHeII " << dnr_cn.dnHeII[im] << endl;
        dnr_cn.dxnHeIII[im] = (nsend.xnHeIII[im] - nr.xnHeIII[im]) / dtnew; // cout << "im = " << im << ", dnHeIII " << dnr_cn.dnHeIII[im] << endl;
        dnr_cn.dxnCI[im] = (nsend.xnCI[im] - nr.xnCI[im]) / dtnew;          // cout << "im = " << im << ", dnCI " << dnr_cn.dnCI[im] << endl;
        dnr_cn.dxnCII[im] = (nsend.xnCII[im] - nr.xnCII[im]) / dtnew;       // cout << "im = " << im << ", dnCII " << dnr_cn.dnCII[im] << endl;
        dnr_cn.dxnCIII[im] = (nsend.xnCIII[im] - nr.xnCIII[im]) / dtnew;    // cout << "im = " << im << ", dnCIII " << dnr_cn.dnCIII[im] << endl;
        dnr_cn.dxnCIV[im] = (nsend.xnCIV[im] - nr.xnCIV[im]) / dtnew;       // cout << "im = " << im << ", dnCIII " << dnr_cn.dnCIII[im] << endl;
        dnr_cn.dxnCV[im] = (nsend.xnCV[im] - nr.xnCV[im]) / dtnew;          // cout << "im = " << im << ", dnCIII " << dnr_cn.dnCIII[im] << endl;
        dEr_cn.dxEH[im] = (Esend.xEH[im] - Er.xEH[im]) / dtnew;             // cout << "im = " << im << ", dEH " << dEr_cn.dEH[im] << endl;
        dEr_cn.dxEH2[im] = (Esend.xEH2[im] - Er.xEH2[im]) / dtnew;          // cout << "im = " << im << ", dEH2 " << dEr_cn.dEH2[im] << endl;
        dEr_cn.dxEHi[im] = (Esend.xEHi[im] - Er.xEHi[im]) / dtnew;          // cout << "im = " << im << ", dEHi " << dEr_cn.dEHi[im] << endl;
        dEr_cn.dxEH2i[im] = (Esend.xEH2i[im] - Er.xEH2i[im]) / dtnew;       // cout << "im = " << im << ", dEH2i " << dEr_cn.dEH2i[im] << endl;
        dEr_cn.dxEH3i[im] = (Esend.xEH3i[im] - Er.xEH3i[im]) / dtnew;       // cout << "im = " << im << ", dEH3i " << dEr_cn.dEH3i[im] << endl;
        dEr_cn.dxEHeI[im] = (Esend.xEHeI[im] - Er.xEHeI[im]) / dtnew;       // cout << "im = " << im << ", dEHeI " << dEr_cn.dEHeI[im] << endl;
        dEr_cn.dxEHeII[im] = (Esend.xEHeII[im] - Er.xEHeII[im]) / dtnew;    // cout << "im = " << im << ", dEHeII " << dEr_cn.dEHeII[im] << endl;
        dEr_cn.dxEHeIII[im] = (Esend.xEHeIII[im] - Er.xEHeIII[im]) / dtnew; // cout << "im = " << im << ", dEHeIII " << dEr_cn.dEHeIII[im] << endl;

        dnr_cn.dxne[im] = (nsend.xne[im] - nr.xne[im]) / dtnew; // cout << "im = " << im << ", dne " << dnr_cn.dne[im] << endl;
        dEr_cn.dxEe[im] = (Esend.xEe[im] - Er.xEe[im]) / dtnew; // cout << "im = " << im << ", dEe " << dEr_cn.dEe[im] << endl;
    }
    //////////////////////////////////////////////////////
    //// update with old time step
    //////////////////////////////////////////////////////
    /*nr.ne     = nr_trHi.ne;         Er.Ee     = Er_trHi.Ee;
    nr.nH     = nr_trHyd.nH;          Er.EH     = Er_trHyd.EH;
    nr.nH2    = nr_trHyd.nH2;         Er.EH2    = Er_trHyd.EH2;
    nr.nHi    = nr_trHi.nHi;          Er.EHi    = Er_trHi.EHi;
    nr.nH2i   = nr_trH2i.nH2i;        Er.EH2i   = Er_trH2i.EH2i;
    nr.nH3i   = nr_trH3i.nH3i;        Er.EH3i   = Er_trH3i.EH3i;
    nr.nHeI   = nr_trHeI.nHeI;        Er.EHeI   = Er_trHeI.EHeI;
    nr.nHeII  = nr_trHeII.nHeII;      Er.EHeII  = Er_trHeII.EHeII;
    nr.nHeIII = nr_trHeIII.nHeIII;    Er.EHeIII = Er_trHeIII.EHeIII;
    nr.nCI    = nr_trCI.nCI;
    nr.nCII   = nr_trCII.nCII;
    nr.nCIII  = nr_trCIII.nCIII;*/
}

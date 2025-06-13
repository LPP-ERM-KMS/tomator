#include "coupledpower.h"

void bgray_func() {
    // double powe=0.0;
    pecabs = 0.0;
    double eprof[NMESHP], inteprof = 0.0;

    int id = 0;
    while (aR[id] < Rdep)
        ++id;
    // cout << aR[id] << endl;
    double necalc = max(nr.ne[id], nr.ne[id + 1]); // cout << id << endl; //nr.ne[90]; if (necalc>1.0e13) {necalc=1.0e13;}
    double Tecalc = max(Tr.Te[id], Tr.Te[id + 1]); // Tr.Te[90];
    if (Tecalc >= 5000.0) {
        Tecalc = 5000.0;
    }
    if (Tecalc >= 50.0) {
        pecabs = interpolate2D<sizeof(Mg1[0]) / sizeof(double), sizeof(Mg1) / sizeof(Mg1[0])>(Mne, MTe, Mg1, necalc * 1.0e6, Tecalc);
    } else {
        double pecabs50 = interpolate2D<sizeof(Mg1[0]) / sizeof(double), sizeof(Mg1) / sizeof(Mg1[0])>(Mne, MTe, Mg1, necalc * 1.0e6, 50.0); // cout << pecabs50 << "   " ;
        double pecabs55 = interpolate2D<sizeof(Mg1[0]) / sizeof(double), sizeof(Mg1) / sizeof(Mg1[0])>(Mne, MTe, Mg1, necalc * 1.0e6, 54.8); // cout << pecabs55 << endl ;
        double dTe = log10(54.8 / 50.0);
        double dpabs = log10(pecabs55 / pecabs50);
        double cst = log10(50.0) * dpabs / dTe - log10(pecabs50);
        pecabs = pow(10.0, log10(Tecalc) * dpabs / dTe - cst); // cout << " ne = " << necalc << "    Te = " << Tecalc << "     pecabs = " << pecabs << endl;
    }
    pecabs = pecabs / (1 - ((1 - pecabs) * (1 - muw)));

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof[id] = aR[id] * exp(-pow((aR[id] - Rdep) / widthech, 2.0));
        eprof[id] += echbackground * aR[id] * exp(-pow((aR[id] - Rdep) / a, 2.0)); // background
    }
    for (int id = 0; id < NMESHP - 1; ++id) {                                                                     // line integrated density
        inteprof += 2.0 * b * pi * 0.5 * (eprof[id] + eprof[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)); // 20181016
    }
    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof[id] = eprof[id] / inteprof;
        // cout << id << "  " << eprof[id]*(pow(aR[NMESHP-1],2.0)-pow(aR[0],2.0))/(pow(aR[id+1],2.0)-pow(aR[id],2.0)) << endl;
    }
#pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        PRFe_array[id] = pecabs * 6.24e18 * eprof[id] * (Prf * 1e3); // 20181016
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
}
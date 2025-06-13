#include "coupledpower.h"

void bram_func() {
    // #include "ram_lookuptable_v3.__.h"
    //
    // double powe=0.0;
    pecabs = 0.15;
    double eprof[NMESHP], inteprof = 0.0;

    int id = 0;
    while (aR[id] < Rdep)
        ++id;
    // cout << aR[id] << endl;
    //  double necalc=max(nr.ne[id],nr.ne[id+1]); //cout << id << endl; //nr.ne[90]; if (necalc>1.0e13) {necalc=1.0e13;}
    //  double Tecalc=max(Tr.Te[id],Tr.Te[id+1]);//Tr.Te[90];
    //  if (Tecalc>=5000.0) {Tecalc = 5000.0;}
    //  pecabs = interpolate2D(Mne, MTe, Mg1, necalc*1.0e6, Tecalc);

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
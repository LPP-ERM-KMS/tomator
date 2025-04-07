#include "coupledpower.h"

void bfixpowerfrac_func() {
    pecabs = pecabs0;
    PIerrorP = 0.0;
    double eprof[NMESHP];
    double inteprof = 0.0;

    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof[id] = aR[id] * exp(-pow((aR[id] - Rdep) / widthech, 2.0));
        eprof[id] += echbackground * aR[id] * exp(-pow((aR[id] - Rdep) / a, 2.0)); // background
    }
    for (int id = 0; id < NMESHP - 1; ++id) {                                                                     // line integrated density
        inteprof += 2.0 * b * pi * 0.5 * (eprof[id] + eprof[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)); // 20181016
    }
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof[id] = eprof[id] / inteprof;
    }
    for (int id = 0; id < NMESHP; ++id) {
        PRFe_array[id] = pecabs * 6.24e18 * eprof[id] * (Prf * 1e3); // 20181016
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
}
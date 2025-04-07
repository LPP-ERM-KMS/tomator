#include "coupledpower.h"

void bnefix_func() {
    pecabs = pecabs0;
    nefact = (nr.ne[ic]) / (necfix);
    double nemax = 0.0;
    double eprof[NMESHP], inteprof = 0.0, inteneprof = 0.0;

    
    for (int id = 0; id < NMESHP - 1; ++id) { // line integrated density
        inteneprof += 2.0 * b * pi * 0.5 * (nr.ne[id] + nr.ne[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    }
    for (int id = 0; id < NMESHP; ++id) { // EC power deposition profile
        eprof[id] = aR[id] * exp(-pow((aR[id] - Rdep) / widthech, 2.0));
        // eprof[id] += echbackground * aR[id] * exp(-pow((aR[id] - Rdep) / a, 2.0)); // background
    }
    for (int id = 0; id < NMESHP - 1; ++id) { // line integrated power deposition profile 
        inteprof += 2.0 * b * pi * 0.5 * (eprof[id] + eprof[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)); 
    }
    for (int id = 0; id < NMESHP; ++id) { // EC power deposition profile + background profile
        eprof[id] = (1.0 - echbackground) * eprof[id] / inteprof + echbackground * nr.ne[id] / inteneprof;
    }
    
    
    int id = 0;
    while (aR[id] < Rdep)
        ++id;
    // cout << aR[id] << endl;
    double necalc = max(nr.ne[id], nr.ne[id + 1]); // cout << id << endl; //nr.ne[90]; if (necalc>1.0e13) {necalc=1.0e13;}
    double Tecalc = max(Tr.Te[id], Tr.Te[id + 1]); // Tr.Te[90];
    if (Tecalc >= 5000.0) {
        Tecalc = 5000.0;
    }

    double dTe = log10(54.8 / 50.0);
    double dpabs = log10(0.109052 / 0.100026);
    double cst = log10(50.0) * dpabs / dTe - log10(0.100026);
    pecabs = pow(10.0, log10(Tecalc) * dpabs / dTe - cst); 
    // cout << " ne = " << necalc << "    Te = " << Tecalc << "     pecabs = " << pecabs << endl;


    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        if (nr.ne[id] > nemax) {
            nemax = nr.ne[id];
        }
    }

    double PIerrorP_old = PIerrorP;

    pecabs = nefact * pecabs;
    // cout << pecabs << " " ;
    // cout << nefact << "    " << pecabs << "    " << PIerrorP <<  "    " ;
    // Bv scan, 1.54T
    PIerrorP += (1.0 - nefact) * dtnew / tauP;
    if (PIerrorP < 1.0 - (pow(1.0 / nefact, 2.0))) {
        PIerrorP = 1.0 - pow(1.0 / nefact, 2.0);
    }
    pecabs = pecabs * max(1.0, (pow(1.0 / nefact, 2.0) + PIerrorP));
    // cout << pecabs << "    " << PIerrorP <<  "    " ;
    if (pecabs > 1.0) {
        pecabs = 1.0;
        PIerrorP = PIerrorP_old;
    }
    if (pecabs < 0.001) {
        pecabs = 0.001;
        PIerrorP = PIerrorP_old;
    }

    // cout << pecabs << endl ;

    // cout << pecabs << "    " << PIerrorP << endl;
    // } else // P scan
    // {      // 143T
    //     PIerrorP += (1.0 - nefact) * dtnew / tauP;
    //     if (PIerrorP < 1.0 - (pow(1.0 / nefact, 2.0))) {
    //         PIerrorP = 1.0 - pow(1.0 / nefact, 2.0);
    //     }
    //     pecabs = 0.333 * pecabs * max(1.0, (pow(1.0 / nefact, 2.0) + PIerrorP));
    //     if (pecabs > 1.0) {
    //         pecabs = 1.0;
    //         PIerrorP = PIerrorP_old;
    //     }
    //     if (pecabs < 0.001) {
    //         pecabs = 0.001;
    //         PIerrorP = PIerrorP_old;
    //     }
    // }

    for (int id = 0; id < NMESHP; ++id) {
        PRFe_array[id] = pecabs * 6.24e18 * eprof[id] * (Prf * 1e3); // 20181016
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
}
#include "coupledpower.h"

void bproptone_func() {
    double lne = 0.0;
    double lnenue = 0.0;
    for (int id = 0; id < NMESHP - 1; ++id) { // line integrated density
        lne += 0.5 * (nr.ne[id] * aR[id] + nr.ne[id + 1] * aR[id + 1]) * (aR[id + 1] - aR[id]);
        lnenue += 0.5 * (nr.ne[id] * colrate.nue[id] * aR[id] + nr.ne[id + 1] * colrate.nue[id + 1] * aR[id + 1]) * (aR[id + 1] - aR[id]);
    }
    lnenue = lnenue / (aR[NMESHP - 1] - aR[0]);
    #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        // PRFe[id] = 6.24e18 * Prf/Vpl  * nr.ne[id] / lne;
        PRFe_array[id] = 6.24e18 * Prf / Vpl * nr.ne[id] * colrate.nue[id] * aR[id] / lnenue;
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
    //        PRFe[0] = 0.0;
    //        PRFHi[0] = 0.0;
}

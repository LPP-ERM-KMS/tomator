#include "coupledpower.h"

void bnopower_func() {
#pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        PRFe_array[id] = 0.0;
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
}
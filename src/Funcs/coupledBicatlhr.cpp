#include "coupledpower.h"

void bicatlhr_func() {
    // pecabs = pecabs0;
    PIerrorP = 0.0;
    double Prest = 1.0; // remaining power after passing lhr
    double Pnorm = 0.0; // to normalize total power to 1 after one forward and one backward pass
    double eprof[NMESHP];
    double eproftot[NMESHP];
    for (int id = 0; id < NMESHP; ++id) { 
        eprof[id]=0.0;
        eproftot[id]=0.0;
    }
    double inteprof = 0.0;
    double wrf2 = pow(2.0*pi*freq*1e6,2.0);
    double Nlhr = 0.0; // counts number LHRs
    double Rlhr[NMESHP]; // locations of LHRs
    for (int id = 0; id < NMESHP; ++id) { 
        Rlhr[id]=0.0; //ifnot zero, then lhr, value corresponds to the fraction of the power coupled at this lhr
    }
    int iAnt = 0;
    while (aR[iAnt]<Rant){
        iAnt++;
    }
    int ilHFS = NMESHP-1;
    while (aR[ilHFS]>lHFS){
        ilHFS--;
    }

    double wce2;
    double wcHi2;
    double wcDi2;
    double wcH2i2;
    double wcHDi2;
    double wcD2i2;
    double wcH3i2;
    double wcHeII2;
    double wcHeIII2;

    double wpe2;
    double wpHi2;
    double wpDi2;
    double wpH2i2;
    double wpHDi2;
    double wpD2i2;
    double wpH3i2;
    double wpHeII2;
    double wpHeIII2;

    double eps_perp[NMESHP];
    double sign_eps_perp;

    // for (int id = 0; id < NMESHP; ++id) {                                                                    
    //     wce2[id] = pow( 1.76e7*Br[id], 2.0);
    //     wcHi2[id] = pow( 9.58e3*Br[id], 2.0);
    //     wcH2i2[id] = pow( 9.58e3*Br[id]/2.0, 2.0);
    //     wcH3i2[id] = pow( 9.58e3*Br[id]/3.0, 2.0);
    //     wcHeII2[id] = pow( 9.58e3*Br[id]/4.0, 2.0);
    //     wcHeIII2[id] = pow( 9.58e3*Br[id]*2.0/4.0, 2.0);

    //     wpe2[id] = pow( 5.64e4*sqrt(nr.ne[id]),2.0);
    //     wpHi2[id] = pow( 1.32e3*sqrt(nr.nHi[id]),2.0);
    //     wpH2i2[id] = pow( 1.32e3*sqrt(nr.nH2i[id]/2.0),2.0);
    //     wpH3i2[id] = pow( 1.32e3*sqrt(nr.nH3i[id]/3.0),2.0);
    //     wpHeII2[id] = pow( 1.32e3*sqrt(nr.nHeII[id]/4.0),2.0);
    //     wpHeIII2[id] = pow( 1.32e3*2.0*sqrt(nr.nHeIII[id]/4.0),2.0);
    // }

    for (int id = 0; id < NMESHP; ++id) {                                                                    
        wce2 = pow( 1.76e7*Br[id]*1e4, 2.0);
        wcHi2 = pow( 9.58e3*Br[id]*1e4, 2.0);
        wcDi2 = pow( 9.58e3*Br[id]*1e4/2.0, 2.0);
        wcH2i2 = pow( 9.58e3*Br[id]*1e4/2.0, 2.0);
        wcHDi2 = pow( 9.58e3*Br[id]*1e4/3.0, 2.0);
        wcD2i2 = pow( 9.58e3*Br[id]*1e4/4.0, 2.0);
        wcH3i2 = pow( 9.58e3*Br[id]*1e4/3.0, 2.0);
        wcHeII2 = pow( 9.58e3*Br[id]*1e4/4.0, 2.0);
        wcHeIII2 = pow( 9.58e3*Br[id]*1e4*2.0/4.0, 2.0);

        wpe2 = pow( 5.64e4*sqrt(nr.ne[id]), 2.0);
        wpHi2 = pow( 1.32e3*sqrt(nr.nHi[id]*HtoHD), 2.0);
        wpDi2 = pow( 1.32e3*sqrt(nr.nHi[id]*(1-HtoHD)/2.0), 2.0);
        wpH2i2 = pow( 1.32e3*sqrt(nr.nH2i[id]*HtoHD*HtoHD/2.0), 2.0);
        wpHDi2 = pow( 1.32e3*sqrt(nr.nH2i[id]*2.0*(1-HtoHD)*HtoHD/3.0), 2.0);
        wpD2i2 = pow( 1.32e3*sqrt(nr.nH2i[id]*(1-HtoHD)*(1-HtoHD)/4.0), 2.0);
        wpH3i2 = pow( 1.32e3*sqrt(nr.nH3i[id]/3.0), 2.0);
        wpHeII2 = pow( 1.32e3*sqrt(nr.nHeII[id]/4.0), 2.0);
        wpHeIII2 = pow( 1.32e3*2.0*sqrt(nr.nHeIII[id]/4.0), 2.0);

        eps_perp[id] = 1.0 - wpe2 / (wrf2 - wce2)
                           - wpHi2 / (wrf2 - wcHi2) 
                           - wpDi2 / (wrf2 - wcDi2) 
                           - wpH2i2 / (wrf2 - wcH2i2) 
                           - wpHDi2 / (wrf2 - wcHDi2) 
                           - wpD2i2 / (wrf2 - wcD2i2) 
                           - wpH3i2 / (wrf2 - wcH3i2) 
                           - wpHeII2 / (wrf2 - wcHeII2)
                           - wpHeIII2 / (wrf2 - wcHeIII2) ;
    }

    // start from antenna and check sign changes for eps-perp, indicating a lhr
    sign_eps_perp = copysign(1.0,eps_perp[iAnt]);
    for (int id = iAnt-1; id > ilHFS; --id) { // forward
        // cout << iAnt << "  " << id << "  " << eps_perp[id] << "  " << copysign(1.0,eps_perp[id]) << endl;
        if (sign_eps_perp != copysign(1.0,eps_perp[id])){
            // cout << "LHR at " << aR[id] << endl;
            ++Nlhr;
            Rlhr[id]=Prest*fraclhr;
            Prest-=Prest*fraclhr;
        }
        sign_eps_perp = copysign(1.0,eps_perp[id]);
    }
    for (int id = ilHFS; id <= iAnt; ++id) { // backward
        if (Rlhr[id]!=0.0){ // then lhr
            Rlhr[id]+=Prest*fraclhr;
            Prest-=Prest*fraclhr;
            Pnorm+=Rlhr[id];
        }
    }
    for (int id = ilHFS; id <= iAnt; ++id) { // normalise to 1
        if (Rlhr[id]!=0.0){ // then lhr
            Rlhr[id]=Rlhr[id]/Pnorm;
        }
    }

    // if no lhr, then couple power at antenna
    if (Nlhr==0.0){
        Nlhr = 1.0;
        Rlhr[iAnt]=1.0;
    }
    // cout << Nlhr << endl;
    // start again from antenna, put power at lhrs
    for (int idlhr = iAnt; idlhr > ilHFS; --idlhr) { 
        if (Rlhr[idlhr]!=0.0){ // then lhr
            // cout << aR[idlhr] << endl;
            for (int id = 0; id < NMESHP; ++id) { // predefined profile at lhr
                eprof[id] = aR[id] * exp(-pow((aR[id] - aR[idlhr]) / widthlhr, 2.0));
                eprof[id] += lhrbackground * aR[id] * exp(-pow((aR[id] - aR[idlhr]) / a, 2.0)); // background
            }
            inteprof=0.0;
            for (int id = 0; id < NMESHP - 1; ++id) { // integrate profile
                inteprof += 2.0 * b * pi * 0.5 * (eprof[id] + eprof[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
            }
            for (int id = 0; id < NMESHP; ++id) { 
                eprof[id] = eprof[id] / inteprof * Rlhr[idlhr]; /// Nlhr; // normalize coupled power to Nlhr
                eproftot[id] += eprof[id];  // add to total coupled power profile
                // cout << aR[id] << "  " << eproftot[id] << endl;
            }
        }
    }

    for (int id = 0; id < NMESHP; ++id) { // put fracpne% on the electrons and the rest on ions
        PRFe_array[id] = 6.24e18 * eproftot[id] * (Prf * 1e3);
        PRFHi_array[id] = (1-fracpne) * PRFe_array[id] * nr.nHi[id]/(nr.nHi[id]+nr.nH2i[id]+nr.nH3i[id]+nr.nHeII[id]+nr.nHeIII[id]);
        PRFH2i_array[id] = (1-fracpne) * PRFe_array[id] * nr.nH2i[id]/(nr.nHi[id]+nr.nH2i[id]+nr.nH3i[id]+nr.nHeII[id]+nr.nHeIII[id]);
        PRFH3i_array[id] = (1-fracpne) * PRFe_array[id] * nr.nH3i[id]/(nr.nHi[id]+nr.nH2i[id]+nr.nH3i[id]+nr.nHeII[id]+nr.nHeIII[id]);
        PRFHeII_array[id] = (1-fracpne) * PRFe_array[id] * nr.nHeII[id]/(nr.nHi[id]+nr.nH2i[id]+nr.nH3i[id]+nr.nHeII[id]+nr.nHeIII[id]);
        PRFHeIII_array[id] = (1-fracpne) * PRFe_array[id] * nr.nHeIII[id]/(nr.nHi[id]+nr.nH2i[id]+nr.nH3i[id]+nr.nHeII[id]+nr.nHeIII[id]);
        PRFe_array[id] = fracpne * PRFe_array[id];
        //cout << aR[id] << "  " << PRFe_array[id] << endl;
    }
    
}
#include "coupledpower.h"

void bTOMAS_func() {
    double eprof1[NMESHP];
    double inteprof1 = 0.0;
    double eprof2[NMESHP];
    double inteprof2 = 0.0;
    double eprof3[NMESHP];
    double inteprof3 = 0.0;
    double wuh[NMESHP];
    double wRcut[NMESHP];
    double wLcut[NMESHP];
    // for (int id=0; id<NMESHP; ++id)
    // 	{PTEST[id]=0;}

    double nemax = nr.ne[0];
    double Temax = Tr.Te[0];
    for (int id = 1; id < NMESHP - 1; ++id) { // line integrated density
        nemax = max(nemax, nr.ne[id]);
        Temax = max(Temax, Tr.Te[id]);
    }

    //------------------//
    // Locate ECR-layer //
    //------------------//

    double Bce = 2.0 * pi * freq * 1e6 * me / qe / harmonic; // value of magnetic field for EC resonance
    int id = 0;
    while (Br[id] > Bce)
        ++id;
    // cout << "id= " << id << endl;
    if ((abs(Br[id] - Bce)) > (abs(Br[id + 1] - Bce))) {
        ++id;
    }
    // cout << "idcorr= " << id << endl;
    double Rdepce = aR[id];
    // double Rdepceid = id; // Meshpoint where ECR is located
    double Teabsce = Tr.Te[id];
    double neabsce = nr.ne[id];
    double wpe = sqrt(neabsce * 1e6 * qe * qe / me / eps0);
    double wce = 2.0 * pi * freq * 1e6 / harmonic;
    double alphaCE = wpe * wpe / wce / wce; // wpe^2/wce^2 @ CE

    //-----------------//
    // Locate R-cutoff //
    //-----------------//

    for (int id = 0; id < NMESHP; ++id) {
        wRcut[id] = 0.5 * Br[id] * qe / me + 0.5 * sqrt(4.0 * nr.ne[id] * 1e6 * qe * qe / me / eps0 + pow(Br[id] * qe / me, 2)); // frequency for UH resonance
    }
    id = 0;
    while (wRcut[id] > wce)
        ++id;

    if ((abs(wRcut[id] - wce)) > (abs(wRcut[id + 1] - wce))) {
        ++id;
    }
    // double RRcut = aR[id]; // position of R-cutoff
    //  double DensRcut = abs(((4*pow(wce-0.5*qe*Bt/me,2)*me*me*eps0/qe/qe)-Bt*Bt*eps0)/4.0/me*1e-6); // density for R-cutoff [cm-3]
    double DensRcut = nr.ne[id];
    // cout << "DensRcut= " << DensRcut << "RRcut= " << RRcut << "Rdepce= " << Rdepce    << endl;

    //-----------------//
    // Locate L-cutoff //
    //-----------------//

    for (int id = 0; id < NMESHP; ++id) {
        wLcut[id] = -0.5 * Br[id] * qe / me + 0.5 * sqrt(4.0 * nr.ne[id] * 1e6 * qe * qe / me / eps0 + pow(Br[id] * qe / me, 2)); // frequency for UH resonance
    }
    id = NMESHP / 2.0 - 1.0;
    while (wLcut[id] > wce)
        ++id;
    if ((abs(wLcut[id] - wce)) > (abs(wLcut[id + 1] - wce))) {
        ++id;
    }
    // double RLcut = aR[id]; // position of L-cutoff
    double DensLcut = nr.ne[id];
    // double DensLcut = abs(((4*pow(wce-0.5*qe*Bt/me,2)*me*me*eps0/qe/qe)+Bt*Bt*eps0)/4.0/me*1e-6); // density for R-cutoff [cm-3]

    //------------------//
    // Locate UHR-layer //
    //------------------//

    double Ruh = 0;
    for (int id = 0; id < NMESHP; ++id) {
        wuh[id] = sqrt(nr.ne[id] * 1e6 * qe * qe / me / eps0 + pow(Br[id] * qe / me, 2)); // frequency for UH resonance
    }
    id = 0;
    while (wuh[id] > wce)
        ++id;
    if ((abs(wuh[id] - wce)) > (abs(wuh[id + 1] - wce))) {
        ++id;
    }
    Ruh = aR[id];
    double Densuh = nr.ne[id];
    // if (aR[id] > R) {Ruh = aR[id];}   // position of UHR-layer if present
    // if (aR[id] > R+a) {Ruh = R+a;}
    // double Densuh = abs((me*eps0/qe/qe*wce*wce - Bt*Bt*eps0/me)*1e-6);		 // density for UHR [cm-3]
    // cout << "Ruh = " << Ruh << "      Densuh= " << Densuh  << endl;

    //-----------------------//
    // Budden Analysis @ UHR //
    //-----------------------//

    double Lne = abs(nr.ne[id] * abs(aR[id + 1] - aR[id - 1]) / 100.0 / (nr.ne[id - 1] - nr.ne[id + 1]));                                                                                                                                     // density scale length [m]
    double LBr = abs(Br[id] * abs(aR[id + 1] - aR[id - 1]) / 100.0 / (Br[id - 1] - Br[id + 1]));                                                                                                                                              // magnetic scal length [m]
    double alphaUHR = sqrt(nr.ne[id] * 1.0e6 * qe * qe / me / eps0) / (Br[id] * qe / me);                                                                                                                                                     // wpe/wce @ UHR
    double Buddenparam = (Br[id] * qe / me) * Lne / c * alphaUHR / sqrt(alphaUHR * alphaUHR + 2.0 * (Lne / LBr)) * pow((sqrt(1.0 + alphaUHR * alphaUHR) - 1.0) / (alphaUHR * alphaUHR + (Lne / LBr) * sqrt(1.0 + alphaUHR * alphaUHR)), 0.5); // Budden parameter
    double Tunnel = exp(-2.0 * Buddenparam);                                                                                                                                                                                                  // Tunneling fraction
    double Reflect = pow((1.0 - Tunnel), 2.0);                                                                                                                                                                                                // Reflected part
    double Convert = 1.0 - Tunnel - Reflect;                                                                                                                                                                                                  // XB converted fraction
    // cout << "T= " << Tunnel << "      R= " << Reflect << "      C= " << Convert << endl;

    //--------------//
    // O absorption //
    //--------------//

    double DensO = wce * wce * me * eps0 / qe / qe * 1e-6;
    double PabsO = (1.0 - exp(-0.25 * pi * R / 100.0 / c / c / c * Teabsce * 11604.45 * kb / me * wpe * wpe / wce * sqrt(1.0 - alphaCE))); // 1 pass aborption for O-wave
    PabsO = PabsO / (1.0 - (1.0 - muw) * (1.0 - PabsO));                                                                                   // multi-pass absorption with losses to the wall (muw)
    if (nemax > DensO) {
        PabsO = 0.0;
    } //
    // cout << "      DensO= " << DensO  << endl;

    //--------------//
    // X absorption //
    //--------------//

    double PabsX = (1.0 - exp(-0.5 * pi * R / 100.0 / c / c / c * Teabsce * 11604.45 * kb / me * wpe * wpe / wce * pow(2.0 - alphaCE, 2.5) * (1.0 + sqrt(alphaCE)))); // * (DensRcut*me/eps0/Bce/Bce-1e-3) / (DensRcut*me/eps0/Bce/Bce) ));  // 1 pass aborption for X-wave below Rcutoff
    DensRcut = 3.13e9;
    Densuh = 6.26e9;
    if (nemax < DensRcut) {
        PabsX = PabsX * abs(DensRcut * 1.0e6 * me / eps0 / Bce / Bce - alphaCE - 1.0e-5) / (DensRcut * 1.0e6 * me / eps0 / Bce / Bce - alphaCE + 1.0e-5); // adaptation for X-wave nearing Rcutoff
        PabsX = PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX));                                                                                              // multi-pass absorption with losses to the wall (muw)
    } else if (nemax < Densuh) {
        // PabsX = (1.0-Reflect)*PabsX + (1.0-Reflect)*(1.0-PabsX)*(1.0-muw)*PabsX;  //Not reflected part tunnels and reaches ECR,, reflecting on wall and coupling at second pass
        // PabsX = PabsX * abs(DensRcut*1.0e6*me/eps0/Bce/Bce-alphaCE-1.0e-5) / (DensRcut*1.0e6*me/eps0/Bce/Bce-alphaCE+1.0e-5) ;  // adaptation for X-wave nearing Rcutoff

        PabsX = (1.0 - Reflect) * PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX)); // multi-pass absorption of the tunneling wave with losses to the wall (muw)
        PabsX = PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX));                   // multi-pass absorption with losses to the wall (muw)
    }
    // Not reflected part tunnels and reaches ECR,, reflecting on wall and coupling at second pass
    else if (nemax < DensLcut) {
        PabsX = (1.0 - PRdep2) * Tunnel * PabsX + (1.0 - PRdep4) * Tunnel * (1.0 - PabsX) * (1.0 - muw) * PabsX;
    }
    // X-wave tunneling through R-cutoff, reflecting on wall and coupling at second pass
    // third pass not possible because of evanescent region, SX-B conversion possible
    else {
        PabsX = 0.0;
    } // All X-wave is converted to B-wave

    //--------------//
    // B absorption //
    //--------------//

    double PabsB = 0.0; // EBW coupling near ECR
    id = 0;
    while (Br[id] > Bce)
        ++id;
    if ((abs(Br[id] - Bce)) > (abs(Br[id + 1] - Bce))) {
        ++id;
    }
    // double RdepB = aR[id]; //no shift to LFS of ECR implemented for the moment, doppler shift to be added ????
    double RdepB = Rdep3; // Hardcoded position of B-coupling

    if (nemax < Densuh) {
        PabsB = 0.0;
    }
    // No B-wave present

    else if (nemax < DensO) {
        PabsB = PRdep2 * Convert * (1.0 - exp(-0.5 * pi * R / 100.0 / c * wce * alphaCE)); // percentage of converted X-wave to B-wave coupling
    }
    // B-wave appears when UHR appears
    else {
        PabsB = PRdep3 * Convert * (1.0 - exp(-0.5 * pi * R / 100.0 / c * wce * alphaCE));
        PabsB += PRdep3 * (1.0 - PRdep4) * (1.0 - PRdep4) * Tunnel;
    }
    // double conversion to B-wave present: 1st converted at UHR, 2nd from tunneling X-wave reflected at Lcutoff and converted to B-wave at UHR

    //--------------//
    // UH Resonance //
    //--------------//

    double PabsU = 0.0; // coupling at the UHR

    if (nemax < Densuh) {
        PabsU = 0.0;
    }
    // No UHR present
    else if (nemax < DensO) {
        PabsU = PRdep4 * (Tunnel); // PRdep4 = percentage of UHR coupling of tunneled X-wave
    } else {
        PabsU = PRdep4 * (Tunnel);                   // 1st resonance at UHR
        PabsU += PRdep4 * (1.0 - PRdep4) * (Tunnel); // 2nd resonance at UHR from tunneled X-wave reflected at Lcutoff
    }

    // double Rdepvar = 0.0;
    //  if (nemax < 1.4e10) {Rdepvar=Rdepce-7.0;}
    //   else {Rdepvar=Rdepce+10.0e-10*(nemax-2.1e10);}

    // cout << "Rdep= " << Rdepvar << endl;

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        //    eprof[id] = PRdep1*aR[id]*pow(1+pow((aR[id]-Rdep1)/(pow(2.0*a/100.0,0.5)*0.5),2.0),-(0.5+1.0)); // 20181016
        //	eprof[id] = PRdep1*(nemax*Temax/2.5e10)*aR[id]*pow(1.0+pow((aR[id]-Rdep1)/(pow(2.0*a/100.0,0.5)*0.5),2.0),-(0.5+1.0));

        eprof1[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - Rdepce) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0)); // normaalverdeling met hoogte 1 voor O/X-wave deposition at ECR
        eprof2[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - RdepB) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0));  // B-wave deposition at ECR+1.5cm (should be moving from 0 to + 2.0 ???)
        eprof3[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - Ruh) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0));    // UH Resonance
                                                                                                                               // eprof3[id]  = 0.01*aR[id]*pow(1.0+pow((aR[id]-RColl)/(pow(2.0*a,0.5)*0.5),2.0),-(0.5+1.0)); // collisional damping (broader profile)
    }

    for (int id = 0; id < NMESHP - 1; ++id) {                                                                        // line integrated density
        inteprof1 += 2.0 * b * pi * 0.5 * (eprof1[id] + eprof1[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)); // 20181016
        inteprof2 += 2.0 * b * pi * 0.5 * (eprof2[id] + eprof2[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)); // 20181016
        inteprof3 += 2.0 * b * pi * 0.5 * (eprof3[id] + eprof3[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)); // 20181016
    }
    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof1[id] = eprof1[id] / inteprof1;
        eprof2[id] = eprof2[id] / inteprof2;
        eprof3[id] = eprof3[id] / inteprof3;

        // cout << id << "  " << eprof[id]*(pow(aR[NMESHP-1],2.0)-pow(aR[0],2.0))/(pow(aR[id+1],2.0)-pow(aR[id],2.0)) << endl;
    }
    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        // PRFe[id] =  6.24e18 * (PRdep1*Pabs*eprof1[id]+ PRdep2*eprof2[id] +PRdep3*eprof3[id] ) *(Prf*1e3);  // 20181016

        Pabs1[id] = (PRdep1 * PabsO * eprof1[id]) * (Prf * 1.0e3); // Power (kW) of O-wave absorption at ECR
        Pabs2[id] = (PRdep2 * PabsX * eprof1[id]) * (Prf * 1.0e3); // Power (kW) of X-wave absorption at ECR (inluding tunneled)
        Pabs3[id] = (PRdep2 * PabsB * eprof2[id]) * (Prf * 1.0e3); // Power (kW) of B-wave absorption near ECR (converted from X-wave => PRdep2)
        Pabs4[id] = (PRdep2 * PabsU * eprof3[id]) * (Prf * 1.0e3); // Power (kW) of UHR absorption at UHR
        PRFe_array[id] = 6.24e18 * Pabs1[id] + 6.24e18 * Pabs2[id] + 6.24e18 * Pabs3[id] + 6.24e18 * Pabs4[id];
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
    // PTEST[0]=RRcut;
    // PTEST[1]=DensRcut;
    // PTEST[2]=Ruh;
    // PTEST[3]=Densuh;
    // PTEST[4]=RLcut;
    // PTEST[5]=DensLcut;
    // PTEST[6]=nemax;
    // PTEST[7]=Temax;
    // PTEST[8]=Reflect;
    // PTEST[9]=Tunnel;
    // PTEST[10]=Convert;
    //
    // PTEST[11]=Rdepce;
    // PTEST[12]=RdepB;
    // PTEST[13]=PabsO;
    // PTEST[14]=PabsX;
    // PTEST[15]=PabsB;
    // PTEST[16]=PabsU;
}
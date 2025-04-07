#include "timeStep.h"

void timeStep() {
    double daccur = accur * 1.1;
    int success = 0;
    int TimeStepCounter = 0;
    double drmaxold = 0;
    double drmax4old = 0;
    bool dtfix = false;

    double w = 0.0;

    double dth1;
    double dth2;
    double neh;
    double tempr = 0.0;
    int counter = 0;

    //double tfact = (1.0 - 0.99 * exp(-tmain / 1e-4));
    double tfact = 1.0;

    while (success == 0) {
        counter++;
        drmaxold = drmax;
        drmax4old = drmax4;
        // cout << " dtnew " << dtnew << endl ;

        w = dtnew / oldtstep; // REFERENCE work Implementation of an Adaptive BDF2 Formula and Comparison with the MATLAB Ode15s

        if (Nit > 10) {
            coef1 = (1 + w) / (1 + 2 * w);
            coef2 = pow((1 + w), 2) / (1 + 2 * w);
            coef3 = pow(w, 2) / (1 + 2 * w);
        }
        else {
            coef1 = 2.0 / 3.0;
            coef2 = 4.0 / 3.0;
            coef3 = 1.0 / 3.0;
        }

        // #include "../Phys/solve_Ab=x_v3.03.h" // solve a first time with the previous time step
        solverAb_x();

        drval = 0.0;
        drmax = 0.0;
        for (int im = 0; im < NMESHP; ++im) { // The change of a state cannot be larger than accur
            drval = (dnr_cn.dne[im]) / nr.ne[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnH[im]) / nr.nH[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnH2[im]) / nr.nH2[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnHi[im]) / nr.nHi[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnH2i[im]) / nr.nH2i[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnH3i[im]) / nr.nH3i[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnHeI[im]) / nr.nHeI[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnHeII[im]) / nr.nHeII[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dnr_cn.dnHeIII[im]) / nr.nHeIII[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            if (bimpur) {
                drval = (dnr_cn.dnCI[im]) / nr.nCI[im];
                if (abs(drmax) < abs(drval)) {
                    drmax = drval;
                }
                drval = (dnr_cn.dnCII[im]) / nr.nCII[im];
                if (abs(drmax) < abs(drval)) {
                    drmax = drval;
                }
                drval = (dnr_cn.dnCIII[im]) / nr.nCIII[im];
                if (abs(drmax) < abs(drval)) {
                    drmax = drval;
                }
                drval = (dnr_cn.dnCIV[im]) / nr.nCIV[im];
                if (abs(drmax) < abs(drval)) {
                    drmax = drval;
                }
                drval = (dnr_cn.dnCV[im]) / nr.nCV[im];
                if (abs(drmax) < abs(drval)) {
                    drmax = drval;
                }
            }
            drval = (dEr_cn.dEe[im]) / Er.Ee[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEH[im]) / Er.EH[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEH2[im]) / Er.EH2[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEHi[im]) / Er.EHi[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEH2i[im]) / Er.EH2i[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEH3i[im]) / Er.EH3i[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEHeI[im]) / Er.EHeI[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEHeII[im]) / Er.EHeII[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
            drval = (dEr_cn.dEHeIII[im]) / Er.EHeIII[im];
            if (abs(drmax) < abs(drval)) {
                drmax = drval;
            }
        }
        drval2 = 0.0;
        drmax2 = 0.0;
        for (int im = 1; im < NMESHP - 1; ++im) { // The change of difference between two neighbouring states cannot be larger than accur
            double dr;
            if (im != NMESHP - 1) {
                dr = aR[im + 1] - aR[im];
            }
            else {
                dr = aR[NMESHP - 1] - aR[NMESHP - 2];
            } // cout << dr << endl;
            drval2 = abs((dnr_cn.dxne[im] * dr / nr.ne[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   ne   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnH[im] * dr / nr.nH[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nH   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnH2[im] * dr / nr.nH2[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nH2   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnHi[im] * dr / nr.nHi[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nHi   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnH2i[im] * dr / nr.nH2i[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nH2i   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnH3i[im] * dr / nr.nH3i[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nH3i   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnHeI[im] * dr / nr.nHeI[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nHeI   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnHeII[im] * dr / nr.nHeII[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   nHeII   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dnr_cn.dxnHeIII[im] * dr / nr.nHeIII[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2;
            }
            if (bimpur) {
                drval2 = abs((dnr_cn.dxnCI[im] * dr / nr.nCI[im]));
                if (abs(drmax2) < abs(drval2)) {
                    drmax2 = drval2; /*cout << "   nCI   dt = " << accur/abs(drmax2) << endl; */
                }
                drval2 = abs((dnr_cn.dxnCII[im] * dr / nr.nCII[im]));
                if (abs(drmax2) < abs(drval2)) {
                    drmax2 = drval2; /*cout << "   nCII   dt = " << accur/abs(drmax2) << endl; */
                }
                drval2 = abs((dnr_cn.dxnCIII[im] * dr / nr.nCIII[im]));
                if (abs(drmax2) < abs(drval2)) {
                    drmax2 = drval2; /*cout << "   nCIII   dt = " << accur/abs(drmax2) << endl; */
                }
                drval2 = abs((dnr_cn.dxnCIV[im] * dr / nr.nCIV[im]));
                if (abs(drmax2) < abs(drval2)) {
                    drmax2 = drval2; /*cout << "   nCIII   dt = " << accur/abs(drmax2) << endl; */
                }
                drval2 = abs((dnr_cn.dxnCV[im] * dr / nr.nCV[im]));
                if (abs(drmax2) < abs(drval2)) {
                    drmax2 = drval2; /*cout << "   nCIII   dt = " << accur/abs(drmax2) << endl; */
                }
            }
            drval2 = abs((dEr_cn.dxEe[im] * dr / Er.Ee[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   Ee   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEH[im] * dr / Er.EH[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EH   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEH2[im] * dr / Er.EH2[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EH2   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEHi[im] * dr / Er.EHi[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EHi   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEH2i[im] * dr / Er.EH2i[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EH2i   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEH3i[im] * dr / Er.EH3i[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EH3i   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEHeI[im] * dr / Er.EHeI[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EHeI   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEHeII[im] * dr / Er.EHeII[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EHeII   dt = " << accur/abs(drmax2) << endl; */
            }
            drval2 = abs((dEr_cn.dxEHeIII[im] * dr / Er.EHeIII[im]));
            if (abs(drmax2) < abs(drval2)) {
                drmax2 = drval2; /*cout << "   EHeIII   dt = " << accur/abs(drmax2) << endl; */
            }
        }
        drval4 = 0.0;
        drmax4 = 0.0;
        for (int im = 0; im < NMESHP; ++im) {
            double dr;
            if (im != NMESHP - 1) {
                dr = aR[im + 1] - aR[im];
            }
            else {
                dr = aR[NMESHP - 1] - aR[NMESHP - 2];
            } // cout << dr << endl;
            drval4 = abs(asin(dnr_cn.dxne[im] * dr * dtnew / nr.ne[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   ne   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnH[im] * dr * dtnew / nr.nH[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nH   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnH2[im] * dr * dtnew / nr.nH2[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nH2   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnHi[im] * dr * dtnew / nr.nHi[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nHi   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnH2i[im] * dr * dtnew / nr.nH2i[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nH2i   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnH3i[im] * dr * dtnew / nr.nH3i[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nH3i   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnHeI[im] * dr * dtnew / nr.nHeI[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nHeI   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnHeII[im] * dr * dtnew / nr.nHeII[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   nHeII   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dnr_cn.dxnHeIII[im] * dr * dtnew / nr.nHeIII[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4;
            }
            if (bimpur) {
                drval4 = abs(asin(dnr_cn.dxnCI[im] * dr * dtnew / nr.nCI[im])) / (pi / 2.0) / dtnew;
                if (abs(drmax4) < abs(drval4)) {
                    drmax4 = drval4; /*cout << "   nCI   dt = " << accur/abs(drmax4) << endl; */
                }
                drval4 = abs(asin(dnr_cn.dxnCII[im] * dr * dtnew / nr.nCII[im])) / (pi / 2.0) / dtnew;
                if (abs(drmax4) < abs(drval4)) {
                    drmax4 = drval4; /*cout << "   nCII   dt = " << accur/abs(drmax4) << endl; */
                }
                drval4 = abs(asin(dnr_cn.dxnCIII[im] * dr * dtnew / nr.nCIII[im])) / (pi / 2.0) / dtnew;
                if (abs(drmax4) < abs(drval4)) {
                    drmax4 = drval4; /*cout << "   nCIII   dt = " << accur/abs(drmax4) << endl; */
                }
                drval4 = abs(asin(dnr_cn.dxnCIV[im] * dr * dtnew / nr.nCIV[im])) / (pi / 2.0) / dtnew;
                if (abs(drmax4) < abs(drval4)) {
                    drmax4 = drval4; /*cout << "   nCIII   dt = " << accur/abs(drmax4) << endl; */
                }
                drval4 = abs(asin(dnr_cn.dxnCV[im] * dr * dtnew / nr.nCV[im])) / (pi / 2.0) / dtnew;
                if (abs(drmax4) < abs(drval4)) {
                    drmax4 = drval4; /*cout << "   nCIII   dt = " << accur/abs(drmax4) << endl; */
                }
            }
            drval4 = abs(asin(dEr_cn.dxEe[im] * dr * dtnew / Er.Ee[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   Ee   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEH[im] * dr * dtnew / Er.EH[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EH   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEH2[im] * dr * dtnew / Er.EH2[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EH2   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEHi[im] * dr * dtnew / Er.EHi[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EHi   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEH2i[im] * dr * dtnew / Er.EH2i[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EH2i   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEH3i[im] * dr * dtnew / Er.EH3i[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EH3i   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEHeI[im] * dr * dtnew / Er.EHeI[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EHeI   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEHeII[im] * dr * dtnew / Er.EHeII[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EHeII   dt = " << accur/abs(drmax4) << endl; */
            }
            drval4 = abs(asin(dEr_cn.dxEHeIII[im] * dr * dtnew / Er.EHeIII[im])) / (pi / 2.0) / dtnew;
            if (abs(drmax4) < abs(drval4)) {
                drmax4 = drval4; /*cout << "   EHeIII   dt = " << accur/abs(drmax4) << endl; */
            }
        }
        drmax4 = drmax4 / 2.0;
        if (dtnew > dtmin * tfact) {
            daccur = max(abs(drmax), abs(drmax4)) * dtnew;
        }
        else {
            daccur = accur / 1.1;
        }

        if (dtfix) {
            dt = oldtstep;
            daccur = accur / 1.1;
        }
        else {
            dt = accur / max(abs(drmax), abs(drmax4)) * 0.5; // ideal time step for this loop
            if (TimeStepCounter > 5) {
                dt *= 0.5;
            }
        }

        if (dtsmooth && !dtfix) {
            if (max(abs(drmax), abs(drmax4)) >= max(abs(drmaxold), abs(drmax4old)) && TimeStepCounter > 0 && daccur > accur) {
                dtfix = true;
                dt = 0.1 * oldtstep;
                // cout << "DISCONTINUITY!!!" << endl;
                // cout << "dt will be " << dt << endl;
            }
        }

        if (daccur < accur) { // calculation was accurate enough
            success = 1;
            if (oldtstep / dtnew * daccur / olddaccur > shokparam && TimeStepCounter == 0 && dtnew > oldtstep) {
                dt = 0.1 * min(dt, dtnew);
                // cout << "!!!!!!!!!!!!!!shock incoming...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
            }
            if (dt > dtnew) {
                dt = min(dt, maxtstepincrement * dtnew);
                // dt = (dt + dtnew) / 2.0;
            }
        }
        else { // calculation was not accurate enough
            success = 0;
            /// cout << " new time step " << endl;
        }

        if (success == 1) {
            nr_m1 = nr;
            Er_m1 = Er;
// dt_m1 = dtnew;
//#pragma omp parallel for private(tempr, dth1, dth2, neh)
            for (int im = 0; im < NMESHP; ++im) {
                // Do time step, take smaller step if accur is not good enough.

                if (((abs(dnr_cn.dne[im] * dtnew) / nr.ne[im]) <= accur) | ((abs(dEr_cn.dEe[im] * dtnew) / Er.Ee[im]) <= accur)) {
                    nr.ne[im] += (dnr_cn.dne[im]) * dtnew;

                    if (fixedTe) {
                        Er.Ee[im] += 3.0 / 2.0 * dnr_cn.dne[im] * Tr.Te[im] * dtnew;
                    }
                    else {
                        Er.Ee[im] += (dEr_cn.dEe[im]) * dtnew;
                    }
                    nr.xne[im] += (dnr_cn.dxne[im]) * dtnew;
                    if (fixedTe) {
                        Er.xEe[im] += 3.0 / 2.0 * dnr_cn.dxne[im] * Tr.Te[im] * dtnew;
                    }
                    else {
                        Er.xEe[im] += (dEr_cn.dxEe[im]) * dtnew;
                    }
                }
                else {
                    dth1 = (nr.ne[im] * accur) / abs(dnr_cn.dne[im] * dtnew) * dtnew;
                    dth2 = (Er.Ee[im] * accur) / abs(dEr_cn.dEe[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.ne[im] += (dnr_cn.dne[im]) * dth1;
                    if (fixedTe) {
                        Er.Ee[im] += 3.0 / 2.0 * dnr_cn.dne[im] * Tr.Te[im] * dth1;
                    }
                    else {
                        Er.Ee[im] += (dEr_cn.dEe[im]) * dth1;
                    }
                    nr.xne[im] += (dnr_cn.dxne[im]) * dth1;
                    if (fixedTe) {
                        Er.xEe[im] += 3.0 / 2.0 * dnr_cn.dxne[im] * Tr.Te[im] * dth1;
                    }
                    else {
                        Er.xEe[im] += (dEr_cn.dxEe[im]) * dth1;
                    }
                }
                if (((abs(dnr_cn.dnH[im] * dtnew) / nr.nH[im]) <= accur) && ((abs(dEr_cn.dEH[im] * dtnew) / Er.EH[im]) <= accur)) {
                    tempr = nr.nH[im];
                    tempr += (dnr_cn.dnH[im]) * dtnew;
                    if (tempr < 0.0) {
                        tempr = nr.nH[im];
                        tempr = nr.nH[im] + (dnr_cn.dnH[im]) * dtnew;
                        nr.nH[im] = tempr;
                    }
                    nr.nH[im] += (dnr_cn.dnH[im]) * dtnew;
                    Er.EH[im] += (dEr_cn.dEH[im]) * dtnew;
                    nr.xnH[im] += (dnr_cn.dxnH[im]) * dtnew;
                    Er.xEH[im] += (dEr_cn.dxEH[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nH[im] * accur) / abs(dnr_cn.dnH[im] * dtnew) * dtnew;
                    dth2 = (Er.EH[im] * accur) / abs(dEr_cn.dEH[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    tempr = nr.nH[im];
                    tempr += (dnr_cn.dnH[im]) * dth1;
                    if (tempr < 0.0) {
                        tempr = nr.nH[im];
                        tempr = nr.nH[im] + (dnr_cn.dnH[im]) * dth1;
                        nr.nH[im] = tempr;
                    }
                    nr.nH[im] += (dnr_cn.dnH[im]) * dtnew;
                    Er.EH[im] += (dEr_cn.dEH[im]) * dth1;
                    nr.xnH[im] += (dnr_cn.dxnH[im]) * dth1;
                    Er.xEH[im] += (dEr_cn.dxEH[im]) * dth1;
                }

                if (((abs(dnr_cn.dnH2[im] * dtnew) / nr.nH2[im]) <= accur) | ((abs(dEr_cn.dEH2[im] * dtnew) / Er.EH2[im]) <= accur)) {
                    nr.nH2[im] += (dnr_cn.dnH2[im]) * dtnew;
                    Er.EH2[im] += (dEr_cn.dEH2[im]) * dtnew;
                    nr.xnH2[im] += (dnr_cn.dxnH2[im]) * dtnew;
                    Er.xEH2[im] += (dEr_cn.dxEH2[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nH2[im] * accur) / abs(dnr_cn.dnH2[im] * dtnew) * dtnew;
                    dth2 = (Er.EH2[im] * accur) / abs(dEr_cn.dEH2[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nH2[im] += (dnr_cn.dnH2[im]) * dth1;
                    Er.EH2[im] += (dEr_cn.dEH2[im]) * dth1;
                    nr.xnH2[im] += (dnr_cn.dxnH2[im]) * dth1;
                    Er.xEH2[im] += (dEr_cn.dxEH2[im]) * dth1;
                }

                if (((abs(dnr_cn.dnHi[im] * dtnew) / nr.nHi[im]) <= accur) | ((abs(dEr_cn.dEHi[im] * dtnew) / Er.EHi[im]) <= accur)) {
                    nr.nHi[im] += (dnr_cn.dnHi[im]) * dtnew;
                    Er.EHi[im] += (dEr_cn.dEHi[im]) * dtnew;
                    nr.xnHi[im] += (dnr_cn.dxnHi[im]) * dtnew;
                    Er.xEHi[im] += (dEr_cn.dxEHi[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nHi[im] * accur) / abs(dnr_cn.dnHi[im] * dtnew) * dtnew;
                    dth2 = (Er.EHi[im] * accur) / abs(dEr_cn.dEHi[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nHi[im] += (dnr_cn.dnHi[im]) * dth1;
                    Er.EHi[im] += (dEr_cn.dEHi[im]) * dth1;
                    nr.xnHi[im] += (dnr_cn.dxnHi[im]) * dth1;
                    Er.xEHi[im] += (dEr_cn.dxEHi[im]) * dth1;
                }

                if (((abs(dnr_cn.dnH2i[im] * dtnew) / nr.nH2i[im]) <= accur) | ((abs(dEr_cn.dEH2i[im] * dtnew) / Er.EH2i[im]) <= accur)) {
                    nr.nH2i[im] += (dnr_cn.dnH2i[im]) * dtnew;
                    Er.EH2i[im] += (dEr_cn.dEH2i[im]) * dtnew;
                    nr.xnH2i[im] += (dnr_cn.dxnH2i[im]) * dtnew;
                    Er.xEH2i[im] += (dEr_cn.dxEH2i[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nH2i[im] * accur) / abs(dnr_cn.dnH2i[im] * dtnew) * dtnew;
                    dth2 = (Er.EH2i[im] * accur) / abs(dEr_cn.dEH2i[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nH2i[im] += (dnr_cn.dnH2i[im]) * dth1;
                    Er.EH2i[im] += (dEr_cn.dEH2i[im]) * dth1;
                    nr.xnH2i[im] += (dnr_cn.dxnH2i[im]) * dth1;
                    Er.xEH2i[im] += (dEr_cn.dxEH2i[im]) * dth1;
                }

                if (((abs(dnr_cn.dnH3i[im] * dtnew) / nr.nH3i[im]) <= accur) | ((abs(dEr_cn.dEH3i[im] * dtnew) / Er.EH3i[im]) <= accur)) {
                    nr.nH3i[im] += (dnr_cn.dnH3i[im]) * dtnew;
                    Er.EH3i[im] += (dEr_cn.dEH3i[im]) * dtnew;
                    nr.xnH3i[im] += (dnr_cn.dxnH3i[im]) * dtnew;
                    Er.xEH3i[im] += (dEr_cn.dxEH3i[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nH3i[im] * accur) / abs(dnr_cn.dnH3i[im] * dtnew) * dtnew;
                    dth2 = (Er.EH3i[im] * accur) / abs(dEr_cn.dEH3i[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nH3i[im] += (dnr_cn.dnH3i[im]) * dth1;
                    Er.EH3i[im] += (dEr_cn.dEH3i[im]) * dth1;
                    nr.xnH3i[im] += (dnr_cn.dxnH3i[im]) * dth1;
                    Er.xEH3i[im] += (dEr_cn.dxEH3i[im]) * dth1;
                }

                if (((abs(dnr_cn.dnHeI[im] * dtnew) / nr.nHeI[im]) <= accur) | ((abs(dEr_cn.dEHeI[im] * dtnew) / Er.EHeI[im]) <= accur)) {
                    nr.nHeI[im] += (dnr_cn.dnHeI[im]) * dtnew;
                    Er.EHeI[im] += (dEr_cn.dEHeI[im]) * dtnew;
                    nr.xnHeI[im] += (dnr_cn.dxnHeI[im]) * dtnew;
                    Er.xEHeI[im] += (dEr_cn.dxEHeI[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nHeI[im] * accur) / abs(dnr_cn.dnHeI[im] * dtnew) * dtnew;
                    dth2 = (Er.EHeI[im] * accur) / abs(dEr_cn.dEHeI[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nHeI[im] += (dnr_cn.dnHeI[im]) * dth1;
                    Er.EHeI[im] += (dEr_cn.dEHeI[im]) * dth1;
                    nr.xnHeI[im] += (dnr_cn.dxnHeI[im]) * dth1;
                    Er.xEHeI[im] += (dEr_cn.dxEHeI[im]) * dth1;
                }

                if (((abs(dnr_cn.dnHeII[im] * dtnew) / nr.nHeII[im]) <= accur) | ((abs(dEr_cn.dEHeII[im] * dtnew) / Er.EHeII[im]) <= accur)) {
                    nr.nHeII[im] += (dnr_cn.dnHeII[im]) * dtnew;
                    Er.EHeII[im] += (dEr_cn.dEHeII[im]) * dtnew;
                    nr.xnHeII[im] += (dnr_cn.dxnHeII[im]) * dtnew;
                    Er.xEHeII[im] += (dEr_cn.dxEHeII[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nHeII[im] * accur) / abs(dnr_cn.dnHeII[im] * dtnew) * dtnew;
                    dth2 = (Er.EHeII[im] * accur) / abs(dEr_cn.dEHeII[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nHeII[im] += (dnr_cn.dnHeII[im]) * dth1;
                    Er.EHeII[im] += (dEr_cn.dEHeII[im]) * dth1;
                    nr.xnHeII[im] += (dnr_cn.dxnHeII[im]) * dth1;
                    Er.xEHeII[im] += (dEr_cn.dxEHeII[im]) * dth1;
                }

                if (((abs(dnr_cn.dnHeIII[im] * dtnew) / nr.nHeIII[im]) <= accur) | ((abs(dEr_cn.dEHeIII[im] * dtnew) / Er.EHeIII[im]) <= accur)) {
                    nr.nHeIII[im] += (dnr_cn.dnHeIII[im]) * dtnew;
                    Er.EHeIII[im] += (dEr_cn.dEHeIII[im]) * dtnew;
                    nr.xnHeIII[im] += (dnr_cn.dxnHeIII[im]) * dtnew;
                    Er.xEHeIII[im] += (dEr_cn.dxEHeIII[im]) * dtnew;
                }
                else {
                    dth1 = (nr.nHeIII[im] * accur) / abs(dnr_cn.dnHeIII[im] * dtnew) * dtnew;
                    dth2 = (Er.EHeIII[im] * accur) / abs(dEr_cn.dEHeIII[im] * dtnew) * dtnew;
                    dth1 = min(dth1, dth2);
                    nr.nHeIII[im] += (dnr_cn.dnHeIII[im]) * dth1;
                    Er.EHeIII[im] += (dEr_cn.dEHeIII[im]) * dth1;
                    nr.xnHeIII[im] += (dnr_cn.dxnHeIII[im]) * dth1;
                    Er.xEHeIII[im] += (dEr_cn.dxEHeIII[im]) * dth1;
                }

                if (bimpur) {
                    nr.nCI[im] += (dnr_cn.dnCI[im]) * dtnew;
                    if (nr.nCI[im] < 0) {
                        cout << " '<0 ts' 10 " << im << " " << mytid << " " << scientific << setprecision(5) << dtnew << endl;
                    }
                    nr.nCII[im] += (dnr_cn.dnCII[im]) * dtnew;
                    if (nr.nCII[im] < 0) {
                        cout << " '<0 ts' 11 " << im << " " << mytid << " " << scientific << setprecision(5) << dtnew << endl;
                    }
                    nr.nCIII[im] += (dnr_cn.dnCIII[im]) * dtnew;
                    if (nr.nCIII[im] < 0) {
                        cout << " '<0 ts' 12 " << im << " " << mytid << " " << scientific << setprecision(5) << dtnew << endl;
                    }
                    nr.nCIV[im] += (dnr_cn.dnCIV[im]) * dtnew;
                    if (nr.nCIV[im] < 0) {
                        cout << " '<0 ts' 13 " << im << " " << mytid << " " << scientific << setprecision(5) << dtnew << endl;
                    }
                    nr.nCV[im] += (dnr_cn.dnCV[im]) * dtnew;
                    if (nr.nCV[im] < 0) {
                        cout << " '<0 ts' 14 " << im << " " << mytid << " " << scientific << setprecision(5) << dtnew << endl;
                    }

                    nr.xnCI[im] += (dnr_cn.dxnCI[im]) * dtnew;
                    nr.xnCII[im] += (dnr_cn.dxnCII[im]) * dtnew;
                    nr.xnCIII[im] += (dnr_cn.dxnCIII[im]) * dtnew;
                    nr.xnCIV[im] += (dnr_cn.dxnCIV[im]) * dtnew;
                    nr.xnCV[im] += (dnr_cn.dxnCV[im]) * dtnew;
                }

                // density cannot be lower than nevac
                if (nr.nHi[im] < nevac) {
                    Er.EHi[im] = Er.EHi[im] * nevac / nr.nHi[im];
                    nr.xnHi[im] = nr.xnHi[im] * nevac / nr.nHi[im];
                    Er.xEHi[im] = Er.xEHi[im] * nevac / nr.nHi[im];
                    nr.nHi[im] = nevac;
                }
                if (nr.nH2i[im] < nevac) {
                    Er.EH2i[im] = Er.EH2i[im] * nevac / nr.nH2i[im];
                    nr.xnH2i[im] = nr.xnH2i[im] * nevac / nr.nH2i[im];
                    Er.xEH2i[im] = Er.xEH2i[im] * nevac / nr.nH2i[im];
                    nr.nH2i[im] = nevac;
                }
                if (nr.nH3i[im] < nevac) {
                    Er.EH3i[im] = Er.EH3i[im] * nevac / nr.nH3i[im];
                    nr.xnH3i[im] = nr.xnH3i[im] * nevac / nr.nH3i[im];
                    Er.xEH3i[im] = Er.xEH3i[im] * nevac / nr.nH3i[im];
                    nr.nH3i[im] = nevac;
                }
                if (nr.nHeII[im] < nevac) {
                    Er.EHeII[im] = Er.EHeII[im] * nevac / nr.nHeII[im];
                    nr.xnHeII[im] = nr.xnHeII[im] * nevac / nr.nHeII[im];
                    Er.xEHeII[im] = Er.xEHeII[im] * nevac / nr.nHeII[im];
                    nr.nHeII[im] = nevac;
                }
                if (nr.nHeIII[im] < nevac) {
                    Er.EHeIII[im] = Er.EHeIII[im] * nevac / nr.nHeIII[im];
                    nr.xnHeIII[im] = nr.xnHeIII[im] * nevac / nr.nHeIII[im];
                    Er.xEHeIII[im] = Er.xEHeIII[im] * nevac / nr.nHeIII[im];
                    nr.nHeIII[im] = nevac;
                }

                neh = nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + 2.0 * nr.nHeIII[im];
                if (bimpur) {
                    neh += nr.nCII[im] + 2.0 * nr.nCIII[im] + 3.0 * nr.nCIV[im] + 4.0 * nr.nCV[im];
                }
                // HERE
                // cout << neh << endl;
                Er.Ee[im] = Er.Ee[im] * neh / nr.ne[im];
                nr.xne[im] = nr.xne[im] * neh / nr.ne[im];
                Er.xEe[im] = Er.xEe[im] * neh / nr.ne[im];
                nr.ne[im] = neh;
            }

            // if (bkipt) {
//#pragma omp parallel for
            for (int im = 1; im < NMESHP - 1; ++im) {
                nr.xne[im] = (nr.ne[im + 1] - nr.ne[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEe[im] = (Er.Ee[im + 1] - Er.Ee[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnH[im] = (nr.nH[im + 1] - nr.nH[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEH[im] = (Er.EH[im + 1] - Er.EH[im - 1]) / (aR[im + 1] - aR[im - 1]);
                //
                nr.xnH2[im] = (nr.nH2[im + 1] - nr.nH2[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEH2[im] = (Er.EH2[im + 1] - Er.EH2[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnHi[im] = (nr.nHi[im + 1] - nr.nHi[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEHi[im] = (Er.EHi[im + 1] - Er.EHi[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnH2i[im] = (nr.nH2i[im + 1] - nr.nH2i[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEH2i[im] = (Er.EH2i[im + 1] - Er.EH2i[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnH3i[im] = (nr.nH3i[im + 1] - nr.nH3i[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEH3i[im] = (Er.EH3i[im + 1] - Er.EH3i[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnHeI[im] = (nr.nHeI[im + 1] - nr.nHeI[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEHeI[im] = (Er.EHeI[im + 1] - Er.EHeI[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnHeII[im] = (nr.nHeII[im + 1] - nr.nHeII[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEHeII[im] = (Er.EHeII[im + 1] - Er.EHeII[im - 1]) / (aR[im + 1] - aR[im - 1]);

                nr.xnHeIII[im] = (nr.nHeIII[im + 1] - nr.nHeIII[im - 1]) / (aR[im + 1] - aR[im - 1]);
                Er.xEHeIII[im] = (Er.EHeIII[im + 1] - Er.EHeIII[im - 1]) / (aR[im + 1] - aR[im - 1]);
            }

            nr.xne[0] = (nr.ne[1] - nr.ne[0]) / (aR[1] - aR[0]);
            Er.xEe[0] = (Er.Ee[1] - Er.Ee[0]) / (aR[1] - aR[0]);
            nr.xne[NMESHP - 1] = (nr.ne[NMESHP - 1] - nr.ne[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
            Er.xEe[NMESHP - 1] = (Er.Ee[NMESHP - 1] - Er.Ee[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

            nr.xnHi[0] = (nr.nHi[1] - nr.nHi[0]) / (aR[1] - aR[0]);
            Er.xEHi[0] = (Er.EHi[1] - Er.EHi[0]) / (aR[1] - aR[0]);
            nr.xnHi[NMESHP - 1] = (nr.nHi[NMESHP - 1] - nr.nHi[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
            Er.xEHi[NMESHP - 1] = (Er.EHi[NMESHP - 1] - Er.EHi[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

            nr.xnH2i[0] = (nr.nH2i[1] - nr.nH2i[0]) / (aR[1] - aR[0]);
            Er.xEH2i[0] = (Er.EH2i[1] - Er.EH2i[0]) / (aR[1] - aR[0]);
            nr.xnH2i[NMESHP - 1] = (nr.nH2i[NMESHP - 1] - nr.nH2i[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
            Er.xEH2i[NMESHP - 1] = (Er.EH2i[NMESHP - 1] - Er.EH2i[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

            nr.xnH3i[0] = (nr.nH3i[1] - nr.nH3i[0]) / (aR[1] - aR[0]);
            Er.xEH3i[0] = (Er.EH3i[1] - Er.EH3i[0]) / (aR[1] - aR[0]);
            nr.xnH3i[NMESHP - 1] = (nr.nH3i[NMESHP - 1] - nr.nH3i[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
            Er.xEH3i[NMESHP - 1] = (Er.EH3i[NMESHP - 1] - Er.EH3i[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

            nr.xnHeII[0] = (nr.nHeII[1] - nr.nHeII[0]) / (aR[1] - aR[0]);
            Er.xEHeII[0] = (Er.EHeII[1] - Er.EHeII[0]) / (aR[1] - aR[0]);
            nr.xnHeII[NMESHP - 1] = (nr.nHeII[NMESHP - 1] - nr.nHeII[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
            Er.xEHeII[NMESHP - 1] = (Er.EHeII[NMESHP - 1] - Er.EHeII[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

            nr.xnHeIII[0] = (nr.nHeIII[1] - nr.nHeIII[0]) / (aR[1] - aR[0]);
            Er.xEHeIII[0] = (Er.EHeIII[1] - Er.EHeIII[0]) / (aR[1] - aR[0]);
            nr.xnHeIII[NMESHP - 1] = (nr.nHeIII[NMESHP - 1] - nr.nHeIII[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
            Er.xEHeIII[NMESHP - 1] = (Er.EHeIII[NMESHP - 1] - Er.EHeIII[NMESHP - 2]) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

//#pragma omp parallel for
            for (int im = 0; im < NMESHP; ++im) {
                // temperature cannot be lower than vessel temperature (Ta0)
                if (Er.Ee[im] / (3.0 / 2.0 * nr.ne[im]) < 0.1) {
                    Er.xEe[im] = Er.xEe[im] * (3.0 / 2.0 * nr.ne[im]) * Ta0 / Er.Ee[im];
                    Er.Ee[im] = (3.0 / 2.0 * nr.ne[im]) * 0.1;
                }
                if (Er.EH[im] / (3.0 / 2.0 * nr.nH[im]) < 0.1) {
                    Er.xEH[im] = Er.xEH[im] * (3.0 / 2.0 * nr.nH[im]) * Ta0 / Er.EH[im];
                    Er.EH[im] = (3.0 / 2.0 * nr.nH[im]) * 0.1;
                }
                if (Er.EH2[im] / (3.0 / 2.0 * nr.nH2[im]) < Ta0) {
                    Er.xEH2[im] = Er.xEH2[im] * (3.0 / 2.0 * nr.nH2[im]) * Ta0 / Er.EH2[im];
                    Er.EH2[im] = (3.0 / 2.0 * nr.nH2[im]) * Ta0;
                }
                if (Er.EHi[im] / (3.0 / 2.0 * nr.nHi[im]) < 0.1) {
                    Er.xEHi[im] = Er.xEHi[im] * (3.0 / 2.0 * nr.nHi[im]) * Ta0 / Er.EHi[im];
                    Er.EHi[im] = (3.0 / 2.0 * nr.nHi[im]) * 0.1;
                }
                if (Er.EH2i[im] / (3.0 / 2.0 * nr.nH2i[im]) < 0.1) {
                    Er.xEH2i[im] = Er.xEH2i[im] * (3.0 / 2.0 * nr.nH2i[im]) * Ta0 / Er.EH2i[im];
                    Er.EH2i[im] = (3.0 / 2.0 * nr.nH2i[im]) * 0.1;
                }
                if (Er.EH3i[im] / (3.0 / 2.0 * nr.nH3i[im]) < 0.1) {
                    Er.xEH3i[im] = Er.xEH3i[im] * (3.0 / 2.0 * nr.nH3i[im]) * Ta0 / Er.EH3i[im];
                    Er.EH3i[im] = (3.0 / 2.0 * nr.nH3i[im]) * 0.1;
                }
                if (Er.EHeI[im] / (3.0 / 2.0 * nr.nHeI[im]) < Ta0) {
                    Er.xEHeI[im] = Er.xEHeI[im] * (3.0 / 2.0 * nr.nHeI[im]) * Ta0 / Er.EHeI[im];
                    Er.EHeI[im] = (3.0 / 2.0 * nr.nHeI[im]) * Ta0;
                }
                if (Er.EHeII[im] / (3.0 / 2.0 * nr.nHeII[im]) < 0.1) {
                    Er.xEHeII[im] = Er.xEHeII[im] * (3.0 / 2.0 * nr.nHeII[im]) * Ta0 / Er.EHeII[im];
                    Er.EHeII[im] = (3.0 / 2.0 * nr.nHeII[im]) * 0.1;
                }
                if (Er.EHeIII[im] / (3.0 / 2.0 * nr.nHeIII[im]) < 0.1) {
                    Er.xEHeIII[im] = Er.xEHeIII[im] * (3.0 / 2.0 * nr.nHeIII[im]) * Ta0 / Er.EHeIII[im];
                    Er.EHeIII[im] = (3.0 / 2.0 * nr.nHeIII[im]) * 0.1;
                }
                // temperature cannot be higher than ...
                if (Er.Ee[im] / (3.0 / 2.0 * nr.ne[im]) > 1e3) {
                    Er.xEe[im] = Er.xEe[im] * (3.0 / 2.0 * nr.ne[im]) * Ta0 / Er.Ee[im];
                    Er.Ee[im] = (3.0 / 2.0 * nr.ne[im]) * 1e3;
                }
                if (Er.EH[im] / (3.0 / 2.0 * nr.nH[im]) > 1e3) {
                    Er.xEH[im] = Er.xEH[im] * (3.0 / 2.0 * nr.nH[im]) * Ta0 / Er.EH[im];
                    Er.EH[im] = (3.0 / 2.0 * nr.nH[im]) * 1e3;
                }
                if (Er.EHi[im] / (3.0 / 2.0 * nr.nHi[im]) > 1e3) {
                    Er.xEHi[im] = Er.xEHi[im] * (3.0 / 2.0 * nr.nHi[im]) * Ta0 / Er.EHi[im];
                    Er.EHi[im] = (3.0 / 2.0 * nr.nHi[im]) * 1e3;
                }
                if (Er.EH2i[im] / (3.0 / 2.0 * nr.nH2i[im]) > 1e3) {
                    Er.xEH2i[im] = Er.xEH2i[im] * (3.0 / 2.0 * nr.nH2i[im]) * Ta0 / Er.EH2i[im];
                    Er.EH2i[im] = (3.0 / 2.0 * nr.nH2i[im]) * 1e3;
                }
                if (Er.EH3i[im] / (3.0 / 2.0 * nr.nH3i[im]) > 1e3) {
                    Er.xEH3i[im] = Er.xEH3i[im] * (3.0 / 2.0 * nr.nH3i[im]) * Ta0 / Er.EH3i[im];
                    Er.EH3i[im] = (3.0 / 2.0 * nr.nH3i[im]) * 1e3;
                }
                if (Er.EHeII[im] / (3.0 / 2.0 * nr.nHeII[im]) > 1e3) {
                    Er.xEHeII[im] = Er.xEHeII[im] * (3.0 / 2.0 * nr.nHeII[im]) * Ta0 / Er.EHeII[im];
                    Er.EHeII[im] = (3.0 / 2.0 * nr.nHeII[im]) * 1e3;
                }
                if (Er.EHeIII[im] / (3.0 / 2.0 * nr.nHeIII[im]) > 1e3) {
                    Er.xEHeIII[im] = Er.xEHeIII[im] * (3.0 / 2.0 * nr.nHeIII[im]) * Ta0 / Er.EHeIII[im];
                    Er.EHeIII[im] = (3.0 / 2.0 * nr.nHeIII[im]) * 1e3;
                }
            }

            dtout += dtnew / Nloopsave; // calculate average time step for output to screen

            tmain += dtnew;

            AvTimeStep = (Nit * AvTimeStep + (dtnew)) / (Nit + 1);
            oldtstep = dtnew;
            olddaccur = daccur;

            nr_p1 = nr;
            Er_p1 = Er;
        }

        dtcalc = dt; //
        if (dtcalc < dtmin * tfact) {
            dtcalc = dtmin * tfact;
        }
        if (dtcalc > dtmax) {
            dtcalc = dtmax;
        } // (in case time step is larger than max allowed step, does happen...)

        if (tmain + dtcalc > tmainend) {
            dtcalc = tmainend - tmain;
        }

        dtnew = dtcalc;
        ++TimeStepCounter;

        if (isnan(nr.ne[0])) {
            tmain = tmainend;
        }
    }
}

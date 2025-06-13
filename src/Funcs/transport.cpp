#include "transport.h"

void transpCoef() { // called in Tomator1D.cpp
    #ifdef debug
    	double Dmin = 10000000.0;
        double Dmax = -1.0;
	double mfpiV[NMESHP] = {0.0};
	double griV[NMESHP] = {0.0};
	double nuiV[NMESHP] = {0.0};
    #endif
    //double Z = 1.0;  // charge state
    if (btunedv) { // calculate prefactor Dfsave and Vfsave
        nefact = (nr.ne[ic]) / (necfix);

        if (btunevleft)                                                              // 1 for Bv scan
        {                                                                            // high density
            PIerrorV += 1.0 * (nr.ne[il] / nefact - nelfix) / nelfix * dtnew / tauV; // 0.0005
            PIerrorD += 1.0 * (nerfix - nr.ne[ir] / nefact) / nerfix * dtnew / tauD; // 0.005
        } else                                                                       // 0 for power scan
        {                                                                            // low density
            PIerrorV -= 1.0 * (nr.ne[ir] / nefact - nerfix) / nerfix * dtnew / tauV; // 0.0005
            PIerrorD += 1.0 * (nelfix - nr.ne[il] / nefact) / nelfix * dtnew / tauD; // 0.005
        }

        if (PIerrorV < -4.25) /////////////////
        {
            PIerrorV = -4.25;
        }
        if (PIerrorV > 20.0) {
            PIerrorV = 20.0;
        }
        if (PIerrorD < -0.99) /////////////////
        {
            PIerrorD = -0.99;
        }
        if (PIerrorD > 10.0) // 1.0*exp(-tmain/0.003)+0.2)
        {
            PIerrorD = 10.0; // *exp(-tmain/0.003)+0.2;
        }

        if (btunevleft) // Bv scan
        {
            Dfsave = 1.0 + max(-0.99, (0.2 * (nerfix - nr.ne[ir] / nefact) / nerfix + PIerrorD));
            Vfsave = (4.5 + PIerrorV);
            Dfsave = 0.3 * Dfsave;
            Vfsave = 0.3 * Vfsave;
        } else // power scan
        {
            Dfsave = max(0.1, (1.0 + PIerrorD) * (1.0 + 0.2 * (nelfix - nr.ne[il] / nefact) / nelfix));
            Vfsave = max(0.3, (4.5 + PIerrorV) * (1.0 - 0.5 * (nr.ne[ir] / nefact - nerfix) / nerfix));
            Dfsave = 0.3 * Dfsave;
            Vfsave = 0.3 * Vfsave;
        }
    }
    else // case of bDscaling or dDbohm
    { 
        PIerrorV = 0.0; 
        PIerrorD = 0.0; 
        Dfsave = Dfact;
        Vfsave = Vfact;
    }

    if (bVscaling) { // also used for btunedv
        for (int im = 0; im < NMESHP; ++im) {
            Vionh[im] = max(nr.ne[im] * Tr.Te[im] + nr.nHi[im] * Tr.THi[im] + nr.nH2i[im] * Tr.TH2i[im] + nr.nH3i[im] * Tr.TH3i[im] + nr.nHeII[im] * Tr.THeII[im] + nr.nHeIII[im] * Tr.THeIII[im], nr.ne[im] * 0.05) //(Ta0*2.0)
                        / nr.ne[im];
        }
    }

    if (bDscaling) { // also used for btunedv
        for (int im = 0; im < NMESHP; ++im) {
            mfpi = nr.nHi[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.THi[im] * nr.nHi[im] / nr.ne[im]) / 1.0) / max(colrate.nuHi[im], 1e4);
            gri = nr.nHi[im] * 1.02e2 / 1.0 * sqrt(1.0 * Tr.THi[im]) / Br[im] / 1e4; // cm
            nui = nr.nHi[im] * max(colrate.nuHi[im], 1e4);
            // ai = nr.nHi[im]*  8.254e3/1.0/(aR[im]/100.0)*(Tr.Te[im]+Tr.THi[im]*nr.nHi[im]/nr.ne[im])*11600.0 ;
            // if (im==90 & mytid==2) {cout << " 1.  " << mfpi/nr.nHi[im] << "   " << gri/nr.nHi[im] << "   " << nui << endl;}
            mfpi += nr.nH2i[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.TH2i[im] * nr.nH2i[im] / nr.ne[im]) / 2.0) / max(colrate.nuH2i[im], 1e4);
            gri += nr.nH2i[im] * 1.02e2 / 1.0 * sqrt(2.0 * Tr.TH2i[im]) / Br[im] / 1e4; // cm
            nui += nr.nH2i[im] * max(colrate.nuH2i[im], 1e4);
            // ai += nr.nH2i[im]*  8.254e3/2.0/(aR[im]/100.0)*(Tr.Te[im]+Tr.TH2i[im]*nr.nH2i[im]/nr.ne[im])*11600.0 ;
            // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
            mfpi += nr.nH3i[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.TH3i[im] * nr.nH3i[im] / nr.ne[im]) / 3.0) / max(colrate.nuH3i[im], 1e4);
            gri += nr.nH3i[im] * 1.02e2 / 1.0 * sqrt(3.0 * Tr.TH3i[im]) / Br[im] / 1e4; // cm
            nui += nr.nH3i[im] * max(colrate.nuH3i[im], 1e4);
            // ai += nr.nH3i[im]*  8.254e3/3.0/(aR[im]/100.0)*(Tr.Te[im]+Tr.TH3i[im]*nr.nH3i[im]/nr.ne[im])*11600.0 ;
            // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
            mfpi += nr.nHeII[im] * 9.79e5 * sqrt((Tr.Te[im] + Tr.THeII[im] * nr.nHeII[im] / nr.ne[im]) / 4.0) / max(colrate.nuHeII[im], 1e4);
            gri += nr.nHeII[im] * 1.02e2 / 1.0 * sqrt(4.0 * Tr.THeII[im]) / Br[im] / 1e4; // cm
            nui += nr.nHeII[im] * max(colrate.nuHeII[im], 1e4);
            // ai += nr.nHeII[im]*  8.254e3/4.0/(aR[im]/100.0)*(Tr.Te[im]+Tr.THeII[im]*nr.nHeII[im]/nr.ne[im])*11600.0 ;
            // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
            mfpi += nr.nHeIII[im] * 9.79e5 * sqrt(2.0 * (Tr.Te[im] + Tr.THeIII[im] * nr.nHeIII[im] / nr.ne[im]) / 4.0) / max(colrate.nuHeIII[im], 1e4);
            gri += nr.nHeIII[im] * 1.02e2 / 2.0 * sqrt(4.0 * Tr.THeIII[im]) / Br[im] / 1e4; // cm
            nui += nr.nHeIII[im] * max(colrate.nuHeIII[im], 1e4);
            // ai += nr.nHeIII[im]*  8.254e3/4.0/(aR[im]/100.0)*(Tr.Te[im]+Tr.THeIII[im]*nr.nHeIII[im]/nr.ne[im])*11600.0 ;
            // if (im==90 & mytid==2) {cout << "     " << mfpi << "   " << gri << "   " << nu << endl;}
            mfpi = mfpi / (nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + nr.nHeIII[im]);
            gri = gri / (nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + nr.nHeIII[im]);
            nui = nui / (nr.nHi[im] + nr.nH2i[im] + nr.nH3i[im] + nr.nHeII[im] + nr.nHeIII[im]);
            // ai=ai/(nr.nHi[im]+nr.nH2i[im]+nr.nH3i[im]+nr.nHeII[im]+nr.nHeIII[im]);

            Dionh[im] = 0.333 * nui * mfpi * (gri + mfpi * Bh / Bt);

    	    #ifdef debug
		//cout << "MFPi = " << mfpi/100 << "[m]" << endl;
		mfpiV[im] = mfpi/100;
		//cout << "GRi = " << gri/100 << "[m]" << endl;
		griV[im] = gri/100;
		//cout << "NUi = " << nui << "[Hz]" << endl;
		nuiV[im] = nui;
		cout << "Dionh = " << Dionh[im]/10000 << "[m2/s]" << endl;
		cout << "Br = " << Br[im] << "[T]" << endl;
		cout << "Bh = " << Bh << "[T]" << endl;
		cout << "Bt = " << Bt << "[T]" << endl;
	    #endif
        }
    }

    // Dion rescaling
    if (bDfix) {
        for (int im = 0; im < NMESHP; ++im) {
            Dion[im] = Dfix;
            #ifdef debug
            if(Dion[im]>Dmax){
                Dmax = Dion[im];
            }
            if(Dion[im]<Dmin){
                Dmin = Dion[im];
            }
            #endif
        }
    } // END bDfix
    else if (bDscaling) { // also used for btunedv
        for (int im = 0; im < NMESHP; ++im) {
            Dion[im] = max(1e2, Dfsave * Dionh[im]);
            #ifdef debug
            if(Dion[im]>Dmax){
                Dmax = Dion[im];
            }
            if(Dion[im]<Dmin){
                Dmin = Dion[im];
            }
            #endif
        }
    } // END bDscaling
    else if (bDbohm) {
        // #pragma omp parallel for private(mfpi, gri, nui)
        for (int im = 0; im < NMESHP; ++im) {
            Dionh[im] = max(nr.ne[im] * Tr.Te[im] + nr.nHi[im] * Tr.THi[im] + nr.nH2i[im] * Tr.TH2i[im] + nr.nH3i[im] * Tr.TH3i[im] + nr.nHeII[im] * Tr.THeII[im] + nr.nHeIII[im] * Tr.THeIII[im], nr.ne[im] * 0.05) //(Ta0*2.0)
                        / nr.ne[im];
            #ifdef debug
            if(Dion[im]>Dmax){
                Dmax = Dion[im];
            }
            if(Dion[im]<Dmin){
                Dmin = Dion[im];
            }
            #endif
        }
        for (int im = 0; im < NMESHP; ++im) {
            Dion[im] = max(1e2, Dfsave * Dionh[im] / Br[im]);
            #ifdef debug
            if(Dion[im]>Dmax){
                Dmax = Dion[im];
            }
            if(Dion[im]<Dmin){
                Dmin = Dion[im];
            }
            #endif
        }
    } // END bDbohm

    // Vion rescaling
    if (bVfix) {
        for (int im = 0; im < NMESHP; ++im) {
            Vionh[im] = Vfix;
        }
        // cout << Vfix << endl;
    } // END bVfix
    else if (bVscaling) { // also used for btunedv
        // Vfsave = (1.0-1.0*exp(-tmain/5e-3))*Vfact*2.0/R;
        for (int im = 0; im < NMESHP; ++im) {
            if (veq == 8) { // Eq (8)
                Vionh[im] = Vfsave * 100.0 * Vionh[im] / Br[im];
            } else if (veq == 9) { // Eq (9)
                Vfsave = Vfsave * 2.0 / R;
                Vionh[im] = Vfsave * Dion[im];
            }
        }
    } // END bVscaling
    else { // Te dependency only --> there is no bool switch for this (yet)
        for (int im = 0; im < NMESHP; ++im) {
            Vionh[im] = Vfsave * 100.0 * Vionh[im];
        }
    }
    #ifdef debug
	cout << "Dfsave = " << Dfsave << "[-]" << endl;
	cout << "Min. diff. coeff. radial = " << Dmin/10000 << "[m2/s]; max. = " << Dmax/10000 << "[m2/s]" << endl;
    #endif

    if (bVscaling) { // also used for btunedv
        // weighted average over radius for Vion
        double nemax = nr.ne[0];
        double Temax = Tr.Te[0];
        for (int id = 1; id < NMESHP - 1; ++id) { // line integrated density
            nemax = max(nemax, nr.ne[id]);
            Temax = max(Temax, Tr.Te[id]);
        }
        Vionh[0] = 0.0;                           // Vionh[0]*nr.ne[0]*Tr.Te[0]*aR[0]*(aR[1]-aR[0])+Vionh[1]*nr.ne[1]*Tr.Te[1]*aR[1]*(aR[1]-aR[0]);//+Vionh[2]*nr.ne[2]*aR[2]*(aR[2]-aR[1]);
        Vionh[1] = 0.0;                           // nr.ne[0]*Tr.Te[0]*aR[0]*(aR[1]-aR[0]) + Tr.Te[1]*Tr.Te[1]*aR[1]*(aR[1]-aR[0]);// + nr.ne[2]*aR[2]*(aR[2]-aR[1]);
        for (int id = 2; id < NMESHP - 1; ++id) { // line integrated density
            if (nr.ne[id] > nemax / 5.0) {
                Vionh[0] += Vionh[id] * nr.ne[id] * Tr.Te[id] * aR[id];
                Vionh[1] += nr.ne[id] * Tr.Te[id] * aR[id]; // inteprof += 2.0*b*pi*0.5*(eprof[id]+eprof[id+1])*(pow(aR[id+1],2.0)-pow(aR[id],2.0)); // 20181016
            }
        }
        Vionh[0] = (Vionh[0] / Vionh[1]);
        for (int id = 1; id < NMESHP - 1; ++id) { // line integrated density
            Vionh[id] = Vionh[0];
        }
    }
    // no convection if density is below the vacuum density
    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        // Vion[im] = Vionh[im] * pow(1.0 + (nevac / nr.ne[im]), -0.66) * pow(1.0 + (nevac / bnion(im * 2)), -0.66); // switch of convection if density ne or nion is low
        Vion[im] = Vionh[im] * pow(1.0 + (nevac / nr.ne[im]), -0.66); // removed bnion as it is very low since transpcoeffs is now called only once
    }

    // consider to remove this two lines below as they do not seem to do much... 
    //Vion[0] = 0.0;
    //Vion[1] = Vion[1] / 2.0;

    // cout << "Dion[cc] = "  << Dion[cc] << ",  Vion[cc] = " << Vion[cc] << endl;
}

void transportions(double tstep) { // called in solver.cpp within timeStep.cpp in Tomator1D.cpp
    
#pragma omp parallel for collapse(2)
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Vs.coeffRef(2 * i, 2 * j) = VLfbva[i] * Vion[i - 1] + VLfbvb[i] * Vion[i] + VRfbvc[i] * Vion[i + 1] + VRfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVLfbva[i] * Vion[i - 1] + dVLfbvb[i] * Vion[i] + dVRfbvc[i] * Vion[i + 1] + dVRfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VLdfbva[i] * Vion[i - 1] + VLdfbvb[i] * Vion[i] + VRdfbvc[i] * Vion[i + 1] + VRdfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVLdfbva[i] * Vion[i - 1] + dVLdfbvb[i] * Vion[i] + dVRdfbvc[i] * Vion[i + 1] + dVRdfbvb[i] * Vion[i];
            } else if (j == i - 1) {
                Vs.coeffRef(2 * i, 2 * j) = VLfava[i] * Vion[i - 1] + VLfavb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVLfava[i] * Vion[i - 1] + dVLfavb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VLdfava[i] * Vion[i - 1] + VLdfavb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVLdfava[i] * Vion[i - 1] + dVLdfavb[i] * Vion[i];
            } else if (j == i + 1) {
                Vs.coeffRef(2 * i, 2 * j) = VRfcvc[i] * Vion[i + 1] + VRfcvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVRfcvc[i] * Vion[i + 1] + dVRfcvb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VRdfcvc[i] * Vion[i + 1] + VRdfcvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVRdfcvc[i] * Vion[i + 1] + dVRdfcvb[i] * Vion[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Vs.coeffRef(2 * i, 2 * j) = VRfbvc[i] * Vion[i + 1] + VRfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVRfbvc[i] * Vion[i + 1] + dVRfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VRdfbvc[i] * Vion[i + 1] + VRdfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVRdfbvc[i] * Vion[i + 1] + dVRdfbvb[i] * Vion[i];
            } else if (j == i + 1) {
                Vs.coeffRef(2 * i, 2 * j) = VRfcvc[i] * Vion[i + 1] + VRfcvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVRfcvc[i] * Vion[i + 1] + dVRfcvb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VRdfcvc[i] * Vion[i + 1] + VRdfcvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVRdfcvc[i] * Vion[i + 1] + dVRdfcvb[i] * Vion[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Vs.coeffRef(2 * i, 2 * j) = VLfbva[i] * Vion[i - 1] + VLfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVLfbva[i] * Vion[i - 1] + dVLfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VLdfbva[i] * Vion[i - 1] + VLdfbvb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVLdfbva[i] * Vion[i - 1] + dVLdfbvb[i] * Vion[i];
            } else if (j == i - 1) {
                Vs.coeffRef(2 * i, 2 * j) = VLfava[i] * Vion[i - 1] + VLfavb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j) = dVLfava[i] * Vion[i - 1] + dVLfavb[i] * Vion[i];
                Vs.coeffRef(2 * i, 2 * j + 1) = VLdfava[i] * Vion[i - 1] + VLdfavb[i] * Vion[i];
                Vs.coeffRef(2 * i + 1, 2 * j + 1) = dVLdfava[i] * Vion[i - 1] + dVLdfavb[i] * Vion[i];
            }
        }
    }
    // BCs ions
    if (1) {
	if(fixBCs) {
	lambda = lam_dec_length_ions;
	bnion(0 + 1) = bnion(0) / lambda;
	lambda = lam_dec_length_energy;
        bEion(0 + 1) = bEion(0) / lambda;
	bEelec(0 + 1) = bnion(0) * Z * 3.0 / 2.0 * Tr.Te[0] /lambda;
	}
	else {
        lambda = sqrt(Dion[0] * (0.66 * 2.0 * pi * aR[0] / nlimiters) /
                      (9.79e5 * sqrt(Z / mu * (Tr.Te[0] + bnion(0) / nr.ne[0] * bEion(0) / (3.0 / 2.0 * bnion(0)))) * pow(1.0 + (nevac / nr.ne[0]), -0.66) * pow(1.0 + (nevac / bnion(0)), -0.66)));
        if (lambda < 0.0) {
            cout << "lambda<0.0 HFS after" << endl;
	}
        bnion(0 + 1) = bnion(0) / lambda;
        bEion(0 + 1) = bEion(0) / (lambda / sqrt(gEe) * sqrt(gEd));
        // bEion(0+1)            = bnion(0)/lambda*3.0/2.0*Tion[0]   ;//                                                          * pow(1.0+(nevac/nr.ne[0]),-0.66) * pow(1.0+(nevac/bEion(0)),-0.66) * pow(1.0+Ta0/Tion[0],-0.66) ;
        bEelec(0 + 1) = bnion(0) / (lambda / sqrt(gEe) * sqrt(gEd)) * Z * 3.0 / 2.0 * Tr.Te[0];
	}

	if(fixBCs) {
	lambda = lam_dec_length_ions;
	bnion((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) / lambda;
	lambda = lam_dec_length_energy;
        bEion((NMESHP - 1) * 2 + 1) = -bEion((NMESHP - 1) * 2) / lambda;
	bEelec((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) * Z * 3.0 / 2.0 * Tr.Te[NMESHP - 1] / lambda;
	}
	else {
        lambda = sqrt(Dion[NMESHP - 1] * (0.66 * 2.0 * pi * aR[NMESHP - 1] / nlimiters) /
                      (9.79e5 * sqrt(Z / mu * (Tr.Te[NMESHP - 1] + bnion((NMESHP - 1) * 2) / nr.ne[NMESHP - 1] * (bEion((NMESHP - 1) * 2) / (3.0 / 2.0 * bnion((NMESHP - 1) * 2))))) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66) * pow(1.0 + (nevac / bnion((NMESHP - 1) * 2)), -0.66)));
        if (lambda < 0.0) {
            cout << "lambda<0.0 LFS after" << endl;
        }
        bnion((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) / lambda;
        bEion((NMESHP - 1) * 2 + 1) = -bEion((NMESHP - 1) * 2) / (lambda / sqrt(gEe) * sqrt(gEd));
        // bEion((NMESHP-1)*2+1) = -bnion((NMESHP-1)*2)/lambda*3.0/2.0*Tion[NMESHP-1]    ;//                                             * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) * pow(1.0+(nevac/bEion((NMESHP-1)*2)),-0.66) * pow(1.0+Ta0/Tion[NMESHP-1],-0.66);
        bEelec((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) / (lambda / sqrt(gEe) * sqrt(gEd)) * Z * 3.0 / 2.0 * Tr.Te[NMESHP - 1];
	}

	if (fixBCs) {
	lambda = lam_dec_length_ions;
	bnion1(0 + 1) = bnion1(0) / lambda;
	lambda = lam_dec_length_energy;
        bEion1(0 + 1) = bEion1(0) / lambda;
	bEelec1(0 + 1) = bnion1(0) / lambda * Z * 3.0 / 2.0 * Tr.Te[0];
	}
	else {
        lambda = sqrt(Dion[0] * (0.66 * 2.0 * pi * aR[0] / nlimiters) /
                      (9.79e5 * sqrt(Z / mu * (Tr.Te[0] + bnion1(0) / nr.ne[0] * bEion1(0) / (3.0 / 2.0 * bnion1(0)))) * pow(1.0 + (nevac / nr.ne[0]), -0.66) * pow(1.0 + (nevac / bnion1(0)), -0.66)));
        if (lambda < 0.0) {
            cout << "lambda<0.0 HFS after" << endl;
        }
        bnion1(0 + 1) = bnion1(0) / lambda;
        bEion1(0 + 1) = bEion1(0) / (lambda / sqrt(gEe) * sqrt(gEd));
        // bEion1(0+1)            = bnion1(0)/lambda*3.0/2.0*Tion[0]   ;//                                                          * pow(1.0+(nevac/nr.ne[0]),-0.66) * pow(1.0+(nevac/bEion1(0)),-0.66) * pow(1.0+Ta0/Tion[0],-0.66) ;
        bEelec1(0 + 1) = bnion1(0) / (lambda / sqrt(gEe) * sqrt(gEd)) * Z * 3.0 / 2.0 * Tr.Te[0];
	}

	if(fixBCs) {
	lambda = lam_dec_length_ions;
	bnion1((NMESHP - 1) * 2 + 1) = -bnion1((NMESHP - 1) * 2) / lambda;
	lambda = lam_dec_length_energy;
        bEion1((NMESHP - 1) * 2 + 1) = -bEion1((NMESHP - 1) * 2) / lambda;
	bEelec1((NMESHP - 1) * 2 + 1) = -bnion1((NMESHP - 1) * 2) / lambda * Z * 3.0 / 2.0 * Tr.Te[NMESHP - 1];
	}
	else {
        lambda = sqrt(Dion[NMESHP - 1] * (0.66 * 2.0 * pi * aR[NMESHP - 1] / nlimiters) /
                      (9.79e5 * sqrt(Z / mu * (Tr.Te[NMESHP - 1] + bnion1((NMESHP - 1) * 2) / nr.ne[NMESHP - 1] * (bEion1((NMESHP - 1) * 2) / (3.0 / 2.0 * bnion1((NMESHP - 1) * 2))))) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66) * pow(1.0 + (nevac / bnion1((NMESHP - 1) * 2)), -0.66)));
        if (lambda < 0.0) {
            cout << "lambda<0.0 LFS after" << endl;
        }
        bnion1((NMESHP - 1) * 2 + 1) = -bnion1((NMESHP - 1) * 2) / lambda;
        bEion1((NMESHP - 1) * 2 + 1) = -bEion1((NMESHP - 1) * 2) / (lambda / sqrt(gEe) * sqrt(gEd));
        // bEion1((NMESHP-1)*2+1) = -bnion1((NMESHP-1)*2)/lambda*3.0/2.0*Tion[NMESHP-1]    ;//                                             * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) * pow(1.0+(nevac/bEion1((NMESHP-1)*2)),-0.66) * pow(1.0+Ta0/Tion[NMESHP-1],-0.66);
        bEelec1((NMESHP - 1) * 2 + 1) = -bnion1((NMESHP - 1) * 2) / (lambda / sqrt(gEe) * sqrt(gEd)) * Z * 3.0 / 2.0 * Tr.Te[NMESHP - 1];
	}
    }

    if (1) {
        if (bnion(0) <= nevac) {
            bnion(0 + 1) = 0.0;
        }
        if (bnion1(0) <= nevac) {
            bnion1(0 + 1) = 0.0;
        }
        if (bnion((NMESHP - 1) * 2) <= nevac) {
            bnion((NMESHP - 1) * 2 + 1) = 0.0;
        }
        if (bnion1((NMESHP - 1) * 2) <= nevac) {
            bnion1((NMESHP - 1) * 2 + 1) = 0.0;
        }

        if (bnion(0 + 1) < 0.0) {
            bnion(0 + 1) = 0.0;
        }
        if (bnion1(0 + 1) < 0.0) {
            bnion1(0 + 1) = 0.0;
        }
        if (bnion((NMESHP - 1) * 2 + 1) > 0.0) {
            bnion((NMESHP - 1) * 2 + 1) = 0.0;
        }
        if (bnion1((NMESHP - 1) * 2 + 1) > 0.0) {
            bnion1((NMESHP - 1) * 2 + 1) = 0.0;
        }

        if (bEion(0) <= 3.0 / 2.0 * bnion(0) * 0.1) {
            bEion(0 + 1) = 0.0;
        }
        if (bEion1(0) <= 3.0 / 2.0 * bnion1(0) * 0.1) {
            bEion1(0 + 1) = 0.0;
        }
        if (bEion((NMESHP - 1) * 2) <= 3.0 / 2.0 * bnion((NMESHP - 1) * 2) * 0.1) {
            bEion((NMESHP - 1) * 2 + 1) = 0.0;
        }
        if (bEion1((NMESHP - 1) * 2) <= 3.0 / 2.0 * bnion1((NMESHP - 1) * 2) * 0.1) {
            bEion1((NMESHP - 1) * 2 + 1) = 0.0;
        }

        if (bEion(0 + 1) < 0.0) {
            bEion(0 + 1) = 0.0;
        }
        if (bEion1(0 + 1) < 0.0) {
            bEion1(0 + 1) = 0.0;
        }
        if (bEion((NMESHP - 1) * 2 + 1) > 0.0) {
            bEion((NMESHP - 1) * 2 + 1) = 0.0;
        }
        if (bnion1((NMESHP - 1) * 2 + 1) > 0.0) {
            bnion1((NMESHP - 1) * 2 + 1) = 0.0;
        }
    }
    // check if bnion has inf
    for (int im = 0; im < NMESHP; ++im) {
        if (isinf(bnion(im * 2)) || isinf(bnion(im * 2 + 1))) {
            // Error: print red message
            cout << "\033[1;31m"
                 << "[Error]: bnion has inf at im = " << im << endl;
            cout << "\033[0m";
            exit(EXIT_FAILURE);
        }
    }
    if (btranspions) {
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                A = Ts - coef1 * Ds * tstep + coef1 * Vs * tstep;
                bnion = (Ts)*bnion * coef2 - Ts * bnion1 * coef3 + coef1 * Ts * Snion * tstep; // - Es*Ssnion*tstep;

                solver.compute(A);
                x = solver.solveWithGuess(bnion, bnion);
                bnion = x;
            }
            #pragma omp section
            {
                B = Ts - coef1 * gEd * Ds * tstep + coef1 * gEv * Vs * tstep;
                bEion = (Ts)*bEion * coef2 - Ts * bEion1 * coef3 + coef1 * Ts * SEion * tstep; // - gEv*Es*SsEion*tstep; // piecewise linear source
                solver2.compute(B);
                x2 = solver2.solveWithGuess(bEion, bEion);
                bEion = x2;
            }
            #pragma omp section
            {
                C = Ts - coef1 * gEd * Ds * tstep + coef1 * gEv * Vs * tstep;
                bEelec = (Ts)*bEelec * coef2 - Ts * bEelec1 * coef3 + coef1 * Ss * SEelec * tstep; // - gEv*Es*SsEelec*tstep ;
                solver3.compute(C);
                x3 = solver3.solveWithGuess(bEelec, bEelec);
                bEelec = x3;
            }
        }
    }

    // #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {
        if (bEion(im * 2) / (3.0 / 2.0 * bnion(im * 2)) < Ta0) {
            bEion(im * 2) = 3.0 / 2.0 * bnion(im * 2) * Ta0;
        }
    }

    if(fixBCs) {
     lambda = lam_dec_length_ions;
     bnion(0 + 1) = bnion(0) / lambda;
     lambda = lam_dec_length_energy;                       
     bEion(0 + 1) = bEion(0) / lambda; 
     bEelec(0 + 1) = bnion(0) / lambda * Z * 3.0 / 2.0 * Tr.Te[0];
     }
    else {
    lambda = sqrt(Dion[0] * (0.66 * 2.0 * pi * aR[0] / nlimiters) /
                  (9.79e5 * sqrt(Z / mu * (Tr.Te[0] + bnion(0) / nr.ne[0] * bEion(0) / (3.0 / 2.0 * bnion(0)))) * pow(1.0 + (nevac / nr.ne[0]), -0.66) * pow(1.0 + (nevac / bnion(0)), -0.66)));
    if (lambda < 0.0) {
        cout << "lambda<0.0 HFS after" << endl;
    }
    bnion(0 + 1) = bnion(0) / lambda;                           //                                                          * pow(1.0+(nevac/nr.ne[0]),-0.66) * pow(1.0+(nevac/bnion(0)),-0.66) ;
    bEion(0 + 1) = bEion(0) / (lambda / sqrt(gEe) * sqrt(gEd)); //                                                          * pow(1.0+(nevac/nr.ne[0]),-0.66) * pow(1.0+(nevac/bEion(0)),-0.66) * pow(1.0+Ta0/Tion[0],-0.66) ;
    // bEion(0+1)            = bnion(0)/lambda*3.0/2.0*Tion[0]   ;//                                                          * pow(1.0+(nevac/nr.ne[0]),-0.66) * pow(1.0+(nevac/bEion(0)),-0.66) * pow(1.0+Ta0/Tion[0],-0.66) ;
    bEelec(0 + 1) = bnion(0) / (lambda / sqrt(gEe) * sqrt(gEd)) * Z * 3.0 / 2.0 * Tr.Te[0]; //                                       * pow(1.0+(nevac/nr.ne[0]),-0.66);

    }

    if(fixBCs) {
	lambda = lam_dec_length_ions;
	bnion((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) / lambda;
	lambda = lam_dec_length_energy;
        bEion((NMESHP - 1) * 2 + 1) = -bEion((NMESHP - 1) * 2) / lambda;
	bEelec((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) * Z * 3.0 / 2.0 * Tr.Te[NMESHP - 1] / lambda;
	}
    else {
    lambda = sqrt(Dion[NMESHP - 1] * (0.66 * 2.0 * pi * aR[NMESHP - 1] / nlimiters) /
                  (9.79e5 * sqrt(Z / mu * (Tr.Te[NMESHP - 1] + bnion((NMESHP - 1) * 2) / nr.ne[NMESHP - 1] * (bEion((NMESHP - 1) * 2) / (3.0 / 2.0 * bnion((NMESHP - 1) * 2))))) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66) * pow(1.0 + (nevac / bnion((NMESHP - 1) * 2)), -0.66)));
    if (lambda < 0.0) {
        cout << "lambda<0.0 LFS after" << endl;
    }
    bnion((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) / lambda;                           //                                             * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) * pow(1.0+(nevac/bnion((NMESHP-1)*2)),-0.66);
    bEion((NMESHP - 1) * 2 + 1) = -bEion((NMESHP - 1) * 2) / (lambda / sqrt(gEe) * sqrt(gEd)); //                                             * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) * pow(1.0+(nevac/bEion((NMESHP-1)*2)),-0.66) * pow(1.0+Ta0/Tion[NMESHP-1],-0.66);
    // bEion((NMESHP-1)*2+1) = -bnion((NMESHP-1)*2)/lambda*3.0/2.0*Tion[NMESHP-1]    ;//                                             * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) * pow(1.0+(nevac/bEion((NMESHP-1)*2)),-0.66) * pow(1.0+Ta0/Tion[NMESHP-1],-0.66);
    bEelec((NMESHP - 1) * 2 + 1) = -bnion((NMESHP - 1) * 2) / (lambda / sqrt(gEe) * sqrt(gEd)) * Z * 3.0 / 2.0 * Tr.Te[NMESHP - 1]; //                    * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66);
    }
}

void transpH(double tstep) { // called in solver.cpp within timeStep.cpp in Tomator1D.cpp
    // Moved to precompiler
    #ifdef debug
    	double DminH = 10000000.0;
        double DmaxH = -1.0;
    	double DminH2 = 10000000.0;
        double DmaxH2 = -1.0;
    #endif
    DionHFS = Dion[0];
    DionLFS = Dion[NMESHP - 1];
    VionLFS = Vion[NMESHP - 1];

    // #pragma omp parallel for private(vth, mfp)
    for (int id = 0; id < NMESHP; ++id) {
        // Neutral H diffusion coefficient
        if (bDfix_neutr) {
        	DH[id] = Dfix_neut;
        	DH2[id] = Dfix_neut;
        }
        else {
        	vth = sqrt((kb * TH_array[id] * 11600.0) / ma) * 100.0;
        	mfp = vth / nuH_array[id];
        	DH[id] = 1.0 / 3.0 * vth * 1.0 / (1.0 / mfp + 1.0 / (a / 2.0)); // DH[id] =  vth*a/2.0 ; // cm2/s
        	#ifdef debug
        		if (DH[id]>DmaxH){
				DmaxH = DH[id];
			}
			if(DH[id]<DminH){
				DminH = DH[id];
			}
		#endif
        	// DH[id] =  1.0/3.0 * vth * a/2.0 ; // DH[id] =  vth*a/2.0 ; // cm2/s
        	// cout << "   " << vth <<  "   " << mfp << "   " << a << "   " << DH[id] << endl;
        
        	// Neutral H2 diffusion coefficient
        	vth = sqrt((kb * TH2_array[id] * 11600.0) / (2.0 * ma)) * 100.0;
        	mfp = vth / nuH2_array[id];
        	DH2[id] = 1.0 / 3.0 * vth * 1.0 / (1.0 / mfp + 1.0 / (a / 2.0)); // DH[id] =  vth*a/2.0 ; // cm2/s
        	// DH2[id] =  1.0/3.0 * vth * a/2.0 ; // DH[id] =  vth*a/2.0 ; // cm2/s // cout << "H2   " << mfp << "   " <<  DH2[id] << endl;
        	#ifdef debug
        		if (DH2[id]>DmaxH2){
				DmaxH2 = DH2[id];
			}
			if(DH2[id]<DminH2){
				DminH2 = DH2[id];
			}
		#endif
    	}
    }
    #ifdef debug
    	cout << "DHmin = " << DminH/1e4 << " [m2/s]; DHmax = " << DmaxH/1e4 << endl;
    	cout << "DH2min = " << DminH2/1e4 << " [m2/s]; DH2max = " << DmaxH2/1e4 << endl;
    #endif

 #pragma omp parallel for collapse(2)
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                DHs.coeffRef(2 * i, 2 * j) = DLfbba[i] * DH[i - 1] + DLfbbb[i] * DH[i] + DRfbbc[i] * DH[i + 1] + DRfbbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DH[i - 1] + dDLfbbb[i] * DH[i] + dDRfbbc[i] * DH[i + 1] + dDRfbbb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DH[i - 1] + DLdfbbb[i] * DH[i] + DRdfbbc[i] * DH[i + 1] + DRdfbbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DH[i - 1] + dDLdfbbb[i] * DH[i] + dDRdfbbc[i] * DH[i + 1] + dDRdfbbb[i] * DH[i];
            } else if (j == i - 1) {
                DHs.coeffRef(2 * i, 2 * j) = DLfaba[i] * DH[i - 1] + DLfabb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DH[i - 1] + dDLfabb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DH[i - 1] + DLdfabb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DH[i - 1] + dDLdfabb[i] * DH[i];
            } else if (j == i + 1) {
                DHs.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DH[i + 1] + DRfcbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DH[i + 1] + dDRfcbb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DH[i + 1] + DRdfcbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DH[i + 1] + dDRdfcbb[i] * DH[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                DHs.coeffRef(2 * i, 2 * j) = DRfbbc[i] * DH[i + 1] + DRfbbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDRfbbc[i] * DH[i + 1] + dDRfbbb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DRdfbbc[i] * DH[i + 1] + DRdfbbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfbbc[i] * DH[i + 1] + dDRdfbbb[i] * DH[i];
            } else if (j == i + 1) {
                DHs.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DH[i + 1] + DRfcbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DH[i + 1] + dDRfcbb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DH[i + 1] + DRdfcbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DH[i + 1] + dDRdfcbb[i] * DH[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                DHs.coeffRef(2 * i, 2 * j) = DLfbba[i] * DH[i - 1] + DLfbbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DH[i - 1] + dDLfbbb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DH[i - 1] + DLdfbbb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DH[i - 1] + dDLdfbbb[i] * DH[i];
            } else if (j == i - 1) {
                DHs.coeffRef(2 * i, 2 * j) = DLfaba[i] * DH[i - 1] + DLfabb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DH[i - 1] + dDLfabb[i] * DH[i];
                DHs.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DH[i - 1] + DLdfabb[i] * DH[i];
                DHs.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DH[i - 1] + dDLdfabb[i] * DH[i];
            }
        }
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                DH2s.coeffRef(2 * i, 2 * j) = DLfbba[i] * DH2[i - 1] + DLfbbb[i] * DH2[i] + DRfbbc[i] * DH2[i + 1] + DRfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DH2[i - 1] + dDLfbbb[i] * DH2[i] + dDRfbbc[i] * DH2[i + 1] + dDRfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DH2[i - 1] + DLdfbbb[i] * DH2[i] + DRdfbbc[i] * DH2[i + 1] + DRdfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DH2[i - 1] + dDLdfbbb[i] * DH2[i] + dDRdfbbc[i] * DH2[i + 1] + dDRdfbbb[i] * DH2[i];
            } else if (j == i - 1) {
                DH2s.coeffRef(2 * i, 2 * j) = DLfaba[i] * DH2[i - 1] + DLfabb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DH2[i - 1] + dDLfabb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DH2[i - 1] + DLdfabb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DH2[i - 1] + dDLdfabb[i] * DH2[i];
            } else if (j == i + 1) {
                DH2s.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DH2[i + 1] + DRfcbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DH2[i + 1] + dDRfcbb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DH2[i + 1] + DRdfcbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DH2[i + 1] + dDRdfcbb[i] * DH2[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                DH2s.coeffRef(2 * i, 2 * j) = DRfbbc[i] * DH2[i + 1] + DRfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDRfbbc[i] * DH2[i + 1] + dDRfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DRdfbbc[i] * DH2[i + 1] + DRdfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfbbc[i] * DH2[i + 1] + dDRdfbbb[i] * DH2[i];
            } else if (j == i + 1) {
                DH2s.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DH2[i + 1] + DRfcbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DH2[i + 1] + dDRfcbb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DH2[i + 1] + DRdfcbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DH2[i + 1] + dDRdfcbb[i] * DH2[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                DH2s.coeffRef(2 * i, 2 * j) = DLfbba[i] * DH2[i - 1] + DLfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DH2[i - 1] + dDLfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DH2[i - 1] + DLdfbbb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DH2[i - 1] + dDLdfbbb[i] * DH2[i];
            } else if (j == i - 1) {
                DH2s.coeffRef(2 * i, 2 * j) = DLfaba[i] * DH2[i - 1] + DLfabb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DH2[i - 1] + dDLfabb[i] * DH2[i];
                DH2s.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DH2[i - 1] + DLdfabb[i] * DH2[i];
                DH2s.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DH2[i - 1] + dDLdfabb[i] * DH2[i];
            }
        }
    }
    // BCs neutrals, H
    // H Edge condition, diffusion
    vth = sqrt((kb * TH_array[0] * 11600.0) / ma) * 100.0;
    bnH(0 + 1) = +1.0 / 2.0 * vth / DH[0] * bnH(0) * (1.0 - RH);
    bEH(0 + 1) = +1.0 / 2.0 * gEdn * vth / (DH[0]) * bnH(0) * 3.0 / 2.0 * TH_array[0] * (1.0 - REH); //*(1.0-RH * (1.0 - ((1.0-REH) * (1.0-Ta0/TH_array[0]))));

    vth = sqrt((kb * TH_array[NMESHP - 1] * 11600.0) / ma) * 100.0;
    bnH((NMESHP - 1) * 2 + 1) = -1.0 / 2.0 * vth / DH[NMESHP - 1] * bnH((NMESHP - 1) * 2) * (1.0 - RH);
    bEH((NMESHP - 1) * 2 + 1) = -1.0 / 2.0 * gEdn * vth / (DH[NMESHP - 1]) * bnH((NMESHP - 1) * 2) * 3.0 / 2.0 * TH_array[NMESHP - 1] * (1.0 - REH); //*(1.0-RH * (1.0 - ((1.0-REH) * (1.0-Ta0/TH_array[NMESHP-1]))));

    vth = sqrt((kb * TH_array[0] * 11600.0) / ma) * 100.0;
    bn1H(0 + 1) = +1.0 / 2.0 * vth / DH[0] * bn1H(0) * (1.0 - RH);
    bE1H(0 + 1) = +1.0 / 2.0 * gEdn * vth / (DH[0]) * bn1H(0) * 3.0 / 2.0 * TH_array[0] * (1.0 - REH); //*(1.0-RH * (1.0 - ((1.0-REH) * (1.0-Ta0/TH_array[0]))));

    vth = sqrt((kb * TH_array[NMESHP - 1] * 11600.0) / ma) * 100.0;
    bn1H((NMESHP - 1) * 2 + 1) = -1.0 / 2.0 * vth / DH[NMESHP - 1] * bn1H((NMESHP - 1) * 2) * (1.0 - RH);
    bE1H((NMESHP - 1) * 2 + 1) = -1.0 / 2.0 * gEdn * vth / (DH[NMESHP - 1]) * bn1H((NMESHP - 1) * 2) * 3.0 / 2.0 * TH_array[NMESHP - 1] * (1.0 - REH); //*(1.0-RH * (1.0 - ((1.0-REH) * (1.0-Ta0/TH_array[NMESHP-1]))));
    // H2 edge conditions
    // if (pcst) {
    //     bEH2(0)              = 3.0/2.0*nH20*Ta0;
    //     bEH2((NMESHP-1)*2)   = 3.0/2.0*nH20*Ta0;
    //     bEH2(0+1)            = (bEH2(2)-bEH2(0))/(aR[1]-aR[0]);
    //     bEH2((NMESHP-1)*2+1) = (bEH2((NMESHP-1)*2)-bEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     bnH2(0)              = nH20;
    //     bnH2((NMESHP-1)*2)   = nH20;
    //     bnH2(0+1)            = (bnH2(2)-bnH2(0))/(aR[1]-aR[0]);
    //     bnH2((NMESHP-1)*2+1) = (bnH2((NMESHP-1)*2)-bnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     SnH2(0)              = 0.0;
    //     SnH2((NMESHP-1)*2)   = 0.0;
    //     SnH2(0+1)            = (SnH2(2)-SnH2(0))/(aR[1]-aR[0]);
    //     SnH2((NMESHP-1)*2+1) = (SnH2((NMESHP-1)*2)-SnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     SEH2(0)              = 0.0;
    //     SEH2((NMESHP-1)*2)   = 0.0;
    //     SEH2(0+1)            = (SEH2(2)-SEH2(0))/(aR[1]-aR[0]);
    //     SEH2((NMESHP-1)*2+1) = (SEH2((NMESHP-1)*2)-SEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    // }
    // else if (quick)
    // {
    //     bEH2(0)              = 3.0/2.0*nH20*TH2_array[0];
    //     bEH2((NMESHP-1)*2)   = 3.0/2.0*nH20*TH2_array[NMESHP-1];
    //     bEH2(0+1)            = (bEH2(2)-bEH2(0))/(aR[1]-aR[0]);
    //     bEH2((NMESHP-1)*2+1) = (bEH2((NMESHP-1)*2)-bEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     bnH2(0)              = nH20;
    //     bnH2((NMESHP-1)*2)   = nH20;
    //     bnH2(0+1)            = (bnH2(2)-bnH2(0))/(aR[1]-aR[0]);
    //     bnH2((NMESHP-1)*2+1) = (bnH2((NMESHP-1)*2)-bnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     SnH2(0)              = 0.0;
    //     SnH2((NMESHP-1)*2)   = 0.0;
    //     SnH2(0+1)            = (SnH2(2)-SnH2(0))/(aR[1]-aR[0]);
    //     SnH2((NMESHP-1)*2+1) = (SnH2((NMESHP-1)*2)-SnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     SEH2(0)              = 0.0;
    //     SEH2((NMESHP-1)*2)   = 0.0;
    //     SEH2(0+1)            = (SEH2(2)-SEH2(0))/(aR[1]-aR[0]);
    //     SEH2((NMESHP-1)*2+1) = (SEH2((NMESHP-1)*2)-SEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    // }
    // else {
    //     if (ncst)
    //     {
    // bEH2(0)              = bEH2(0)/bnH2(0)*nH20; // PSI
    // bEH2((NMESHP-1)*2)   = bEH2((NMESHP-1)*2)/bnH2((NMESHP-1)*2)*nH20; // PSI
    if (nH20 > bnH2(0)) {
        bEH2(0) += 3.0 / 2.0 * Ta0 * (nH20 - bnH2(0));
    } // 20220531
    else {
        bEH2(0) = bEH2(0) / bnH2(0) * nH20;
    }
    if (nH20 > bnH2((NMESHP - 1) * 2)) {
        bEH2((NMESHP - 1) * 2) += 3.0 / 2.0 * Ta0 * (nH20 - bnH2((NMESHP - 1) * 2));
    } // 20220531
    else {
        bEH2((NMESHP - 1) * 2) = bEH2((NMESHP - 1) * 2) / bnH2((NMESHP - 1) * 2) * nH20;
    }
    bEH2(0 + 1) = (bEH2(2) - bEH2(0)) / (aR[1] - aR[0]);
    bEH2((NMESHP - 1) * 2 + 1) = (bEH2((NMESHP - 1) * 2) - bEH2((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    bnH2(0) = nH20;
    bnH2((NMESHP - 1) * 2) = nH20;
    bnH2(0 + 1) = (bnH2(2) - bnH2(0)) / (aR[1] - aR[0]);
    bnH2((NMESHP - 1) * 2 + 1) = (bnH2((NMESHP - 1) * 2) - bnH2((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    // bE1H2(0)              = bE1H2(0)/bn1H2(0)*nH20; // PSI
    // bE1H2((NMESHP-1)*2)   = bE1H2((NMESHP-1)*2)/bn1H2((NMESHP-1)*2)*nH20; // PSI
    if (nH20 > bn1H2(0)) {
        bE1H2(0) += 3.0 / 2.0 * Ta0 * (nH20 - bn1H2(0));
    } // 20220531
    else {
        bE1H2(0) = bE1H2(0) / bn1H2(0) * nH20;
    }
    if (nH20 > bn1H2((NMESHP - 1) * 2)) {
        bE1H2((NMESHP - 1) * 2) += 3.0 / 2.0 * Ta0 * (nH20 - bn1H2((NMESHP - 1) * 2));
    } // 20220531
    else {
        bE1H2((NMESHP - 1) * 2) = bE1H2((NMESHP - 1) * 2) / bn1H2((NMESHP - 1) * 2) * nH20;
    }
    bE1H2(0 + 1) = (bE1H2(2) - bE1H2(0)) / (aR[1] - aR[0]);
    bE1H2((NMESHP - 1) * 2 + 1) = (bE1H2((NMESHP - 1) * 2) - bE1H2((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    bn1H2(0) = nH20;
    bn1H2((NMESHP - 1) * 2) = nH20;
    bn1H2(0 + 1) = (bn1H2(2) - bn1H2(0)) / (aR[1] - aR[0]);
    bn1H2((NMESHP - 1) * 2 + 1) = (bn1H2((NMESHP - 1) * 2) - bn1H2((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    // SnH2(0)              = 0.0;
    // SnH2((NMESHP-1)*2)   = 0.0;
    // SnH2(0+1)            = (SnH2(2)-SnH2(0))/(aR[1]-aR[0]);
    // SnH2((NMESHP-1)*2+1) = (SnH2((NMESHP-1)*2)-SnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    // SEH2(0)              = 0.0;
    // SEH2((NMESHP-1)*2)   = 0.0;
    // SEH2(0+1)            = (SEH2(2)-SEH2(0))/(aR[1]-aR[0]);
    // SEH2((NMESHP-1)*2+1) = (SEH2((NMESHP-1)*2)-SEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //     }
    //     // Hi Edge conditions, diffusion
    //     lambda = pow( DionHFS*(2.0*pi*aR[0]/nlimiters) / ( 9.79e5*sqrt(1.0/1.0*Tr.Te[0]) * pow(1.0+(nevac/nr.ne[0]),-0.66) ) , 0.5);
    //     bnHi(0+1)            = +bnHi(0)/lambda;
    //
    //     lambda = pow( DionLFS*(2.0*pi*aR[NMESHP-1]/nlimiters) / ( 9.79e5*sqrt(1.0/1.0*Tr.Te[NMESHP-1]) * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) ) , 0.5);
    //     bnHi((NMESHP-1)*2+1)  = -bnHi((NMESHP-1)*2)/lambda;
    //
    //     // H2i Edge conditions, diffusion
    //     lambda = pow( DionHFS*(2.0*pi*aR[0]/nlimiters) / ( 9.79e5*sqrt(1.0/2.0*Tr.Te[0]) * pow(1.0+(nevac/nr.ne[0]),-0.66) ) , 0.5);
    //     bnH2i(0+1)            = +bnH2i(0)/lambda;
    //
    //     lambda = pow( DionLFS*(2.0*pi*aR[NMESHP-1]/nlimiters) / ( 9.79e5*sqrt(1.0/2.0*Tr.Te[NMESHP-1]) * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) ) , 0.5);
    //     bnH2i((NMESHP-1)*2+1)  = -bnH2i((NMESHP-1)*2)/lambda;
    //
    //     // H3i Edge conditions, diffusion
    //     lambda = pow( DionHFS*(2.0*pi*aR[0]/nlimiters) / ( 9.79e5*sqrt(1.0/3.0*Tr.Te[0]) * pow(1.0+(nevac/nr.ne[0]),-0.66) ) , 0.5);
    //     bnH3i(0+1)            = +bnH3i(0)/lambda;
    //
    //     lambda = pow( DionLFS*(2.0*pi*aR[NMESHP-1]/nlimiters) / ( 9.79e5*sqrt(1.0/3.0*Tr.Te[NMESHP-1]) * pow(1.0+(nevac/nr.ne[NMESHP-1]),-0.66) ) , 0.5);
    //     bnH3i((NMESHP-1)*2+1)  = -bnH3i((NMESHP-1)*2)/lambda;
    //
    //     SsnH2(0)            -= bnH2(0+1)                * DH2[0]        ;
    //     SsEH2(0)            -= gEd*bEH2(0+1)            * DH2[0]        ;
    //     SsnH2((NMESHP-1)*2) -= bnH2((NMESHP-1)*2+1)     * DH2[NMESHP-1] ;
    //     SsEH2((NMESHP-1)*2) -= gEd*bEH2((NMESHP-1)*2+1) * DH2[NMESHP-1] ;
    //
    //     vth = scale * sqrt((kb * TH2_array[0]        * 11600.0)/(2.0*ma))*100.0;
    //     SsnH2(0)            += 1.0/2.0*vth/3.0/DH2[0]*bnH2(0)                       * DH2[0]        * (1.0-RH2);
    //     SsEH2(0)            += gEd*1.0/2.0*vth/3.0/DH2[0]*bEH2(0)                   * DH2[0]        * (1.0-RH2 * (1.0 - ((1.0-REH2) * (1.0-Ta0/TH2_array[0]))));
    //     vth = scale * sqrt((kb * TH2_array[NMESHP-1] * 11600.0)/(2.0*ma))*100.0;
    //     SsnH2((NMESHP-1)*2) -= 1.0/2.0*vth/3.0/DH2[NMESHP-1]*bnH2((NMESHP-1)*2)     * DH2[NMESHP-1] * (1.0-RH2);
    //     SsEH2((NMESHP-1)*2) -= gEd*1.0/2.0*vth/3.0/DH2[NMESHP-1]*bEH2((NMESHP-1)*2) * DH2[NMESHP-1] * (1.0-RH2 * (1.0 - ((1.0-REH2) * (1.0-Ta0/TH2_array[NMESHP-1]))));
    //
    //     SsnH2(0)            -= 0.5 * bnH(0+1)                 * DH[0]; // Reflection coefficient already included in slope
    //     SsEH2(0)            -= gEd*0.5 * bnH(0+1)             * DH[0] * Ta0 * 3.0/2.0;
    //     SsnH2((NMESHP-1)*2) -= 0.5 * bnH((NMESHP-1)*2+1)      * DH[NMESHP-1];
    //     SsEH2((NMESHP-1)*2) -= gEd*0.5 * bnH((NMESHP-1)*2+1)  * DH[NMESHP-1] * Ta0 * 3.0/2.0;
    //
    //     SsnH2(0)            -= 0.5 * bnHi(0+1)                * DionHFS                   * RHi;
    //     SsEH2(0)            -= gEd*0.5 * bnHi(0+1)            * DionHFS * Ta0 * 3.0/2.0   * RHi;
    //     SsnH2((NMESHP-1)*2) -= 0.5 * bnHi((NMESHP-1)*2+1)     * DionLFS                   * RHi;
    //     SsEH2((NMESHP-1)*2) -= gEd*0.5 * bnHi((NMESHP-1)*2+1) * DionLFS * Ta0 * 3.0/2.0   * RHi;
    //
    //     SsnH2(0)            -= 1.0 * bnH2i(0+1)            * DionHFS                      * RHi;
    //     SsEH2(0)            -= gEd*1.0 * bnH2i(0+1)            * DionHFS * Ta0 * 3.0/2.0  * RHi;
    //     SsnH2((NMESHP-1)*2) -= 1.0 * bnH2i((NMESHP-1)*2+1) * DionLFS                      * RHi;
    //     SsEH2((NMESHP-1)*2) -= gEd*1.0 * bnH2i((NMESHP-1)*2+1) * DionLFS * Ta0 * 3.0/2.0  * RHi;
    //
    //     SsnH2(0)            -= 1.5 * bnH3i(0+1)            * DionHFS                      * RHi;
    //     SsEH2(0)            -= gEd*1.5 * bnH3i(0+1)            * DionHFS * Ta0 * 3.0/2.0  * RHi;
    //     SsnH2((NMESHP-1)*2) -= 1.5 * bnH3i((NMESHP-1)*2+1) * DionLFS                      * RHi;
    //     SsEH2((NMESHP-1)*2) -= gEd*1.5 * bnH3i((NMESHP-1)*2+1) * DionLFS * Ta0 * 3.0/2.0  * RHi;
    //
    //
    //     SsnH2((NMESHP-1)*2) += 0.5 * bnHi((NMESHP-1)*2) * VionLFS                   * RHi;
    //     SsEH2((NMESHP-1)*2) += 0.5 * bnHi((NMESHP-1)*2) * VionLFS * Ta0 * 3.0/2.0   * RHi;
    //
    //     SsnH2((NMESHP-1)*2) += 1.0 * bnH2i((NMESHP-1)*2) * VionLFS                  * RHi;
    //     SsEH2((NMESHP-1)*2) += 1.0 * bnH2i((NMESHP-1)*2) * VionLFS * Ta0 * 3.0/2.0  * RHi;
    //
    //     SsnH2((NMESHP-1)*2) += 1.5 * bnH3i((NMESHP-1)*2) * VionLFS                  * RHi;
    //     SsEH2((NMESHP-1)*2) += 1.5 * bnH3i((NMESHP-1)*2) * VionLFS * Ta0 * 3.0/2.0  * RHi;
    //
    //     //////////////////////////////////////////////////////
    //     //// Gasinjection at LFS & HFS + Pumping at LFS & HFS
    //     ////////////////////////////////////////////////////////
    //     if (binj) {
    //     SsnH2(0)            -= nH20*Vpl/(2*pi*aR[0]*150)        / (taupumpH2*2.0);
    //     SsEH2(0)            -= nH20*Vpl/(2*pi*aR[0]*150)        / (taupumpH2*2.0) * Ta0 * 3.0/2.0;
    //     SsnH2((NMESHP-1)*2) += nH20*Vpl/(2*pi*aR[NMESHP-1]*150) / (taupumpH2*2.0);
    //     SsEH2((NMESHP-1)*2) += nH20*Vpl/(2*pi*aR[NMESHP-1]*150) / (taupumpH2*2.0) * Ta0 * 3.0/2.0;
    //     }
    //     if (bpump) {
    //     SsnH2(0)            += bnH2(0) * Vpl/(2*pi*aR[0]*150) / (taupumpH2*2.0) * TH2_array[0]/Ta0;
    //     SsEH2(0)            += bEH2(0) * Vpl/(2*pi*aR[0]*150) / (taupumpH2*2.0) * TH2_array[0]/Ta0;
    //     SsnH2((NMESHP-1)*2) -= bnH2((NMESHP-1)*2) * Vpl/(2*pi*aR[NMESHP-1]*150) / (taupumpH2*2.0) * TH2_array[NMESHP-1]/Ta0;
    //     SsEH2((NMESHP-1)*2) -= bEH2((NMESHP-1)*2) * Vpl/(2*pi*aR[NMESHP-1]*150) / (taupumpH2*2.0) * TH2_array[NMESHP-1]/Ta0;
    //     }
    //
    // }

    // solve
    if (btranspneut) {
#pragma omp parallel sections
        {
#pragma omp section
            {
                A = Ts - coef1 * DHs * tstep; // SnH2  = -0.5*0.5*(SnHs*bnH);    SEH2  = -0.5*0.5*Ta0*3.0/2.0*(SnHs*bnH);
                //        bnH = (Ts + 0.5*DHs*tstep)*bnH + Ss*SnH*tstep; // piecewise linear source
                bnH = (Ts)*bnH * coef2 - Ts * bn1H * coef3 + coef1 * Ts * SnH * tstep;
                solver.compute(A);
                x = solver.solveWithGuess(bnH, bnH);
                bnH = x; // SnH2 += -0.5*0.5*(SnHs*bnH);    SEH2 += -0.5*0.5*Ta0*3.0/2.0*(SnHs*bnH);
            }
#pragma omp section
            {
                B = Ts - coef1 * gEdn * DHs * tstep;
                //        bEH = (Ts + 0.5*DHs*tstep)*bEH + Ss*SEH*tstep; // piecewise linear source
                bEH = (Ts)*bEH * coef2 - Ts * bE1H * coef3 + coef1 * Ts * SEH * tstep;
                solver2.compute(B);
                x2 = solver2.solveWithGuess(bEH, bEH);
                bEH = x2;
            }
#pragma omp section
            {
                C = Ts - coef1 * DH2s * tstep;
                //        bnH2 = (Ts + 0.5*DH2s*tstep)*bnH2 + Ss*SnH2*tstep - Es*SsnH2*tstep; // piecewise linear source
                bnH2 = (Ts)*bnH2 * coef2 - Ts * bn1H2 * coef3 + coef1 * Ts * SnH2 * tstep - coef1 * Es * SsnH2 * tstep;
                solver3.compute(C);
                x3 = solver3.solveWithGuess(bnH2, bnH2);
                bnH2 = x3;
            }
#pragma omp section
            {
                D = Ts - coef1 * gEdn * DH2s * tstep;
                //        bEH2 = (Ts + gEdn*0.5*DH2s*tstep)*bEH2 + Ss*SEH2*tstep - Es*SsEH2*tstep; // piecewise linear source
                bEH2 = (Ts)*bEH2 * coef2 - Ts * bE1H2 * coef3 + coef1 * Ts * SEH2 * tstep - coef1 * Es * SsEH2 * tstep;
                solver4.compute(D);
                x4 = solver4.solveWithGuess(bEH2, bEH2);
                bEH2 = x4;
            }
        }
    }
    // cout << "after:  " << 0.5*bnH(0+1) + 0.5*bnHi(0+1) + 1.0*bnH2i(0+1) + 1.5*bnH3i(0+1) << "   " << bnH2(0+1) << endl;

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        if (bEH(id * 2) / (3.0 / 2.0 * bnH(id * 2)) < Ta0) {
            bEH(id * 2) = 3.0 / 2.0 * bnH(id * 2) * Ta0;
        }
    }
    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        if (bEH2(id * 2) / (3.0 / 2.0 * bnH2(id * 2)) < Ta0) {
            bEH2(id * 2) = 3.0 / 2.0 * bnH2(id * 2) * Ta0;
        }
    }

    // if (pcst) {
    //     bEH2(0)              = 3.0/2.0*nH20*Ta0;
    //     bEH2((NMESHP-1)*2)   = 3.0/2.0*nH20*Ta0;
    //     bEH2(0+1)            = (bEH2(2)-bEH2(0))/(aR[1]-aR[0]);
    //     bEH2((NMESHP-1)*2+1) = (bEH2((NMESHP-1)*2)-bEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     bnH2(0)              = nH20;
    //     bnH2((NMESHP-1)*2)   = nH20;
    //     bnH2(0+1)            = (bnH2(2)-bnH2(0))/(aR[1]-aR[0]);
    //     bnH2((NMESHP-1)*2+1) = (bnH2((NMESHP-1)*2)-bnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    // }
    // else if (quick)
    // {
    //     // bEH2(0)              = 3.0/2.0*nH20*TH2_array[0];
    //     // bEH2((NMESHP-1)*2)   = 3.0/2.0*nH20*TH2_array[NMESHP-1];
    //     // bEH2(0+1)            = (bEH2(2)-bEH2(0))/(aR[1]-aR[0]);
    //     // bEH2((NMESHP-1)*2+1) = (bEH2((NMESHP-1)*2)-bEH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //     //
    //     // bnH2(0)              = nH20;
    //     // bnH2((NMESHP-1)*2)   = nH20;
    //     // bnH2(0+1)            = (bnH2(2)-bnH2(0))/(aR[1]-aR[0]);
    //     // bnH2((NMESHP-1)*2+1) = (bnH2((NMESHP-1)*2)-bnH2((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    // }
    // else if (ncst)
    // {
    // bEH2(0)              = bEH2(0)/bnH2(0)*nH20; // PSI
    // bEH2((NMESHP-1)*2)   = bEH2((NMESHP-1)*2)/bnH2((NMESHP-1)*2)*nH20; // PSI
    if (nH20 > bnH2(0)) {
        bEH2(0) += 3.0 / 2.0 * Ta0 * (nH20 - bnH2(0));
    } // 20220531
    else {
        bEH2(0) = bEH2(0) / bnH2(0) * nH20;
    }
    if (nH20 > bnH2((NMESHP - 1) * 2)) {
        bEH2((NMESHP - 1) * 2) += 3.0 / 2.0 * Ta0 * (nH20 - bnH2((NMESHP - 1) * 2));
    } // 20220531
    else {
        bEH2((NMESHP - 1) * 2) = bEH2((NMESHP - 1) * 2) / bnH2((NMESHP - 1) * 2) * nH20;
    }
    bEH2(0 + 1) = (bEH2(2) - bEH2(0)) / (aR[1] - aR[0]);                                                                // 20220531
    bEH2((NMESHP - 1) * 2 + 1) = (bEH2((NMESHP - 1) * 2) - bEH2((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]); // 20220531

    bnH2(0) = nH20;
    bnH2((NMESHP - 1) * 2) = nH20;
    bnH2(0 + 1) = (bnH2(2) - bnH2(0)) / (aR[1] - aR[0]);
    bnH2((NMESHP - 1) * 2 + 1) = (bnH2((NMESHP - 1) * 2) - bnH2((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);
    // }
}

void transpHe(double tstep) { // called in solver.cpp within timeStep.cpp in Tomator1D.cpp

    double DHeI[NMESHP];
    #ifdef debug
    	double DminHe = 10000000.0;
        double DmaxHe = -1.0;
    #endif

    DionHFS = Dion[0];
    DionLFS = Dion[NMESHP - 1];
    VionLFS = Vion[NMESHP - 1];

    // #pragma omp parallel for private(vth, mfp)
    for (int id = 0; id < NMESHP; ++id) {
        // double mfp;
        if (bDfix_neutr) {
        	DHeI[id] = Dfix_neut;
        }
        else {
        	vth = sqrt((kb * THeI_array[id] * 11600.0) / (4.0 * ma)) * 100.0;
        	mfp = vth / nuHeI_array[id];
        	DHeI[id] = 1.0 / 3.0 * vth * 1.0 / (1.0 / mfp + 1.0 / (a / 2.0)); // vth*a/2.0 ; // cm2/s
                                                                          // DHeI[id] = 1.0/3.0 * vth * a/2.0 ; // vth*a/2.0 ; // cm2/s
        	// cout << "He   " << mfp << "   " <<  DHeI[id] << endl;
        	#ifdef debug
        		if (DHeI[id]>DmaxHe){
				DmaxHe = DH[id];
			}
			if(DHeI[id]<DminHe){
				DminHe = DH[id];
			}
        	#endif
    	}
    }
    #ifdef debug
    	cout << "DHe_min = " << DminHe/1e4 << " [m2/s]; DHe_max =" << DmaxHe/1e4 << " [m2/s]" << endl;
    #endif

#pragma omp parallel for collapse(2)
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DLfbba[i] * DHeI[i - 1] + DLfbbb[i] * DHeI[i] + DRfbbc[i] * DHeI[i + 1] + DRfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DHeI[i - 1] + dDLfbbb[i] * DHeI[i] + dDRfbbc[i] * DHeI[i + 1] + dDRfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DHeI[i - 1] + DLdfbbb[i] * DHeI[i] + DRdfbbc[i] * DHeI[i + 1] + DRdfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DHeI[i - 1] + dDLdfbbb[i] * DHeI[i] + dDRdfbbc[i] * DHeI[i + 1] + dDRdfbbb[i] * DHeI[i];
            } else if (j == i - 1) {
                Ds.coeffRef(2 * i, 2 * j) = DLfaba[i] * DHeI[i - 1] + DLfabb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DHeI[i - 1] + dDLfabb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DHeI[i - 1] + DLdfabb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DHeI[i - 1] + dDLdfabb[i] * DHeI[i];
            } else if (j == i + 1) {
                Ds.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DHeI[i + 1] + DRfcbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DHeI[i + 1] + dDRfcbb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DHeI[i + 1] + DRdfcbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DHeI[i + 1] + dDRdfcbb[i] * DHeI[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DRfbbc[i] * DHeI[i + 1] + DRfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfbbc[i] * DHeI[i + 1] + dDRfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfbbc[i] * DHeI[i + 1] + DRdfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfbbc[i] * DHeI[i + 1] + dDRdfbbb[i] * DHeI[i];
            } else if (j == i + 1) {
                Ds.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DHeI[i + 1] + DRfcbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DHeI[i + 1] + dDRfcbb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DHeI[i + 1] + DRdfcbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DHeI[i + 1] + dDRdfcbb[i] * DHeI[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DLfbba[i] * DHeI[i - 1] + DLfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DHeI[i - 1] + dDLfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DHeI[i - 1] + DLdfbbb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DHeI[i - 1] + dDLdfbbb[i] * DHeI[i];
            } else if (j == i - 1) {
                Ds.coeffRef(2 * i, 2 * j) = DLfaba[i] * DHeI[i - 1] + DLfabb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DHeI[i - 1] + dDLfabb[i] * DHeI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DHeI[i - 1] + DLdfabb[i] * DHeI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DHeI[i - 1] + dDLdfabb[i] * DHeI[i];
            }
        }
    }

    if (nHeI0 > bnHeI(0)) {
        bEHeI(0) += 3.0 / 2.0 * Ta0 * (nHeI0 - bnHeI(0));
    } // 20220531
    else {
        bEHeI(0) = bEHeI(0) / bnHeI(0) * nHeI0;
    }
    if (nHeI0 > bnHeI((NMESHP - 1) * 2)) {
        bEHeI((NMESHP - 1) * 2) += 3.0 / 2.0 * Ta0 * (nHeI0 - bnHeI((NMESHP - 1) * 2));
    } // 20220531
    else {
        bEHeI((NMESHP - 1) * 2) = bEHeI((NMESHP - 1) * 2) / bnHeI((NMESHP - 1) * 2) * nHeI0;
    }
    bEHeI(0 + 1) = (bEHeI(2) - bEHeI(0)) / (aR[1] - aR[0]);
    bEHeI((NMESHP - 1) * 2 + 1) = (bEHeI((NMESHP - 1) * 2) - bEHeI((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    bnHeI(0) = nHeI0;
    bnHeI((NMESHP - 1) * 2) = nHeI0;
    bnHeI(0 + 1) = (bnHeI(2) - bnHeI(0)) / (aR[1] - aR[0]);
    bnHeI((NMESHP - 1) * 2 + 1) = (bnHeI((NMESHP - 1) * 2) - bnHeI((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    // bE1HeI(0)              = bE1HeI(0)/bn1HeI(0)*nHeI0;
    // bE1HeI((NMESHP-1)*2)   = bE1HeI((NMESHP-1)*2)/bn1HeI((NMESHP-1)*2)*nHeI0; // PSI
    if (nHeI0 > bn1HeI(0)) {
        bE1HeI(0) += 3.0 / 2.0 * Ta0 * (nHeI0 - bn1HeI(0));
    } // 20220531
    else {
        bE1HeI(0) = bE1HeI(0) / bn1HeI(0) * nHeI0;
    }
    if (nHeI0 > bn1HeI((NMESHP - 1) * 2)) {
        bE1HeI((NMESHP - 1) * 2) += 3.0 / 2.0 * Ta0 * (nHeI0 - bn1HeI((NMESHP - 1) * 2));
    } // 20220531
    else {
        bE1HeI((NMESHP - 1) * 2) = bE1HeI((NMESHP - 1) * 2) / bn1HeI((NMESHP - 1) * 2) * nHeI0;
    }
    bE1HeI(0 + 1) = (bE1HeI(2) - bE1HeI(0)) / (aR[1] - aR[0]);
    bE1HeI((NMESHP - 1) * 2 + 1) = (bE1HeI((NMESHP - 1) * 2) - bE1HeI((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    bn1HeI(0) = nHeI0;
    bn1HeI((NMESHP - 1) * 2) = nHeI0;
    bn1HeI(0 + 1) = (bn1HeI(2) - bn1HeI(0)) / (aR[1] - aR[0]);
    bn1HeI((NMESHP - 1) * 2 + 1) = (bn1HeI((NMESHP - 1) * 2) - bn1HeI((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    if (btranspneut) {
#pragma omp parallel sections
        {
#pragma omp section
            {
                A = Ts - coef1 * Ds * tstep;
                bnHeI = (Ts)*bnHeI * coef2 - Ts * bn1HeI * coef3 + coef1 * Ts * SnHeI * tstep - coef1 * Es * SsnHeI * tstep;
                solver.compute(A);
                x = solver.solveWithGuess(bnHeI, bnHeI);
                bnHeI = x;
            }
#pragma omp section
            {
                B = Ts - coef1 * gEdn * Ds * tstep;
                bEHeI = (Ts)*bEHeI * coef2 - Ts * bE1HeI * coef3 + coef1 * Ts * SEHeI * tstep - coef1 * Es * SsEHeI * tstep;
                solver2.compute(B);
                x2 = solver2.solveWithGuess(bEHeI, bEHeI);
                bEHeI = x2;
            }
        }
    }

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        if (bEHeI(id * 2) / (3.0 / 2.0 * bnHeI(id * 2)) < Ta0) {
            bEHeI(id * 2) = 3.0 / 2.0 * bnHeI(id * 2) * Ta0;
        }
    }

    // if (pcst) {
    //     bEHeI(0)              = 3.0/2.0*nHeI0*Ta0;
    //     bEHeI((NMESHP-1)*2)   = 3.0/2.0*nHeI0*Ta0;
    //     bEHeI(0+1)            = (bEHeI(2)-bEHeI(0))/(aR[1]-aR[0]);
    //     bEHeI((NMESHP-1)*2+1) = (bEHeI((NMESHP-1)*2)-bEHeI((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //
    //     bnHeI(0)              = nHeI0;
    //     bnHeI((NMESHP-1)*2)   = nHeI0;
    //     bnHeI(0+1)            = (bnHeI(2)-bnHeI(0))/(aR[1]-aR[0]);
    //     bnHeI((NMESHP-1)*2+1) = (bnHeI((NMESHP-1)*2)-bnHeI((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    // }
    // else if (quick) {
    //     // bEHeI(0)              = 3.0/2.0*nHeI0*THeI_array[0];
    //     // bEHeI((NMESHP-1)*2)   = 3.0/2.0*nHeI0*THeI_array[NMESHP-1];
    //     // bEHeI(0+1)            = (bEHeI(2)-bEHeI(0))/(aR[1]-aR[0]);
    //     // bEHeI((NMESHP-1)*2+1) = (bEHeI((NMESHP-1)*2)-bEHeI((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    //     //
    //     // bnHeI(0)              = nHeI0;
    //     // bnHeI((NMESHP-1)*2)   = nHeI0;
    //     // bnHeI(0+1)            = (bnHeI(2)-bnHeI(0))/(aR[1]-aR[0]);
    //     // bnHeI((NMESHP-1)*2+1) = (bnHeI((NMESHP-1)*2)-bnHeI((NMESHP-2)*2))/(aR[NMESHP-1]-aR[NMESHP-2]);
    // }
    // else if (ncst) {

    // bEHeI(0)              = bEHeI(0)/bnHeI(0)*nHeI0; // PSI
    // bEHeI((NMESHP-1)*2)   = bEHeI((NMESHP-1)*2)/bnHeI((NMESHP-1)*2)*nHeI0; // PSI
    if (nHeI0 > bnHeI(0)) {
        bEHeI(0) += 3.0 / 2.0 * Ta0 * (nHeI0 - bnHeI(0));
    } // 20220531
    else {
        bEHeI(0) = bEHeI(0) / bnHeI(0) * nHeI0;
    }
    if (nHeI0 > bnHeI((NMESHP - 1) * 2)) {
        bEHeI((NMESHP - 1) * 2) += 3.0 / 2.0 * Ta0 * (nHeI0 - bnHeI((NMESHP - 1) * 2));
    } // 20220531
    else {
        bEHeI((NMESHP - 1) * 2) = bEHeI((NMESHP - 1) * 2) / bnHeI((NMESHP - 1) * 2) * nHeI0;
    }
    bEHeI(0 + 1) = (bEHeI(2) - bEHeI(0)) / (aR[1] - aR[0]);
    bEHeI((NMESHP - 1) * 2 + 1) = (bEHeI((NMESHP - 1) * 2) - bEHeI((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    bnHeI(0) = nHeI0;
    bnHeI((NMESHP - 1) * 2) = nHeI0;
    bnHeI(0 + 1) = (bnHeI(2) - bnHeI(0)) / (aR[1] - aR[0]);
    bnHeI((NMESHP - 1) * 2 + 1) = (bnHeI((NMESHP - 1) * 2) - bnHeI((NMESHP - 2) * 2)) / (aR[NMESHP - 1] - aR[NMESHP - 2]);

    // }
}

void transpC(double tstep) { // called in solver.cpp within timeStep.cpp in Tomator1D.cpp

    double DCI[NMESHP];

    DionHFS = Dion[0];
    DionLFS = Dion[NMESHP - 1];
    VionLFS = Vion[NMESHP - 1];

    // #pragma omp parallel for private(vth, mfp)
    for (int id = 0; id < NMESHP; ++id) {
        // double mfp;
        vth = sqrt((kb * 0.1 * 11600.0) / (12.0 * ma)) * 100.0; // tau = a/vth, D=a^2/2/tau
        mfp = vth / nuCI_array[id];
        DCI[id] = 1.0 / 3.0 * vth * 1.0 / (1.0 / mfp + 1.0 / (a / 2.0)); // vth*a/2.0 ; // cm2/s
                                                                         // DCI[id] = 1.0/3.0 * vth * a/2.0 ; // vth*a/2.0 ; // cm2/s
    }

#pragma omp parallel for collapse(2)
    for (int i = 1; i < NMESHP - 1; ++i) { // row
        for (int j = 0; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DLfbba[i] * DCI[i - 1] + DLfbbb[i] * DCI[i] + DRfbbc[i] * DCI[i + 1] + DRfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DCI[i - 1] + dDLfbbb[i] * DCI[i] + dDRfbbc[i] * DCI[i + 1] + dDRfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DCI[i - 1] + DLdfbbb[i] * DCI[i] + DRdfbbc[i] * DCI[i + 1] + DRdfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DCI[i - 1] + dDLdfbbb[i] * DCI[i] + dDRdfbbc[i] * DCI[i + 1] + dDRdfbbb[i] * DCI[i];
            } else if (j == i - 1) {
                Ds.coeffRef(2 * i, 2 * j) = DLfaba[i] * DCI[i - 1] + DLfabb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DCI[i - 1] + dDLfabb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DCI[i - 1] + DLdfabb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DCI[i - 1] + dDLdfabb[i] * DCI[i];
            } else if (j == i + 1) {
                Ds.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DCI[i + 1] + DRfcbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DCI[i + 1] + dDRfcbb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DCI[i + 1] + DRdfcbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DCI[i + 1] + dDRdfcbb[i] * DCI[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < 1; ++i) {     // row
        for (int j = 0; j < 2; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DRfbbc[i] * DCI[i + 1] + DRfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfbbc[i] * DCI[i + 1] + dDRfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfbbc[i] * DCI[i + 1] + DRdfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfbbc[i] * DCI[i + 1] + dDRdfbbb[i] * DCI[i];
            } else if (j == i + 1) {
                Ds.coeffRef(2 * i, 2 * j) = DRfcbc[i] * DCI[i + 1] + DRfcbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDRfcbc[i] * DCI[i + 1] + dDRfcbb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DRdfcbc[i] * DCI[i + 1] + DRdfcbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDRdfcbc[i] * DCI[i + 1] + dDRdfcbb[i] * DCI[i];
            }
        }
    }
#pragma omp parallel for collapse(2)
    for (int i = NMESHP - 1; i < NMESHP; ++i) {     // row
        for (int j = NMESHP - 2; j < NMESHP; ++j) { // column
            if (i == j) {
                Ds.coeffRef(2 * i, 2 * j) = DLfbba[i] * DCI[i - 1] + DLfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfbba[i] * DCI[i - 1] + dDLfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfbba[i] * DCI[i - 1] + DLdfbbb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfbba[i] * DCI[i - 1] + dDLdfbbb[i] * DCI[i];
            } else if (j == i - 1) {
                Ds.coeffRef(2 * i, 2 * j) = DLfaba[i] * DCI[i - 1] + DLfabb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j) = dDLfaba[i] * DCI[i - 1] + dDLfabb[i] * DCI[i];
                Ds.coeffRef(2 * i, 2 * j + 1) = DLdfaba[i] * DCI[i - 1] + DLdfabb[i] * DCI[i];
                Ds.coeffRef(2 * i + 1, 2 * j + 1) = dDLdfaba[i] * DCI[i - 1] + dDLdfabb[i] * DCI[i];
            }
        }
    }

    // CII
    lambda = pow(DionHFS * (2.0 * pi * aR[0] / nlimiters) / (9.79e5 * sqrt(1.0 / 12.0 * Tr.Te[0]) * pow(1.0 + (nevac / nr.ne[0]), -0.66)), 0.5);
    bnCII(0 + 1) = +bnCII(0) / lambda;

    lambda = pow(DionLFS * (2.0 * pi * aR[NMESHP - 1] / nlimiters) / (9.79e5 * sqrt(1.0 / 12.0 * Tr.Te[NMESHP - 1]) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66)), 0.5);
    bnCII((NMESHP - 1) * 2 + 1) = -bnCII((NMESHP - 1) * 2) / lambda;

    // CIII
    lambda = pow(DionHFS * (2.0 * pi * aR[0] / nlimiters) / (9.79e5 * sqrt(2.0 / 12.0 * Tr.Te[0]) * pow(1.0 + (nevac / nr.ne[0]), -0.66)), 0.5);
    bnCIII(0 + 1) = +bnCIII(0) / lambda;

    lambda = pow(DionLFS * (2.0 * pi * aR[NMESHP - 1] / nlimiters) / (9.79e5 * sqrt(2.0 / 12.0 * Tr.Te[NMESHP - 1]) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66)), 0.5);
    bnCIII((NMESHP - 1) * 2 + 1) = -bnCIII((NMESHP - 1) * 2) / lambda;

    // CIV
    lambda = pow(DionHFS * (2.0 * pi * aR[0] / nlimiters) / (9.79e5 * sqrt(2.0 / 12.0 * Tr.Te[0]) * pow(1.0 + (nevac / nr.ne[0]), -0.66)), 0.5);
    bnCIV(0 + 1) = +bnCIV(0) / lambda;

    lambda = pow(DionLFS * (2.0 * pi * aR[NMESHP - 1] / nlimiters) / (9.79e5 * sqrt(2.0 / 12.0 * Tr.Te[NMESHP - 1]) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66)), 0.5);
    bnCIV((NMESHP - 1) * 2 + 1) = -bnCIV((NMESHP - 1) * 2) / lambda;
    // CIV
    lambda = pow(DionHFS * (2.0 * pi * aR[0] / nlimiters) / (9.79e5 * sqrt(2.0 / 12.0 * Tr.Te[0]) * pow(1.0 + (nevac / nr.ne[0]), -0.66)), 0.5);
    bnCV(0 + 1) = +bnCV(0) / lambda;

    lambda = pow(DionLFS * (2.0 * pi * aR[NMESHP - 1] / nlimiters) / (9.79e5 * sqrt(2.0 / 12.0 * Tr.Te[NMESHP - 1]) * pow(1.0 + (nevac / nr.ne[NMESHP - 1]), -0.66)), 0.5);
    bnCV((NMESHP - 1) * 2 + 1) = -bnCV((NMESHP - 1) * 2) / lambda;

    SsnCI(0) -= bnCI(0 + 1) * DCI[0];
    SsnCI((NMESHP - 1) * 2) -= bnCI((NMESHP - 1) * 2 + 1) * DCI[NMESHP - 1];

    SsnCI(0) -= 1.0 * bnCII(0 + 1) * DionHFS;
    SsnCI((NMESHP - 1) * 2) -= 1.0 * bnCII((NMESHP - 1) * 2 + 1) * DionLFS;

    SsnCI(0) -= 1.0 * bnCIII(0 + 1) * DionHFS;
    SsnCI((NMESHP - 1) * 2) -= 1.0 * bnCIII((NMESHP - 1) * 2 + 1) * DionLFS;

    SsnCI(0) -= 1.0 * bnCIV(0 + 1) * DionHFS;
    SsnCI((NMESHP - 1) * 2) -= 1.0 * bnCIV((NMESHP - 1) * 2 + 1) * DionLFS;

    SsnCI(0) -= 1.0 * bnCV(0 + 1) * DionHFS;
    SsnCI((NMESHP - 1) * 2) -= 1.0 * bnCV((NMESHP - 1) * 2 + 1) * DionLFS;

    // SsnCI(0) -= 1.0 * bnCII(0) * Vion;
    //
    // SsnCI(0) -= 1.0 * bnCIII(0) * Vion;

    SsnCI((NMESHP - 1) * 2) += 1.0 * bnCII((NMESHP - 1) * 2) * VionLFS;

    SsnCI((NMESHP - 1) * 2) += 1.0 * bnCIII((NMESHP - 1) * 2) * VionLFS;

    SsnCI((NMESHP - 1) * 2) += 1.0 * bnCIV((NMESHP - 1) * 2) * VionLFS;

    SsnCI((NMESHP - 1) * 2) += 1.0 * bnCV((NMESHP - 1) * 2) * VionLFS;

    if (btranspneut) {
        A = Ts - coef1 * Ds * tstep;
        //        bnCI = (Ts + 0.5*Ds*tstep)*bnCI + Ss*SnCI*tstep - Es*SsnCI*tstep; // piecewise linear source
        bnCI = (Ts)*bnCI * coef2 - Ts * bn1CI * coef3 + coef1 * Ts * SnCI * tstep - coef1 * Es * SsnCI * tstep;
        solver.compute(A);
        x = solver.solveWithGuess(bnCI, bnCI);
        bnCI = x;
    }
}

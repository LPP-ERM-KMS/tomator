#include "collisions.h"

void collisions() {

    double k = 0.0;
    double knn = 0.0;
    double L2 = 0.0;

    #pragma omp parallel for
    for (int im = 0; im < NMESHP; ++im) {

        dnr.dne[im] = dnr.dnH[im] = dnr.dnH2[im] = dnr.dnHi[im] = dnr.dnH2i[im] = dnr.dnH3i[im] = 0.0;
        dnr.dnHeI[im] = dnr.dnHeII[im] = dnr.dnHeIII[im] = 0.0;
        dnr.dnCI[im] = dnr.dnCII[im] = dnr.dnCIII[im] = dnr.dnCIV[im] = dnr.dnCV[im] = 0.0;

        dEr.dEe[im] = dEr.dEH[im] = dEr.dEH2[im] = dEr.dEHi[im] = dEr.dEH2i[im] = dEr.dEH3i[im] = 0.0;
        dEr.dEHeI[im] = dEr.dEHeII[im] = dEr.dEHeIII[im] = 0.0;

        colrateRF.nue[im] = colrateRF.nuHi[im] = colrateRF.nuH2i[im] = colrateRF.nuH3i[im] = colrateRF.nuHeII[im] = colrateRF.nuHeIII[im] = 0.0;

        // if(nr.ne[im] < 0.0)     {nr.ne[im] = 0.0; cout << "cor ne " << im << endl;}
        // if(nr.nH[im] < 0.0)     {nr.nH[im] = 0.0; cout << "cor nH " << im << endl;}
        // if(nr.nH2[im] < 0.0)    {nr.nH2[im] = 0.0; cout << "cor nH2 " << im << endl;}
        // if(nr.nHi[im] < 0.0)    {nr.nHi[im] = 0.0; cout << "cor nHi " << im << endl;}
        // if(nr.nH2i[im] < 0.0)   {nr.nH2i[im] = 0.0; cout << "cor nH2i " << im << endl;}
        // if(nr.nH3i[im] < 0.0)   {nr.nH3i[im] = 0.0; cout << "cor nH3i " << im << endl;}
        // if(nr.nHeI[im] < 0.0)   {nr.nHeI[im] = 0.0; cout << "cor nHeI " << im << endl;}
        // if(nr.nHeII[im] < 0.0)  {nr.nHeII[im] = 0.0; cout << "cor nHeII " << im << endl;}
        // if(nr.nHeIII[im] < 0.0) {nr.nHeIII[im] = 0.0; cout << "cor nHeIII " << im << endl;}
        // if(nr.nCI[im] < 0.0)    {nr.nCI[im] = 0.0; cout << "cor nCI " << im << endl;}
        // if(nr.nCII[im] < 0.0)   {nr.nCII[im] = 0.0; cout << "cor nCII " << im << endl;}
        // if(nr.nCIII[im] < 0.0)  {nr.nCIII[im] = 0.0; cout << "cor nCIII " << im << endl;}
        // if(nr.nCIV[im] < 0.0)   {nr.nCIV[im] = 0.0; cout << "cor nCIV " << im << endl;}
        // if(nr.nCV[im] < 0.0)    {nr.nCV[im] = 0.0; cout << "cor nCV " << im << endl;}

        //        if(Tr.Te[im] < 0.1)     {Tr.Te[im] = 0.1;}
        //        if(Tr.TH[im] < 0.1)     {Tr.TH[im] = 0.1;}
        //        if(Tr.TH2[im] < 0.1)    {Tr.TH2[im] = 0.1;}
        //        if(Tr.THi[im] < 0.1)    {Tr.THi[im] = 0.1;}
        //        if(Tr.TH2i[im] < 0.1)   {Tr.TH2i[im] = 0.1;}
        //        if(Tr.TH3i[im] < 0.1)   {Tr.TH3i[im] = 0.1;}
        //        if(Tr.THeI[im] < 0.1)   {Tr.THeI[im] = 0.1;}
        //        if(Tr.THeII[im] < 0.1)  {Tr.THeII[im] = 0.1;}
        //        if(Tr.THeIII[im] < 0.1) {Tr.THeIII[im] = 0.1;}
        //        if(Tr.Te[im] > 2e4)     {Tr.Te[im] = 2e4;}
        //        if(Tr.TH[im] > 2e4)     {Tr.TH[im] = 2e4;}
        //        if(Tr.TH2[im] > 2e4)    {Tr.TH2[im] = 2e4;}
        //        if(Tr.THi[im] > 2e4)    {Tr.THi[im] = 2e4;}
        //        if(Tr.TH2i[im] > 2e4)   {Tr.TH2i[im] = 2e4;}
        //        if(Tr.TH3i[im] > 2e4)   {Tr.TH3i[im] = 2e4;}
        //        if(Tr.THeI[im] > 2e4)   {Tr.THeI[im] = 2e4;}
        //        if(Tr.THeII[im] > 2e4)  {Tr.THeII[im] = 2e4;}
        //        if(Tr.THeIII[im] > 2e4) {Tr.THeIII[im] = 2e4;}
        if (Tr.Te[im] < 0.0) {
            Tr.Te[im] = 0.0;
            cout << "cor Te" << endl;
        }
        if (Tr.TH[im] < 0.0) {
            Tr.TH[im] = 0.0;
            cout << "cor TH" << endl;
        }
        if (Tr.TH2[im] < 0.0) {
            Tr.TH2[im] = 0.0;
            cout << "cor TH2" << endl;
        }
        if (Tr.THi[im] < 0.0) {
            Tr.THi[im] = 0.0;
            cout << "cor THi" << endl;
        }
        if (Tr.TH2i[im] < 0.0) {
            Tr.TH2i[im] = 0.0;
            cout << "cor TH2i" << endl;
        }
        if (Tr.TH3i[im] < 0.0) {
            Tr.TH3i[im] = 0.0;
            cout << "cor TH3i" << endl;
        }
        if (Tr.THeI[im] < 0.0) {
            Tr.THeI[im] = 0.0;
            cout << "cor THeI" << endl;
        }
        if (Tr.THeII[im] < 0.0) {
            Tr.THeII[im] = 0.0;
            cout << "cor THeII" << endl;
        }
        if (Tr.THeIII[im] < 0.0) {
            Tr.THeIII[im] = 0.0;
            cout << "cor THeIII" << endl;
        }
        //        if(Tr.Te[im] > 2e4)     {Tr.Te[im] = 2e4; cout << "cor Te >" << endl;}
        //        if(Tr.TH[im] > 2e4)     {Tr.TH[im] = 2e4; cout << "cor TH >" << endl;}
        //        if(Tr.TH2[im] > 2e4)    {Tr.TH2[im] = 2e4; cout << "cor TH2 >" << endl;}
        //        if(Tr.THi[im] > 2e4)    {Tr.THi[im] = 2e4; cout << "cor THi >" << endl;}
        //        if(Tr.TH2i[im] > 2e4)   {Tr.TH2i[im] = 2e4; cout << "cor TH2i >" << endl;}
        //        if(Tr.TH3i[im] > 2e4)   {Tr.TH3i[im] = 2e4; cout << "cor TH3i >" << endl;}
        //        if(Tr.THeI[im] > 2e4)   {Tr.THeI[im] = 2e4; cout << "cor THeI >" << endl;}
        //        if(Tr.THeII[im] > 2e4)  {Tr.THeII[im] = 2e4; cout << "cor THeII >" << endl;}
        //        if(Tr.THeIII[im] > 2e4) {Tr.THeIII[im] = 2e4; cout << "cor THeIII >" << endl;}
    }

    double ne, nH, nH2, nHi, nH2i, nH3i, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, nCV;
    double Te, TH, TH2, THi, TH2i, TH3i, THeI, THeII, THeIII;
    double dne, dnH, dnH2, dnHi, dnH2i, dnH3i, dnHeI, dnHeII, dnHeIII, dnCI, dnCII, dnCIII, dnCIV, dnCV;
    double dEe, dEH, dEH2, dEHi, dEH2i, dEH3i, dEHeI, dEHeII, dEHeIII;
    double nue, nuH, nuH2, nuHi, nuH2i, nuH3i, nuHeI, nuHeII, nuHeIII;

    /*#pragma omp parallel for private(k, knn,                                                                                           \
                                         ne, nH, nH2, nHi, nH2i, nH3i, nHeI, nHeII, nHeIII, nCI, nCII, nCIII, nCIV, nCV,               \
                                         Te, TH, TH2, THi, TH2i, TH3i, THeI, THeII, THeIII,                                            \
                                         dne, dnH, dnH2, dnHi, dnH2i, dnH3i, dnHeI, dnHeII, dnHeIII, dnCI, dnCII, dnCIII, dnCIV, dnCV, \
                                         dEe, dEH, dEH2, dEHi, dEH2i, dEH3i, dEHeI, dEHeII, dEHeIII,                                   \
                                         nue, nuH, nuH2, nuHi, nuH2i, nuH3i, nuHeI, nuHeII, nuHeIII)*/
    for (int im = 0; im < NMESHP; ++im) {
        k = 0.0;

        ne = nr.ne[im];
        Te = Tr.Te[im];
        nH = nr.nH[im];
        TH = Tr.TH[im];
        nH2 = nr.nH2[im];
        TH2 = Tr.TH2[im];
        nHi = nr.nHi[im];
        THi = Tr.THi[im];
        nH2i = nr.nH2i[im];
        TH2i = Tr.TH2i[im];
        nH3i = nr.nH3i[im];
        TH3i = Tr.TH3i[im];
        nHeI = nr.nHeI[im];
        THeI = Tr.THeI[im];
        nHeII = nr.nHeII[im];
        THeII = Tr.THeII[im];
        nHeIII = nr.nHeIII[im];
        THeIII = Tr.THeIII[im];
        nCI = nr.nCI[im];
        nCII = nr.nCII[im];
        nCIII = nr.nCIII[im];
        nCIV = nr.nCIV[im];
        nCV = nr.nCV[im];

        dne = 0.0;
        dEe = 0.0;
        nue = 0.0;
        dnH = 0.0;
        dEH = 0.0;
        nuH = 0.0;
        dnH2 = 0.0;
        dEH2 = 0.0;
        nuH2 = 0.0;
        dnHi = 0.0;
        dEHi = 0.0;
        nuHi = 0.0;
        dnH2i = 0.0;
        dEH2i = 0.0;
        nuH2i = 0.0;
        dnH3i = 0.0;
        dEH3i = 0.0;
        nuH3i = 0.0;
        dnHeI = 0.0;
        dEHeI = 0.0;
        nuHeI = 0.0;
        dnHeII = 0.0;
        dEHeII = 0.0;
        nuHeII = 0.0;
        dnHeIII = 0.0;
        dEHeIII = 0.0;
        nuHeIII = 0.0;
        dnCI = 0.0;
        dnCII = 0.0;
        dnCIII = 0.0;
        dnCIV = 0.0;
        dnCV = 0.0;

        // We shoud add everything inside each if to different functions but there are too many variables to pass...
        /////////////////////
        // Electron collisions with H and Hi (bH)
        /////////////////////
        if (bH) {
            if (nH2 > 0.0) // Go through this only when Hydrogen is present
            {
                // 1.) Reaction 2.1.1-2.1.4b Excitation (H_exc)
                if (Te > 0.6) {
                    k = 9.70346e-8 * pow((10.2 / Te), 0.92457) * exp(-10.2 / Te) / (0.01351 + (10.2 / Te)); // ok 17/02/2010
                    dEe += -k * ne * nH * 10.2;                                                             // ok 17/02/2010

                    nue += +k * nH;
                    if (nue < 0.0) {
                        cout << "nue 1 " << endl;
                    }
                }
                // 2.) Reaction 2.1.5-2.1.7 Ionization (H_ion)
                if (Te > 0.6) {
                    k = 2.91e-8 * pow((13.6 / Te), 0.39) * exp(-13.6 / Te) / (0.232 + 13.6 / Te) * ndamp(nH);
                    knn = k * ne * nH;
                    dnH += -knn;
                    dnHi += +knn;
                    dne += +knn;

                    dEH += -knn * TH;
                    dEHi += +knn * TH;
                    dEe += -knn * 13.6;

                    nuH += +k * ne;
                    nue += +k * nH; // if(nue < 0.0)     {cout << "nue 2 " << endl;}
                    // 3 body recombination           [170.8160 = 4.0*3.14*13.6]
                    k = 1.4804e-25 * pow((170.8160 / Te), 1.5) * exp(13.6 / Te) * k * ndamp(nHi);
                    knn = k * ne * ne * nHi;
                    dnHi += -knn;
                    dnH += +knn;
                    dne += -knn;

                    dEHi += -knn * THi;
                    dEH += +knn * THi;
                    dEe += +knn * (13.6 - Te / 3.0);

                    nuHi += +k * ne * ne;
                    nue += +k * ne * nHi; // if(nue < 0.0)     {cout << "nue 3 " << endl;}
                }
                // Reaction 2.1.8a Radiative recombination (H_rec)
                if (Te < 1000.0) {
                    k = 7.982e-11 / (sqrt(Te / 2.713e-4) * pow((1 + sqrt(Te / 2.713e-4)), (1.0 - 0.7480)) * pow((1 + sqrt(Te / 60.631)), (1.0 + 0.7480))) * ndamp(nHi);
                    knn = k * ne * nHi;
                    // k += 7.1e-20 * pow(Te/20,-15/3); // added 20/02/2014, estimated correction
                    dnHi += -knn;
                    dne += -knn;
                    dnH += +knn;

                    dEHi += -knn * THi;
                    dEH += +knn * THi;
                    // dEe    +=   - knn*Te/3.0; //CHECK diff with realTe
                    dEe += -knn * Te * 0.667 * (3.0 / 2.0 + // 20220530 hydhel -. Te * (3/2 + dln(k)/dnl(Te)) which is approximately Te for Te<10
                                                (log(k) - log(7.982e-11 / (sqrt(Te * 0.99 / 2.713e-4) * pow((1 + sqrt(Te * 0.99 / 2.713e-4)), (1.0 - 0.7480)) * pow((1 + sqrt(Te * 0.99 / 60.631)), (1.0 + 0.7480))) * ndamp(nHi))) / (log(Te) - log(0.99 * Te)));
                    nue += +k * nHi;
                    if (nue < 0.0) {
                        cout << "nue 4 " << endl;
                    }
                }
            }
        }
        // END bH

        /////////////////////
        // Electron collisions with H2, H2i and H3i (bH2)
        /////////////////////
        if (bH2) {
            if (nH2 > 0.0) {
                //                cout << "  " << dEe ;
                k = RRH2(ELAS, ne, Te); // 2016/02/24
                knn = k * ne * nH2;
                L2 = 4.0 * me * (2.0 * mi) / pow(me + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
                dEe += +knn * L2 * (TH2 - Te);                          // Yoon 2008
                dEH2 += +knn * L2 * (Te - TH2);                         // corrected? 2016/04/19
                nuH2 += +k * ne * sqrt((2 * me) / (me + (2.0 * mi)));   // almost no change in momentum for H2
                nue += +k * nH2;
                //                cout << "  " << dEe << endl ;

                // Reaction 2.2.1-2.2.4 (H2_exc)
                if (Te > 0.1) {
                    k = RR(REAC221a, 0.0, Te);
                    dEe += -k * ne * nH2 * 0.5;
                    nue += +k * nH2;
                    nuH2 += +k * ne * sqrt((2 * me) / (me + (2.0 * mi))); // almost no change in momentum for H2
                    k = RR(REAC221b, 0.0, Te);
                    dEe += -k * ne * nH2 * 1.0;
                    nue += +k * nH2;
                    nuH2 += +k * ne * sqrt((2 * me) / (me + (2.0 * mi))); // almost no change in momentum for H2
                    k = RR(REAC222, 0.0, Te);
                    dEe += -k * ne * nH2 * 12.1;
                    nue += +k * nH2;
                    nuH2 += +k * ne * sqrt((2 * me) / (me + (2.0 * mi))); // almost no change in momentum for H2
                    k = RR(REAC223, 0.0, Te);
                    dEe += -k * ne * nH2 * 12.4;
                    nue += +k * nH2;
                    nuH2 += +k * ne * sqrt((2 * me) / (me + (2.0 * mi))); // almost no change in momentum for H2
                    k = RR(REAC224, 0.0, Te);
                    dEe += -k * ne * nH2 * 12.7; // Check with Dirk
                    nue += +k * nH2;
                    nuH2 += +k * ne * sqrt((2 * me) / (me + (2.0 * mi))); // almost no change in momentum for H2
                }
                // Reaction 2.2.5-2.2.8 (H2_dis)
                k = RRH2(DISS, ne, Te) * ndamp(nH2);
                knn = k * ne * nH2;

                dnH2 += -knn;
                dnH += +knn * 2.0; // ok 17/02/2010

                dEe += -knn * 10.5;
                dEH2 += -knn * (TH2);
                dEH += +knn * (TH2 + 2.0 * 3.0); // ok 17/02/2010

                nuH2 += +k * ne;
                nue += +k * nH2;
                // Reaction 2.2.9 (H2_ion)
                k = RRH2(IONI, ne, Te) * ndamp(nH2);
                knn = k * ne * nH2;
                dnH2 += -knn;
                dne += +knn;
                dnH2i += +knn; // ok 17/02/2010

                dEe += -knn * 15.4;
                dEH2 += -knn * TH2;
                dEH2i += +knn * TH2; // ok 17/02/2010

                nuH2 += +k * ne;
                nue += +k * nH2;
                // Reaction recombination (H2i_rec)
                k = RRH2(RECO, ne, Te) * ndamp(nH2i);
                knn = k * ne * nH2i;

                dnH2i += -knn;
                dne += -knn;
                dnH2 += +knn; // ok 17/02/2010

                // dEe    +=   - knn*Te/3.0; // CHECK realTe // The reaction is more probable at lower energies...
                dEe += -knn * Te * 0.667 * 0.89; // 20220530 Hydhel
                dEH2i += -knn * TH2i;
                dEH2 += +knn * TH2i; // ok 4/03/2010

                nuH2i = nuH2i + k * ne;
                nue = nue + k * nH2i;
                // Reaction 2.2.10 (H2_dision)
                k = RR(REAC2210, 0.0, Te) * ndamp(nH2); // ok 23/10/2009
                knn = k * ne * nH2;

                dnH2 += -knn;
                dne += +knn;
                dnHi += +knn;
                dnH += +knn; // ok 24/10/2009

                dEe += -knn * (18.0);
                dEH2 += -knn * TH2;
                dEHi += +knn * (TH2 / 2.0 + 0.1);
                dEH += +knn * (TH2 / 2.0 + 0.1); // ok 24/10/2009

                nuH2 += +k * ne;
                nue += +k * nH2;
                // // Reaction 2.2.11
                //     k=RR('Reac2211',0,Te); // ok 23/10/2009
                //     dnH2i = dnH2i - k*ne*nH2i;
                //     dne   = dne   + k*ne*nH2i;
                //     dnHi  = dnHi  + k*ne*nH2i*2; // ok 24/10/2009
                //
                //     dEe   = dEe   - k*ne*nH2i*(15.5);
                //     dEH2i = dEH2i - k*ne*nH2i*TH2i;
                //     dEHi  = dEHi  + k*ne*nH2i*(TH2i+2*0.4); // ok 24/10/2009
                // Reaction 2.2.12 (k_H2i_dis)
                k = RR(REAC2212, 0.0, Te) * ndamp(nH2i); // ok 23/10/2009
                knn = k * ne * nH2i;
                dnH2i += -knn;
                dnHi += +knn;
                dnH += +knn; // ok 24/10/2009

                dEe += -knn * (10.5);
                dEH2i += -knn * TH2i;
                dEHi += +knn * (TH2i / 2.0 + 4.3);
                dEH += +knn * (TH2i / 2.0 + 4.3); // ok 24/10/2009

                nuH2i += +k * ne;
                nue += +k * nH2i;
                // Reaction 2.2.13 (k_H2i_disexc)
                k = RR(REAC2213, 0.0, Te) * ndamp(nH2i); // CORRECTED 2013!!!
                knn = k * ne * nH2i;
                dnH2i += -knn;
                dnHi += +knn;
                dnH += +knn; // ok 24/10/2009

                dEe += -knn * (17.5);
                dEH2i += -knn * TH2i;
                dEHi += +knn * (TH2i / 2.0 + 1.5);
                dEH += +knn * (TH2i / 2.0 + 1.5); // ok 24/10/2009

                nuH2i += +k * ne;
                nue += +k * nH2i;
                // Reaction 2.2.14 (k_H2i_disrec)
                k = RR(REAC2214, 0.0, Te) * ndamp(nH2i); // ok 23/10/2009
                knn = k * ne * nH2i;
                dne += -knn;
                dnH2i += -knn;
                dnH += +knn;
                dnH += +knn;

                // dEe    +=   - knn*Te/3.0;
                dEe += -knn * Te * 0.667 * (3.0 / 2.0 + // 20220603 hydhel -. based on H+ reco
                                            (log(k) - log(RR(REAC2214, 0.0, 0.99 * Te) * ndamp(nH2i))) / (log(Te) - log(0.99 * Te)));
                dEH2i += -knn * TH2i;
                dEH += +knn * TH2i;

                nuH2i = nuH2i + k * ne;
                nue = nue + k * nH2i;
                // Reaction 2.2.15 a and b (k_H3i_disrec)
                k = RR(REAC2215, 0.0, Te) * ndamp(nH3i); // ok 23/10/2009
                knn = k * ne * nH3i;
                dne += -knn;
                dnH3i += -knn;
                dnH += +knn * 3.0; // ok 24/10/2009

                // dEe    +=   - knn*Te/3.0;
                dEe += -knn * Te * 0.667 * (3.0 / 2.0 + // 20220603 hydhel -. based on H+ reco
                                            (log(k) - log(RR(REAC2215, 0.0, 0.99 * Te) * ndamp(nH3i))) / (log(Te) - log(0.99 * Te)));
                dEH3i += -knn * TH3i;
                dEH += +knn * (TH3i + Te / 3.0); // ok 24/10/2009

                nuH3i += +k * ne;
                nue += +k * nH3i;

                dne += -knn;
                dnH3i += -knn;
                dnH2 += +knn;
                dnH += +knn;

                // dEe    +=   - knn*Te/3.0;
                dEe += -knn * Te * 0.667 * (3.0 / 2.0 + // 20220603 hydhel -. based on H+ reco
                                            (log(k) - log(RR(REAC2215, 0.0, 0.99 * Te) * ndamp(nH3i))) / (log(Te) - log(0.99 * Te)));
                dEH3i += -knn * TH3i;
                dEH2 += +knn * (TH3i + Te / 3.0) * 2.0 / 3.0;
                dEH += +knn * (TH3i + Te / 3.0) * 1.0 / 3.0;

                nuH3i += +k * ne;
                nue += +k * nH3i;
                // Reaction 2.2.16 (k_H3i_dis)
                k = RR(REAC2216, 0.0, Te) * ndamp(nH3i); // ok 23/10/2009
                knn = k * ne * nH3i;
                dnH3i += -knn;
                dnH += +knn * 2.0;
                dnHi += +knn; // ok 24/10/2009

                dEe += -knn * 14.0;
                dEH3i += -knn * TH3i;
                dEH += +knn * (TH3i / 3.0 + 4.33) * 2.0;
                dEHi += +knn * (TH3i / 3.0 + 4.33); // ok 24/10/2009

                nuH3i += +k * ne;
                nue += +k * nH3i;
            }
        }
        // END bH2

        /////////////////////
        // Electron collisions with He, He+ and He++ (bHe)
        /////////////////////
        if (bHe) {
            if (nHeI > 0.0) {
                // Ionization HeI (2.3.9-2.3.12) (k_HeI_ion)
                k = RRHe(IHE1, 1.0e11, Te) * ndamp(nHeI); // cout << "   " << k ;
                knn = k * ne * nHeI;
                dnHeI += -knn;
                dnHeII += +knn;
                dne += +knn;

                dEHeI += -knn * THeI;
                dEHeII += +knn * THeI;

                nuHeI = nuHeI + k * ne;
                nue = nue + k * nHeI;
                // Ionization HeII (2.3.19) (k_HeII_ion)
                k = RRHe(IHE2, 1.0e11, Te) * ndamp(nHeII); // cout << "   " << k ;
                knn = k * ne * nHeII;
                // if (k<1.0e-20) {k=1.0e-20;}
                dnHeII += -knn;
                dnHeIII += +knn;
                dne += +knn;

                dEHeII += -knn * THeII;
                dEHeIII += +knn * THeII; // cm^3/s*1/cm^3*1/cm^3*eV = eV/cm^3/s

                nuHeII = nuHeII + k * ne;
                nue = nue + k * nHeII;
                if (nue < 0.0) {
                    cout << "nue 5 " << endl;
                }
                double Tx, exp1, exp2, exp3;
                if (bADAS) {

		  // new correction Anthony Piras 2025/03/21
                  // Based on SOLPS-ITER cooling rates (ADAS)
                  double MTe[] = {1.0, 1.10000000000000, 1.20000000000000, 1.30000000000000, 1.40000000000000, 1.50000000000000, 1.60000000000000, 1.80000000000000, 1.90000000000000, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 70.0, 100.0, 200.0, 400.0, 600.0, 800.0, 1000.0};
                  double MkHeI[] = {1.984999768100701e-21, 7.830301492706550e-21, 2.873864536275069e-20, 9.069880352562618e-20, 2.054501954928153e-19, 6.042770680912184e-19, 1.175721397793912e-18, 4.425157944950157e-18, 8.611389418838644e-18, 1.298915147023593e-17, 3.479388485970521e-16, 1.672105079422092e-15, 5.449208218367541e-15, 1.073672202663698e-14, 1.835856451483364e-14, 2.656284560319206e-14, 3.658624533759706e-14, 4.847243721543640e-14, 1.114995167681486e-13, 1.784238083376350e-13, 2.359742151261333e-13, 2.931127115669682e-13, 3.837425083558295e-13, 4.687381793070677e-13, 5.793279278464787e-13, 6.824744616682854e-13, 8.125717138145458e-13, 8.569270493201217e-13, 8.748341396129159e-13, 8.846837461834757e-13, 8.914663006219775e-13}; // MkHeI is the cooling rate due to processes starting from He0, units are [eV*m3/s]

                  double MkHeII[] = {1.057806689207504e-19, 1.109186764983655e-19, 1.158351540623689e-19, 1.205827364963467e-19, 1.251777567565165e-19, 1.295908890698588e-19, 1.339708714000346e-19, 1.427997515007597e-19, 1.478234259499439e-19, 1.530086843149849e-19, 1.170661212865412e-18, 1.735023579467549e-17, 1.526412262328358e-16, 5.448771293527139e-16, 1.422986059639174e-15, 2.706764281302123e-15, 4.741812019556716e-15, 7.709836703261552e-15, 2.857333292986805e-14, 5.693878244302156e-14, 8.292193026020404e-14, 1.110297886397282e-13, 1.546108925301941e-13, 1.982101620805143e-13, 2.528525815567847e-13, 3.057663970721758e-13, 3.793260164363811e-13, 4.064949819055127e-13, 4.072248670092556e-13, 3.990528169230037e-13, 3.896830021251033e-13}; // MkHeII is the cooling rate due to processes starting from He+, units are [eV*m3/s]

                  exp1 = interpolate1D(MTe, MkHeI, 31, Te);
                  exp1 = exp1*1.0e6; // conversion to eV*cm3/s
                  exp2 = interpolate1D(MTe, MkHeII, 31, Te);
                  exp2 = exp2*1.0e6; // conversion to eV*cm3/s
                  exp3 = 0.0;

                  dEe += - 2.0/3.0 * (exp1 * nHeI + exp2 * nHeII + exp3 * nHeIII) * ne; // eV*cm3/s * cm-3 * cm-3 = eV/(cm3*s)
                //cout << - (2.0/3.0)*( exp1*nHeI + exp2*nHeII + exp3*nHeIII ) * ne << endl;
                }
                else {
                  // Energy loss rate excitation and ionization electrons
                  // //     //(2.3.1-2.3.7) ()
                  // //     Eel=SaschaEeloss('He1',Te,ne);
                  // //     dEe     = dEe     - Eel*nHeI; // eV/cm^3/s
                  // //     // nue     = nue     ????????
                  // //     //(2.3.14-2.3.15)
                  // //     Eel=SaschaEeloss('He2',Te,ne);
                  // //     dEe     = dEe     - Eel*nHeII; // eV/cm^3/s
                  // //     // nue     = nue     ????????
                  // Cooling rate given in IAEA
                  double HeIa1, HeIa2, HeIa3, HeIa4, HeIa5, HeIa6;
                  double HeIIa1, HeIIa2, HeIIa3, HeIIa4, HeIIa5, HeIIa6; //, HeIIa7, HeIIa8;
                  // double HeIIIa1, HeIIIa2, HeIIIa3, HeIIIa4, HeIIIa5, HeIIIa6;
                  HeIa1 = .6623e3;
                  HeIa2 = .9476e-1;
                  HeIa3 = .7456;
                  HeIa4 = -.2592;
                  HeIa5 = 3.8098;
                  HeIa6 = .4026;
                  HeIIa1 = .3476e3;
                  HeIIa2 = .1214;
                  HeIIa3 = .7974;
                  HeIIa4 = .4819;
                  HeIIa5 = 1.4066;
                  HeIIa6 = -.3639e-2;
                  // HeIIa7=.9720e-3;    HeIIa8=.4078;
                  // HeIIIa1=.2730e-1;   HeIIIa2=.6004;      HeIIIa3=-.2772e-1;
                  // HeIIIa4=.5925;      HeIIIa5=.3060e-2;	HeIIIa6=.3590;
                  Tx = Te / 1000.0;
                  // Energy loss on He
                  exp1 = (HeIa1 * exp(-HeIa2 / pow(Tx, HeIa3)) / (pow(Tx, HeIa4) + HeIa5 * pow(Tx, HeIa6))) * 1.0e-33 * 1.0e6; // Wcm3
                  // Energy loss on He+
                  exp2 = (HeIIa1 * exp(-HeIIa2 / pow(Tx, HeIIa3)) / (pow(Tx, HeIIa4) + HeIIa5 * pow(Tx, HeIIa6))) * 1.0e-33 * 1.0e6; // Wcm3 (IAEA without recombination losses)  // added 20/02/2014, estimated correction
                  // Energy loss on He++
                  exp3 = 0.0; // Wcm3 (IAEA without recombination losses)  // added 20/02/2014, estimated correction
                  //cout << -(2.0/3.0)* (exp1 * nHeI + exp2 * nHeII + exp3 * nHeIII) * ne/qe << endl;
                  dEe += -2.0 / 3.0 * (exp1 * nHeI + exp2 * nHeII + exp3 * nHeIII) * ne / qe;
                  /*// Energy loss on He
                   exp1=( HeIa1*exp(-HeIa2/pow(Tx,HeIa3))/(pow(Tx,HeIa4)+HeIa5*pow(Tx,HeIa6)) )*1.0e-33*1.0e6; //Wcm3
                   // Energy loss on He+
                   exp2=( HeIIa1*exp(-HeIIa2/pow(Tx,HeIIa3))/(pow(Tx,HeIIa4)+HeIIa5*pow(Tx,HeIIa6)) + HeIIa7*pow(Tx,HeIIa8) )*1.0e-33*1.0e6; //Wcm3
                   // Energy loss on He++
exp3=( HeIIIa1*pow(Tx,HeIIIa2)+HeIIIa3*pow(Tx,HeIIIa4)+HeIIIa5*pow(Tx,HeIIIa6) )*1.0e-33*1.0e6; //Wcm3 */

                }
                // Recombination HeII (2.3.13) (k_HeII_rec)
                k = RRHe(RHE2, 1.0e11, Te) * ndamp(nHeII); //  *  pow(1.0+1.0e1/nHeII+1.0e1/ne,-0.66) ; // cout << "   " << k ;//25/04/2017
                knn = k * ne * nHeII;
                dnHeI += +knn;  // cout << "   " << dnHeII ;
                dnHeII += -knn; // cout << "   " << dnHeII ;
                dne += -knn;

                dEHeI += +knn * THeII;
                dEHeII += -knn * THeII;
                // dEe     +=    - knn*Te/3.0; // added 20/02/2014
                dEe += -knn * Te * 0.667 * (3.0 / 2.0 + // 20220603 hydhel -. based on H+ reco
                                            (log(k) - log(RRHe(RHE2, 1.0e11, 0.99 * Te) * ndamp(nHeII))) / (log(Te) - log(0.99 * Te)));

                nuHeII = nuHeII + k * ne;
                nue = nue + k * nHeII;
                // Recombination HeIII (k_HeIII_rec)
                k = RRHe(RHE3, 1.0e11, Te) * ndamp(nHeIII); //  *  pow(1.0+1.0e1/nHeIII+1.0e1/ne,-0.66) ; // cout << "   " << k << endl; //25/04/2017
                knn = k * ne * nHeIII;
                dnHeII += +knn; // cout << "   " << dnHeII ;
                dnHeIII += -knn;
                dne += -knn;

                dEHeII += +knn * THeIII;
                dEHeIII += -knn * THeIII;
                // dEe     +=    - knn*Te/3.0; // added 20/02/2014
                dEe += -knn * Te * 0.667 * (3.0 / 2.0 + // 20220603 hydhel -. based on H+ reco
                                            (log(k) - log(RRHe(RHE3, 1.0e11, 0.99 * Te) * ndamp(nHeIII))) / (log(Te) - log(0.99 * Te)));
                nuHeIII = nuHeIII + k * ne;
                nue = nue + k * nHeIII;
                if (nue < 0.0) {
                    cout << "nue 6 " << endl;
                }
            }
        }
        // END bHe

        /////////////////////
        //// Charge exchange reactions (bcx)
        /////////////////////
        if (bcx) {
            // if p(6)==0  //&& t>4// Include cx reactions
            // Reaction 3.1.8-3.1.11 Charge exchange (k_HiH_cx)
            if ((THi + TH) / 2.0 > 0.1)
                k = 7.829e-9 * pow((THi + TH) / 2.0, 0.41) * ndamp(nHi); // ok 17/02/2010
            else {
                k = 0.0;
            }
            knn = k * nHi * nH;
            dEHi += +knn * (TH - THi);
            dEH += +knn * (THi - TH); // ok 24/10/2009

            nuH = nuH + k * nHi;
            nuHi = nuHi + k * nH;
            // Reaction 3.2.3 (k_HiH2_cx)
            // if p(8)==0 // Include molecular H2
            k = RRion(REAC323, THi, TH2) * ndamp(nHi); // ok 23/10/2009
            knn = k * nHi * nH2;
            dnHi += -knn;
            dnH2 += -knn;
            dnH += +knn;
            dnH2i += +knn; // ok 24/10/2009

            dEHi += -knn * THi;
            dEH2 += -knn * TH2;
            dEH += +knn * (THi - (1.83));
            dEH2i += +knn * TH2; // ok 24/10/2009

            nuH2 = nuH2 + k * nHi;
            nuHi = nuHi + k * nH2;
            // Reaction 4.3.1 (k_H2iH2_cx)
            // if p(8)==0 // Include molecular H2
            k = RRion(REAC431, TH2i, TH2) * ndamp(nH2i); // ok 23/10/2009
            knn = k * nH2i * nH2;
            dEH2i += +knn * (TH2 - TH2i);
            dEH2 += +knn * (TH2i - TH2); // ok 24/10/2009

            nuH2 = nuH2 + k * nH2i;
            nuH2i = nuH2i + k * nH2;
            // He+ + H(1s)  - > He + H+ ???? Reaction (k_HeIIH_cx)
            k = RRCX(CXHe2H, THeII, TH) * ndamp(nHeII); //   *  pow(1.0+1.0e1/nHeII+1.0e1/nH,-0.66) ;
            knn = k * nHeII * nH;
            dnHeII += -knn; // cout << "   " << dnHeII ;
            dnHeI += +knn;
            dnH += -knn;
            dnHi += +knn;

            dEHeII += -knn * THeII;
            dEHeI += +knn * THeII;
            dEH += -knn * TH;
            dEHi += +knn * TH;

            nuH = nuH + k * nHeII;
            nuHeII = nuHeII + k * nH;
            // He+ + He  - > He + He+ Reaction (k_HeIIHeI_cx)
            k = RRion(REAC531, THeII, THeI) * ndamp(nHeII);
            knn = k * nHeII * nHeI;
            dEHeII += -knn * (THeII - THeI);
            dEHeI += +knn * (THeII - THeI);

            nuHeII = nuHeII + k * nHeI;
            nuHeI = nuHeI + k * nHeII;
            // He++ + H(1s) - > He+ + H+ (k_HeIIIH_cx)
            k = RRCX(CXHe3H, THeIII, TH) * ndamp(nHeIII); //    *  pow(1.0+1.0e1/nHeII+1.0e1/nH,-0.66) ;
            knn = k * nHeIII * nH;
            dnHeIII += -knn;
            dnHeII += +knn; // cout << "   " << dnHeII ;
            dnH += -knn;
            dnHi += +knn;

            dEHeIII += -knn * THeIII;
            dEHeII += +knn * THeIII;
            dEH += -knn * TH;
            dEHi += +knn * TH;

            nuH = nuH + k * nHeIII;
            nuHeIII = nuHeIII + k * nH;
            // He++ + He(1s^2) . (He+)  + (He+)  (k_HeIIIHeI_cxa)
            k = RRCX(CXHe3He1, THeIII, THeI) * ndamp(nHeIII); //   *  pow(1.0+1.0e1/nHeIII+1.0e1/nHeI,-0.66) ; // cout << "   " << k << endl;
            knn = k * nHeIII * nHeI;
            dnHeIII += -knn;
            dnHeII += +knn; // cout << "   " << dnHeII ;
            dnHeI += -knn;
            dnHeII += +knn;

            dEHeIII += -knn * THeIII;
            dEHeII += +knn * THeIII;
            dEHeI += -knn * THeI;
            dEHeII += +knn * THeI;

            nuHeI = nuHeI + k * nHeIII;
            nuHeIII = nuHeIII + k * nHeI;
            // He++ + He . He  + He++   (k_HeIIIHeI_cxb)
            k = RRion(REAC631, THeIII, THeI) * ndamp(nHeIII);
            knn = k * nHeIII * nHeI;
            dEHeIII += -knn * (THeIII - THeI);
            dEHeI += +knn * (THeIII - THeI);

            nuHeI = nuHeI + k * nHeIII;
            nuHeIII = nuHeIII + k * nHeI;
        }
        // END bcx

        /////////////////////
        //// Other considered reactions
        /////////////////////
        if (bion) {
            // if p(7)==0  // && t>8// Include other ion reactions
            // Reaction 3.1.1 (k_HiH_exca)
            k = RRion(REAC311, THi, TH) * ndamp(nHi); // ok 23/10/2009

            dEHi += -k * nHi * nH * 10.2; // ok 24/10/2009

            nuH = nuH + k * nHi;
            nuHi = nuHi + k * nH;
            // Reaction 3.1.2 (k_HiH_excb)
            k = RRion(REAC312, THi, TH) * ndamp(nHi); // ok 23/10/2009

            dEHi += -k * nHi * nH * 10.2; // ok 24/10/2009

            nuH = nuH + k * nHi;
            nuHi = nuHi + k * nH;
            // Reaction 3.1.3
            // Reaction 3.2.1 // may be excluded: induces only a small correction in final
            // temperature, but causes huge troubles in breakdown phase of model, and
            // also in ramp down of the power
            k = RRion(REAC321, THi, TH2) * ndamp(nHi);
            if (ne < 1.0e9) {
                k = sin(ne / 1.0e9) * k; // to avoid numerical troubles at breakdown phase
            }
            dEHi += -k * nHi * nH2 * 0.1;

            nuH2 = nuH2 + k * nHi;
            nuHi = nuHi + k * nH2;
            // Reaction 3.2.2 (k_HiH2_excb)
            k = RRion(REAC322, THi, TH2) * ndamp(nHi);

            dEHi += -k * nHi * nH2 * 1.0;

            nuH2 = nuH2 + k * nHi;
            nuHi = nuHi + k * nH2;
            // Reaction 3.1.6 (k_HiH_ion)
            k = RRion(REAC316, THi, TH) * ndamp(nHi); // check 2013
            knn = k * nHi * nH;
            dnH += -knn;
            dnHi += +knn;
            dne += +knn;

            dEH += -knn * TH;
            dEHi += +knn * (TH - 13.6);

            nuH = nuH + k * nHi;
            nuHi = nuHi + k * nH;
            // Reaction 3.3.2 (k_HiHeI_ion)
            k = RRion(REAC332, THi, THeI) * ndamp(nHi); // check 2013
            knn = k * nHi * nHeI;
            dnHeI += -knn;
            dnHeII += +knn; // cout << "   " << dnHeII ;
            dne += +knn;

            dEHeI += -knn * THeI;
            dEHeII += +knn * THeI;
            dEHi += -knn * 24.58;

            nuHeI = nuHeI + k * nHi;
            nuHi = nuHi + k * nHeI;
            // Reaction 3.2.5
            k = RRion(REAC325, THi, TH2) * ndamp(nHi); // 05/2016
            knn = k * nHi * nH2;
            dne = dne + knn;
            dnH2 = dnH2 - knn;   //
            dnH2i = dnH2i + knn; //

            dEHi = dEHi - knn * 15.4;
            dEH2 = dEH2 - knn * TH2;
            dEH2i = dEH2i + knn * TH2; //

            nuHi = nuHi + k * nH2;
            nuH2 = nuH2 + k * nHi;
            // Reaction 3.2.6
            k = RRion(REAC326, THi, TH2i) * ndamp(nHi) * ndamp(nH2i); // 05/2016
            knn = k * nHi * nH2i;
            dnH2i = dnH2i - knn;
            dnH = dnH + knn;
            dnHi = dnHi + knn;

            dEHi = dEHi - knn * 10.5;
            dEHi = dEHi + knn * (TH2i / 2 + 4.5);
            dEH2i = dEH2i - knn * TH2i;
            dEH = dEH + knn * (TH2i / 2 + 4.5);

            nuHi = nuHi + k * nH2;
            nuH2i = nuH2i + k * nHi;
            // Reaction 4.2.1
            //     k=RR('Reac421',TH,TH2i); // ok 23/10/2009
            //     dnH2i = dnH2i - k*nH2i*nH;
            //     dnH   = dnH   + k*nH2i*nH;
            //     dnHi  = dnHi  + k*nH2i*nH; // ok 24/10/2009
            //
            //     dEH   = dEH   - k*nH2i*nH*11;
            //     dEH2i = dEH2i - k*nH2i*nH*TH2i;
            //     dEH   = dEH   + k*nH2i*nH*(TH2i/2+4.5);
            //     dEHi  = dEHi  + k*nH2i*nH*(TH2i/2+4.5); // ok 24/10/2009
            // Reaction 4.3.2
            //     k=RR('Reac432',TH2i,TH2);
            //     dnH2  = dnH2  - k*nH2i*nH2;
            //     dne   = dne   + k*nH2i*nH2;
            //     dnH2i = dnH2i + k*nH2i*nH2;
            //
            //     dEH2  = dEH2  - k*nH2i*nH2*TH2;
            //     dEH2i = dEH2i + k*nH2i*nH2*(Te-15.4);
            // Reaction 4.3.3 (k_H2iH2_H3i)
            // if p(8)==0 // Include molecular H2
            k = RRion(REAC433, TH2i, TH2) * ndamp(nH2i); // ok 23/10/2009
            knn = k * nH2i * nH2;
            dnH2i += -knn;
            dnH2 += -knn;
            dnH3i += +knn;
            dnH += +knn; // ok 24/10/2009

            dEH2i += -knn * TH2i;
            dEH2 += -knn * TH2;
            dEH3i += +knn * ((TH2i + TH2) / 2.0 + 0.585);
            dEH += +knn * ((TH2i + TH2) / 2.0 + 0.585); // ok 24/10/2009

            nuH2 = nuH2 + k * nH2i;
            nuH2i = nuH2i + k * nH2;
            // Reaction 4.4.1 (k_H2iHeI_ion)
            //     if p(8)==0 // Include molecular H2
            //     k=RRion('Reac441',TH2i,THeI);
            //     dnHeI   = dnHeI   - k*nH2i*nHeI;
            //     dnHeII  = dnHeII  + k*nH2i*nHeI;
            //     dne     = dne     + k*nH2i*nHeI;
            //
            //     dEH2i   = dEH2i   - k*nH2i*nHeI*24.58;
            //     dEHeI   = dEHeI   - k*nH2i*nHeI*THeI;
            //     dEHeII  = dEHeII  + k*nH2i*nHeI*THeI;
            //
            //     nuHeI   = nuHeI   + k*nH2i;
            //     nuH2i   = nuH2i   + k*nHeI;
            //     end
            // Reaction 5.2.1 (k_HeIIH2_ion)
            //     if p(8)==0 // Include molecular H2
            //     k=RRion('Reac521',THeII,TH2);
            //     dnH2    = dnH2    - k*nHeII*nH2;
            //     dnH2i   = dnH2i   + k*nHeII*nH2;
            //     dne     = dne     + k*nHeII*nH2;
            //
            //     dEHeII  = dEHeII  - k*nHeII*nH2*15.4;
            //     dEH2    = dEH2    - k*nHeII*nH2*TH2;
            //     dEH2i   = dEH2i   + k*nHeII*nH2*TH2;
            //
            //     nuHeII  = nuHeII  + k*nH2;
            //     nuH2    = nuH2    + k*nHeII;
            //     end
            // Reaction 5.2.3 (k_HeIIH2_cxdis)
            // if p(8)==0 // Include molecular H2
            k = RRion(REAC523, THeII, TH2) * ndamp(nHeII); //   *  pow(1.0+1.0e1/nHeII+1.0e1/nH2,-0.66) ;;
            knn = k * nHeII * nH2;
            dnHeII += -knn; // cout << "   " << dnHeII ; //<< endl;
            dnH2 += -knn;
            dnH += +knn;
            dnHi += +knn;
            dnHeI += +knn;

            dEHeII += -knn * THeII;
            dEH2 += -knn * TH2;
            dEH += +knn * 1.85;
            dEHi += +knn * 1.85;
            dEHeI += +knn * 3.1;

            nuHeII = nuHeII + k * nH2;
            nuH2 = nuH2 + k * nHeII;
            // Reaction 5.3.2 (k_HeIIHeI_ion)
            //     k=RRion('Reac532',THeII,THeI);
            //     dnHeII  = dnHeII  + k*nHeII*nHeI;
            //     dnHeI   = dnHeI   - k*nHeII*nHeI;
            //     dne     = dne     + k*nHeII*nHeI;
            //
            //     dEHeII  = dEHeII  + k*nHeII*nHeI*(THeI-24.58);
            //     dEHeI   = dEHeI   - k*nHeII*nHeI*THeI;
            //
            //     nuHeII  = nuHeII  + k*nHeI;
            //     nuHeI   = nuHeI   + k*nHeII;
            // Reaction 6.2.1 (k_HeIIIH2_ion)
            //     if p(8)==0 // Include molecular H2
            //     k=RRion('Reac621',THeIII,TH2);
            //     dnH2    = dnH2    - k*nHeII*nH2;
            //     dnH2i   = dnH2i   + k*nHeII*nH2;
            //     dne     = dne     + k*nHeII*nH2;
            //
            //     dEHeII  = dEHeII  + k*nHeII*nHeI*15.4;
            //     dEH2    = dEH2    - k*nHeII*nH2*TH2;
            //     dEH2i   = dEH2i   + k*nHeII*nH2*TH2;
            //     end
        }
        // END bion

        /////////////////////
        //// Elastic collisions (belas)
        /////////////////////
        if (belas) {
            // if p(9)==0 // && t>6 // Include elastic collisions
            // Hi + H
            k = RRel(HiH, THi, TH) * ndamp(nHi);
            L2 = 1.0;
            knn = k * nHi * nH;
            dEHi += -knn * L2 * (THi - TH);
            dEH += +knn * L2 * (THi - TH);
            nuHi = nuHi + k * nH;
            nuH = nuH + k * nHi;
            // H2i + H
            k = RRel(H2iH, TH2i, TH) * ndamp(nH2i);
            L2 = 4.0 * mi * (2.0 * mi) / pow(mi + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nH2i * nH;
            dEH2i += -knn * L2 * (TH2i - TH);
            dEH += +knn * L2 * (TH2i - TH);
            nuH2i = nuH2i + k * nH;
            nuH = nuH + k * nH2i;
            // H3i + H
            k = RRel(H3iH, TH3i, TH) * ndamp(nH3i);
            L2 = 4.0 * mi * (3.0 * mi) / pow(mi + (3.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nH3i * nH;
            dEH3i += -knn * L2 * (TH3i - TH);
            dEH += +knn * L2 * (TH3i - TH);
            nuH3i = nuH3i + k * nH;
            nuH = nuH + k * nH3i;
            // HeII + H
            k = RRel(HeIIH, THeII, TH) * ndamp(nHeII);
            L2 = 4.0 * mi * (4.0 * mi) / pow(mi + (4.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nHeII * nH;
            dEHeII += -knn * L2 * (THeII - TH);
            dEH += +knn * L2 * (THeII - TH);
            nuHeII = nuHeII + k * nH;
            nuH = nuH + k * nHeII;
            // Hi + H2
            k = RRel(HiH2, THi, TH2) * ndamp(nHi);
            L2 = 4.0 * mi * (2.0 * mi) / pow(mi + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nHi * nH2;
            dEHi += -knn * L2 * (THi - TH2);
            dEH2 += +knn * L2 * (THi - TH2);
            nuHi = nuHi + k * nH2;
            nuH2 = nuH2 + k * nHi;
            // H2i + H2
            k = RRel(H2iH2, TH2i, TH2) * ndamp(nH2i);
            L2 = 1.0; // Langevin’s energy loss parameter
            knn = k * nH2i * nH2;
            dEH2i += -knn * L2 * (TH2i - TH2);
            dEH2 += +knn * L2 * (TH2i - TH2);
            nuH2i = nuH2i + k * nH2;
            nuH2 = nuH2 + k * nH2i;
            // H3i + H2
            k = RRel(H3iH2, TH3i, TH2) * ndamp(nH3i);
            L2 = 4.0 * (3.0 * mi) * (2.0 * mi) / pow((3.0 * mi) + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nH3i * nH2;
            dEH3i += -knn * L2 * (TH3i - TH2);
            dEH2 += +knn * L2 * (TH3i - TH2);
            nuH3i = nuH3i + k * nH2;
            nuH2 = nuH2 + k * nH3i;
            // HeII + H2
            k = RRel(HeIIH2, THeII, TH2) * ndamp(nHeII);
            L2 = 4.0 * (4.0 * mi) * (2.0 * mi) / pow((4.0 * mi) + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nHeII * nH2;
            dEHeII += -knn * L2 * (THeII - TH2);
            dEH2 += +knn * L2 * (THeII - TH2);
            nuHeII = nuHeII + k * nH2;
            nuH2 = nuH2 + k * nHeII;
            // Hi + HeI
            k = RRel(HHe, TH, THeI);
            L2 = 4.0 * mi * (4.0 * mi) / pow(mi + (4.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nH * nHeI;
            dEH += -knn * L2 * (TH - THeI);
            dEHeI += +knn * L2 * (TH - THeI);
            nuH = nuH + k * nHeI;
            nuHeI = nuHeI + k * nH;
            // Hi + HeI
            k = RRel(HiHe, THi, THeI) * ndamp(nHi);
            L2 = 4.0 * mi * (4.0 * mi) / pow(mi + (4.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nHi * nHeI;
            dEHi += -knn * L2 * (THi - THeI);
            dEHeI += +knn * L2 * (THi - THeI);
            nuHi = nuHi + k * nHeI;
            nuHeI = nuHeI + k * nHi;
            // HeII + HeI
            k = RRel(HeIIHe, THeII, THeI) * ndamp(nHeII);
            L2 = 1.0; // Langevin’s energy loss parameter
            knn = k * nHeII * nHeI;
            dEHeII += -knn * L2 * (THeII - THeI);
            dEHeI += +knn * L2 * (THeII - THeI);
            nuHeII = nuHeII + k * nHeI;
            nuHeI = nuHeI + k * nHeII;
            // H + H2
            k = RRel(HH2, TH, TH2);
            L2 = 4.0 * mi * (2.0 * mi) / pow(mi + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nH * nH2;
            dEH += -knn * L2 * (TH - TH2);
            dEH2 += +knn * L2 * (TH - TH2);
            nuH = nuH + k * nH2;
            nuH2 = nuH2 + k * nH;
            // HeI + H2 // not the correct reaction rate, but info is missing (not including this reaction gives strange temperature for HeI at start of calculation) 4/02/2015
            k = RRel(HH2, THeI, TH2);
            L2 = 4.0 * (4.0 * mi) * (2.0 * mi) / pow((4.0 * mi) + (2.0 * mi), 2.0); // Langevin’s energy loss parameter
            knn = k * nHeI * nH2;
            dEHeI += -knn * L2 * (THeI - TH2);
            dEH2 += +knn * L2 * (THeI - TH2);
            nuHeI = nuHeI + k * nH2;
            nuH2 = nuH2 + k * nHeI;
            // HeI + HeI // 6/6/2017
            k = RRel(HeIHeI, THeI, THeI);
            nuHeI = nuHeI + k * nHeI;
            // H2 + H2 // 6/6/2017
            k = RRel(H2H2, TH2, TH2);
            nuH2 = nuH2 + k * nH2;
            // H + H // 6/6/2017
            k = RRel(HH, TH, TH);
            nuH = nuH + k * nH;
        }
        // END belas

        /////////////////////
        //// Coulomb collisions (bcoulomb)
        /////////////////////
        if (bcoulomb) {
            double lambdaee, vee;
            vee = 0.0;
            double lambdaei;
            double lambda11, lambda12, lambda13, lambda14, lambda15, lambda23, lambda24, lambda25, lambda34, lambda35, lambda45;
            double veHi, veH2i, veH3i, veHeII, veHeIII, QeHi, QeH2i, QeH3i, QeHeII, QeHeIII;
            double veCII, veCIII, veCIV, veCV, QeCII, QeCIII, QeCIV, QeCV;
            double v12, v13, v14, v15, v23, v24, v25, v34, v35, v45;
            double Q12, Q13, Q14, Q15, Q23, Q24, Q25, Q34, Q35, Q45;
            double v11, v22, v33, v44, v55;
            v11 = 0.0;
            v22 = 0.0;
            v33 = 0.0;
            v44 = 0.0;
            v55 = 0.0;

            double mme = 1000.0 * me;
            double mmHi = 1000.0 * mi;
            double mmH2i = 2.0 * 1000.0 * mi;
            double mmH3i = 3.0 * 1000.0 * mi;
            double mmHe = 4.0 * 1000.0 * mi;
            double mmC = 12.0 * 1000.0 * mi;
            double mu1 = 1.0, mu2 = 2.0, mu3 = 3.0, mu4 = 4.0, mu5 = 4.0;
            double Z1 = 1.0, Z2 = 1.0, Z3 = 1.0, Z4 = 1.0, Z5 = 2.0;

            double TTe;
            TTe = 1.602e-12 * Te;
            TTe += 1.602e-12 * 0.01;
            double TTHi;
            TTHi = 1.602e-12 * THi;
            TTHi += 1.602e-12 * 0.01;
            double TTH2i;
            TTH2i = 1.602e-12 * TH2i;
            TTH2i += 1.602e-12 * 0.01;
            double TTH3i;
            TTH3i = 1.602e-12 * TH3i;
            TTH3i += 1.602e-12 * 0.01;
            double TTHeII;
            TTHeII = 1.602e-12 * THeII;
            TTHeII += 1.602e-12 * 0.01;
            double TTHeIII;
            TTHeIII = 1.602e-12 * THeIII;
            TTHeIII += 1.602e-12 * 0.01;
            // cout << "   " << im << "     Te = " << Te << ",  TTe= " << TTe/1.602e-12 <<endl;
            double qeCGS = 4.8032e-10; // esu

            if (ne >= 1.0e1) {
                lambdaee = 23.5 - log(sqrt(ne) * pow(Te, -1.25)) - pow(1e-5 + pow(log(Te) - 2.0, 2.0) / 16.0, 0.5);
                vee = 2.91e-6 * ne * lambdaee * pow(Te, -3.0 / 2.0);
                nue = nue + vee;
                if (!bselfcol) {
                    colrateRF.nue[im] = -vee;
                }
                // if (Te > 10*pow(Z1,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=THi*Z1*me/mi/mu1)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * Z1 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<THi*Z1*me/mi/mu1)
                // {   lambdaei = (30.0 - log( sqrt(nHi) * pow(THi,-3.0/2.0) * pow(Z1,2.0) / mu1 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veHi = nHi * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmHi * pow(TTe / mme + TTHi / mmHi, 1.5));
                QeHi = (Te - THi) * veHi * ne; // 3/2 *
                nue = nue + veHi;
                if (veHi < 0.0) {
                    cout << "eHi " << endl;
                }
                nuHi = nuHi + veHi;

                // if (Te > 10*pow(Z2,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=TH2i*Z2*me/mi/mu2)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * Z2 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<TH2i*Z2*me/mi/mu2)
                // {   lambdaei = (30.0 - log( sqrt(nH2i) * pow(TH2i,-3.0/2.0) * pow(Z2,2.0) / mu2 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veH2i = nH2i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmH2i * pow(TTe / mme + TTH2i / mmH2i, 1.5));
                QeH2i = (Te - TH2i) * veH2i * ne; // 3/2 *
                nue = nue + veH2i;
                nuH2i = nuH2i + veH2i;
                if (veH2i < 0.0) {
                    cout << "eH2i " << endl;
                }

                // if (Te > 10*pow(Z3,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=TH3i*Z3*me/mi/mu3)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * Z3 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<TH3i*Z3*me/mi/mu3)
                // {   lambdaei = (30.0 - log( sqrt(nH3i) * pow(TH3i,-3.0/2.0) * pow(Z3,2.0) / mu3 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veH3i = nH3i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmH3i * pow(TTe / mme + TTH3i / mmH3i, 1.5));
                QeH3i = (Te - TH3i) * veH3i * ne; // 3/2 *
                nue = nue + veH3i;
                nuH3i = nuH3i + veH3i;
                if (veH3i < 0.0) {
                    cout << "eH3i " << endl;
                }

                // if (Te > 10*pow(Z4,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=THeII*Z4*me/mi/mu4)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * Z4 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<THeII*Z4*me/mi/mu4)
                // {   lambdaei = (30.0 - log( sqrt(nHeII) * pow(THeII,-3.0/2.0) * pow(Z4,2.0) / mu4 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veHeII = nHeII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmHe * pow(TTe / mme + TTHeII / mmHe, 1.5));
                QeHeII = (Te - THeII) * veHeII * ne; // 3/2 *
                nue = nue + veHeII;
                nuHeII = nuHeII + veHeII;
                if (veHeII < 0.0) {
                    cout << "eHeII " << endl;
                }

                // if (Te > 10*pow(Z5,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=THeIII*Z5*me/mi/mu5)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * Z5 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<THeIII*Z5*me/mi/mu5)
                // {   lambdaei = (30.0 - log( sqrt(nHeIII) * pow(THeIII,-3.0/2.0) * pow(Z5,2.0) / mu5 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veHeIII = nHeIII * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmHe * pow(TTe / mme + TTHeIII / mmHe, 1.5));
                QeHeIII = (Te - THeIII) * veHeIII * ne; // 3/2 *
                nue = nue + veHeIII;
                nuHeIII = nuHeIII + veHeIII;
                if (veHeIII < 0.0) {
                    cout << "eHeIII " << endl;
                }

                // if (Te > 10*pow(1.0,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=1.0*1.0*me/mi/12.0)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * 1.0 * pow(Te, -3.0 / 2.0)));
                // 	}
                // else if (Te<1.0*1.0*me/mi/12.0)
                // {   lambdaei = (30.0 - log( sqrt(nCII) * pow(1.0,-3.0/2.0) * pow(1.0,2.0) / 12.0 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veCII = nCII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmC * pow(TTe / mme, 1.5));
                QeCII = (Te - 1.0) * veCII * ne; // 3/2 *
                nue = nue + veCII;
                if (veCII < 0.0) {
                    cout << "eCII end " << endl;
                }

                // if (Te > 10*pow(2.0,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=1.0*2.0*me/mi/12.0)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * 2.0 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<1.0*2.0*me/mi/12.0)
                // {   lambdaei = (30.0 - log( sqrt(nCIII) * pow(1.0,-3.0/2.0) * pow(2.0,2.0) / 12.0 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veCIII = nCIII * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmC * pow(TTe / mme, 1.5));
                QeCIII = (Te - 1.0) * veCIII * ne; // 3/2 *
                nue = nue + veCIII;
                if (veCIII < 0.0) {
                    cout << "eCIII end " << endl;
                }

                // if (Te > 10*pow(3.0,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=1.0*3.0*me/mi/12.0)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * 3.0 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<1.0*3.0*me/mi/12.0)
                // {   lambdaei = (30.0 - log( sqrt(nCIV) * pow(1.0,-3.0/2.0) * pow(3.0,2.0) / 12.0 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veCIV = nCIV * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmC * pow(TTe / mme, 1.5));
                QeCIV = (Te - 1.0) * veCIV * ne; // 3/2 *
                nue = nue + veCIV;
                if (veCIV < 0.0) {
                    cout << "eCIV end " << endl;
                }

                // if (Te > 10*pow(4.0,2.0))
                // {   lambdaei = (24.0 - log( sqrt(ne) * pow(Te,-1.0)));  }
                // else if (Te>=1.0*4.0*me/mi/12.0)
                // {
                lambdaei = (23.0 - log(sqrt(ne) * 4.0 * pow(Te, -3.0 / 2.0)));
                // }
                // else if (Te<1.0*4.0*me/mi/12.0)
                // {   lambdaei = (30.0 - log( sqrt(nCV) * pow(1.0,-3.0/2.0) * pow(4.0,2.0) / 12.0 ));	}
                if (lambdaei < 10.0) {
                    lambdaei = 10.0;
                } else if (lambdaei > 20.0) {
                    lambdaei = 20.0;
                }
                veCV = nCV * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambdaei / (3.0 * mme * mmC * pow(TTe / mme, 1.5));
                QeCV = (Te - 1.0) * veCV * ne; // 3/2 *
                nue = nue + veCV;
                if (veCV < 0.0) {
                    cout << "eCV end " << endl;
                }
            } else {
                QeHi = 0.0;
                QeH2i = 0.0;
                QeH3i = 0.0;
                QeHeII = 0.0;
                QeHeIII = 0.0;
                QeCII = 0.0;
                QeCIII = 0.0;
                QeCIV = 0.0;
                QeCV = 0.0;
            }
            /*cout << "QeHeII  " << QeHeII << endl;
             cout << "QeHeIII  " << QeHeIII << endl;
             cout << "QeHi  " << QeHi << endl;
             cout << "QeH2i  " << QeH2i << endl;
             cout << "QeH3i  " << QeH3i << endl;
             cout << "veCII  " << veCII << endl;
             cout << "veCIII  " << veCIII << endl;
             cout << "QeCII  " << QeCII << endl;
             cout << "QeCIII  " << QeCIII << endl;*/

            if ((nHi * nHi) != 0.0) {
                lambda11 = 23.0 - log(Z1 * Z1 * (mu1 + mu1) / (mu1 * THi + mu1 * THi) * pow(pow(Z1, 2.0) * nHi / THi + pow(Z1, 2.0) * nHi / THi, 0.5));
                v11 = nHi * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda11 / (3.0 * mmHi * mmHi * pow((TTHi / mmHi + TTHi / mmHi), 1.5));
                nuHi = nuHi + v11;
            }
            if ((nH2i * nH2i) != 0.0) {
                lambda11 = 23.0 - log(Z2 * Z2 * (mu2 + mu2) / (mu2 * TH2i + mu2 * TH2i) * pow(pow(Z2, 2.0) * nH2i / TH2i + pow(Z2, 2.0) * nH2i / TH2i, 0.5));
                v22 = nH2i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda11 / (3.0 * mmH2i * mmH2i * pow((TTH2i / mmH2i + TTH2i / mmH2i), 1.5));
                nuH2i = nuH2i + v22;
            }
            if ((nH3i * nH3i) != 0.0) {
                lambda11 = 23.0 - log(Z3 * Z3 * (mu3 + mu3) / (mu3 * TH3i + mu3 * TH3i) * pow(pow(Z3, 2.0) * nH3i / TH3i + pow(Z3, 2.0) * nH3i / TH3i, 0.5));
                v33 = nH3i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda11 / (3.0 * mmH3i * mmH3i * pow((TTH3i / mmH3i + TTH3i / mmH3i), 1.5));
                nuH3i = nuH3i + v33;
            }
            if ((nHeII * nHeII) != 0.0) {
                lambda11 = 23.0 - log(Z4 * Z4 * (mu4 + mu4) / (mu4 * THeII + mu4 * THeII) * pow(pow(Z4, 2.0) * nHeII / THeII + pow(Z4, 2.0) * nHeII / THeII, 0.5));
                v44 = nHeII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda11 / (3.0 * mmHe * mmHe * pow((TTHeII / mmHe + TTHeII / mmHe), 1.5));
                nuHeII = nuHeII + v44;
            }
            if ((nHeIII * nHeIII) != 0.0) {
                lambda11 = 23.0 - log(Z5 * Z5 * (mu5 + mu5) / (mu5 * THeIII + mu5 * THeIII) * pow(pow(Z5, 2.0) * nHeIII / THeIII + pow(Z5, 2.0) * nHeIII / THeIII, 0.5));
                v55 = nHeIII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda11 / (3.0 * mmHe * mmHe * pow((TTHeIII / mmHe + TTHeIII / mmHe), 1.5));
                nuHeIII = nuHeIII + v55;
            }
            if (!bselfcol) {
                colrateRF.nuHi[im] = -v11;
                colrateRF.nuH2i[im] = -v22;
                colrateRF.nuH3i[im] = -v33;
                colrateRF.nuHeII[im] = -v44;
                colrateRF.nuHeIII[im] = -v55;
            }
            //////
            if ((nHi * nH2i) != 0.0) {
                lambda12 = 23.0 - log(Z1 * Z2 * (mu1 + mu2) / (mu1 * TH2i + mu2 * THi) * pow(pow(Z1, 2.0) * nHi / THi + pow(Z2, 2.0) * nH2i / TH2i, 0.5));
                v12 = nH2i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda12 / (3.0 * mmHi * mmH2i * pow((TTHi / mmHi + TTH2i / mmH2i), 1.5));
                Q12 = (THi - TH2i) * v12 * nHi; // 3/2 *
                nuHi = nuHi + v12;
                nuH2i = nuH2i + v12;
            } else {
                Q12 = 0.0;
            }

            if ((nHi * nH3i) != 0.0) {
                lambda13 = 23.0 - log(Z1 * Z3 * (mu1 + mu3) / (mu1 * TH3i + mu3 * THi) * pow(pow(Z1, 2.0) * nHi / THi + pow(Z3, 2.0) * nH3i / TH3i, 0.5));
                v13 = nH3i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda13 / (3.0 * mmHi * mmH3i * pow((TTHi / mmHi + TTH3i / mmH3i), 1.5));
                Q13 = (THi - TH3i) * v13 * nHi; // 3/2 *
                nuHi = nuHi + v13;
                nuH3i = nuH3i + v13;
            } else {
                Q13 = 0.0;
            }

            if ((nHi * nHeII) != 0.0) {
                lambda14 = 23.0 - log(Z1 * Z4 * (mu1 + mu4) / (mu1 * THeII + mu4 * THi) * pow(pow(Z1, 2.0) * nHi / THi + pow(Z4, 2.0) * nHeII / THeII, 0.5));
                v14 = nHeII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda14 / (3.0 * mmHi * mmHe * pow((TTHi / mmHi + TTHeII / mmHe), 1.5));
                Q14 = (THi - THeII) * v14 * nHi; // 3/2 *
                nuHi = nuHi + v14;
                nuHeII = nuHeII + v14;
            } else {
                Q14 = 0.0;
            }

            if ((nHi * nHeIII) != 0.0) {
                lambda15 = 23.0 - log(Z1 * Z5 * (mu1 + mu5) / (mu1 * THeIII + mu5 * THi) * pow(pow(Z1, 2.0) * nHi / THi + pow(Z5, 2.0) * nHeIII / THeIII, 0.5));
                v15 = nHeIII * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambda15 / (3.0 * mmHi * mmHe * pow((TTHi / mmHi + TTHeIII / mmHe), 1.5));
                Q15 = (THi - THeIII) * v15 * nHi; // 3/2 *
                nuHi = nuHi + v15;
                nuHeIII = nuHeIII + v15;
            } else {
                Q15 = 0.0;
            }

            if ((nH2i * nH3i) != 0.0) {
                lambda23 = 23.0 - log(Z2 * Z3 * (mu2 + mu3) / (mu2 * TH3i + mu3 * TH2i) * pow(pow(Z2, 2.0) * nH2i / TH2i + pow(Z3, 2.0) * nH3i / TH3i, 0.5));
                v23 = nH3i * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda23 / (3.0 * mmH2i * mmH3i * pow((TTH2i / mmH2i + TTH3i / mmH3i), 1.5));
                Q23 = (TH2i - TH3i) * v23 * nH2i; // 3/2 *
                nuH2i = nuH2i + v23;
                nuH3i = nuH3i + v23;
            } else {
                Q23 = 0.0;
            }

            if ((nH2i * nHeII) != 0.0) {
                lambda24 = 23.0 - log(Z2 * Z4 * (mu2 + mu4) / (mu2 * THeII + mu4 * TH2i) * pow(pow(Z2, 2.0) * nH2i / TH2i + pow(Z4, 2.0) * nHeII / THeII, 0.5));
                v24 = nHeII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda24 / (3.0 * mmH2i * mmHe * pow((TTH2i / mmH2i + TTHeII / mmHe), 1.5));
                Q24 = (TH2i - THeII) * v24 * nH2i; // 3/2 *
                nuH2i = nuH2i + v24;
                nuHeII = nuHeII + v24;
            } else {
                Q24 = 0.0;
            }

            if ((nH2i * nHeIII) != 0.0) {
                lambda25 = 23.0 - log(Z2 * Z5 * (mu2 + mu5) / (mu2 * THeIII + mu5 * TH2i) * pow(pow(Z2, 2.0) * nH2i / TH2i + pow(Z5, 2.0) * nHeIII / THeIII, 0.5));
                v25 = nHeIII * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambda25 / (3.0 * mmH2i * mmHe * pow((TTH2i / mmH2i + TTHeIII / mmHe), 1.5));
                Q25 = (TH2i - THeIII) * v25 * nH2i; // 3/2 *
                nuH2i = nuH2i + v25;
                nuHeIII = nuHeIII + v25;
            } else {
                Q25 = 0.0;
            }

            if ((nH3i * nHeII) != 0.0) {
                lambda34 = 23.0 - log(Z3 * Z4 * (mu3 + mu4) / (mu3 * THeII + mu4 * TH3i) * pow(pow(Z3, 2.0) * nH3i / TH3i + pow(Z4, 2.0) * nHeII / THeII, 0.5));
                v34 = nHeII * pow(qeCGS, 4.0) * 8.0 * sqrt(2.0 * pi) * lambda34 / (3.0 * mmH3i * mmHe * pow((TTH3i / mmH3i + TTHeII / mmHe), 1.5));
                Q34 = (TH3i - THeII) * v34 * nH3i; // 3/2 *
                nuH3i = nuH3i + v34;
                nuHeII = nuHeII + v34;
            } else {
                Q34 = 0.0;
            }

            if ((nH3i * nHeIII) != 0.0) {
                lambda35 = 23.0 - log(Z3 * Z5 * (mu3 + mu5) / (mu3 * THeIII + mu5 * TH3i) * pow(pow(Z3, 2.0) * nH3i / TH3i + pow(Z5, 2.0) * nHeIII / THeIII, 0.5));
                v35 = nHeIII * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambda35 / (3.0 * mmH3i * mmHe * pow((TTH3i / mmH3i + TTHeIII / mmHe), 1.5));
                Q35 = (TH3i - THeIII) * v35 * nH3i; // 3/2 *
                nuH3i = nuH3i + v35;
                nuHeIII = nuHeIII + v35;
            } else {
                Q35 = 0.0;
            }

            if ((nHeII * nHeIII) != 0.0) {
                lambda45 = 23.0 - log(Z4 * Z5 * (mu4 + mu5) / (mu4 * THeIII + mu5 * THeII) * pow(pow(Z4, 2.0) * nHeII / THeII + pow(Z5, 2.0) * nHeIII / THeIII, 0.5));
                v45 = nHeIII * pow(qeCGS, 4.0) * 4.0 * 8.0 * sqrt(2.0 * pi) * lambda45 / (3.0 * mmHe * mmHe * pow((TTHeII / mmHe + TTHeIII / mmHe), 1.5));
                Q45 = (THeII - THeIII) * v45 * nHeII; // 3/2 *
                nuHeII = nuHeII + v45;
                nuHeIII = nuHeIII + v45;
            } else {
                Q45 = 0.0;
            }

            /*cout << "Q12  " << Q12 << endl;
             cout << "Q13  " << Q13 << endl;
             cout << "Q14  " << Q14 << endl;
             cout << "Q15  " << Q15 << endl;
             cout << "Q23  " << Q23 << endl;
             cout << "Q24  " << Q24 << endl;
             cout << "Q25  " << Q25 << endl;
             cout << "Q34  " << Q34 << endl;
             cout << "Q35  " << Q35 << endl;
             cout << "Q45  " << Q45 << endl;*/

            //////// dn and dE for coulomb collisions
            dEe += -QeHi - QeH2i - QeH3i - QeHeII - QeHeIII - QeCII - QeCIII - QeCIV - QeCV;
            dEHi += +QeHi - Q12 - Q13 - Q14 - Q15;
            dEH2i += +QeH2i + Q12 - Q23 - Q24 - Q25;
            dEH3i += +QeH3i + Q13 + Q23 - Q34 - Q35;
            dEHeII += +QeHeII + Q14 + Q24 + Q34 - Q45;
            dEHeIII += +QeHeIII + Q15 + Q25 + Q35 + Q45;
        }
        // END bcoulomb

        /////////////////////
        //// Impurities Te (bimpur)
        /////////////////////
        if (bimpur) // correct Te and reaction rates
        {
            if (nCI > 0.0) {
                double Tx = Te / 1000.0;
                static double exp1I = 0.0;
                static double exp1II = 0.0;
                static double exp2III = 0.0;
                static double exp1IV = 0.0;
                static double exp4V = 0.0;
                //     if Tec<0.1
                //         Tec=0.1;
                //     end
                knn = RRC(IC1, Te) * ne * nCI; // cout << k << "  " ;
                dnCI += -knn;
                dnCII += +knn;
                dne += +knn;
                knn = RRC(IC2, Te) * ne * nCII; // cout << k << "  " ;
                dnCII += -knn;
                dnCIII += +knn;
                dne += +knn;
                knn = RRC(IC3, Te) * ne * nCIII; // cout << k << "  " ;
                dnCIII += -knn;
                dnCIV += +knn;
                dne += +knn;
                knn = RRC(IC4, Te) * ne * nCIV; // cout << k << endl;
                dnCIV += -knn;
                dnCV += +knn;
                dne += +knn;

                // Cooling rates given in IAEA
                // I fit 1
                static double CIa1 = .6448E+03, CIa2 = .3824E-02, CIa3 = 1.1008, CIa4 = .3730, CIa5 = .2766, CIa6 = -.1836;
                // II fit 1
                static double CIIa1 = .4519E+03, CIIa2 = .6083E-01, CIIa3 = .7313, CIIa4 = .3943, CIIa5 = .4318, CIIa6 = -.4402E-01;
                // III fit 2
                static double CIIIa1 = .4368E+03, CIIIa2 = .9100E-04, CIIIa3 = .16421, CIIIa4 = .2334E-01, CIIIa5 = -1.1336, CIIIa6 = 1.1203, CIIIa7 = .4493;
                // IV fit 1
                static double CIVa1 = .1738E+03, CIVa2 = .1640E-03, CIVa3 = 1.5417, CIVa4 = .3624, CIVa5 = .4363, CIVa6 = -.2909;
                // V fit 4
                static double CVa1 = .1173E+03, CVa2 = .2612, CVa3 = .9629, CVa4 = -.4263, CVa5 = .9247, CVa6 = .3279, CVa7 = .6987E-02, CVa8 = .2506;

                exp1I = (CIa1 * exp(-CIa2 / pow(Tx, CIa3)) / (pow(Tx, CIa4) + CIa5 * pow(Tx, CIa6))) * 1.0e-33 * 1.0e6;                            // Wcm3
                exp1II = (CIIa1 * exp(-CIIa2 / pow(Tx, CIIa3)) / (pow(Tx, CIIa4) + CIIa5 * pow(Tx, CIIa6))) * 1.0e-33 * 1.0e6;                     // Wcm3
                exp2III = (CIIIa1 * exp(-CIIIa2 / pow(Tx, CIIIa3)) / (1 + CIIIa4 * pow(Tx, CIIIa5) + CIIIa6 * pow(Tx, CIIIa7))) * 1.0e-33 * 1.0e6; // Wcm3
                exp1IV = (CIVa1 * exp(-CIVa2 / pow(Tx, CIVa3)) / (pow(Tx, CIVa4) + CIVa5 * pow(Tx, CIVa6))) * 1e-33 * 1e6;                         // Wcm3
                exp4V = (CVa1 * exp(-CVa2 / pow(Tx, CVa3)) / (pow(Tx, CVa4) + CVa5 * pow(Tx, CVa6)) + CVa7 * pow(Tx, CVa8)) * 1e-33 * 1e6;         // Wcm3
                dEe += -2.0 / 3.0 * (exp1I * nCI + exp1II * nCII + exp2III * nCIII + exp1IV * nCIV + exp4V * nCV) * ne / qe;
                //- exp1IV  *(nH20+nHeI0)*Imp/100*ne/qe...
                //- exp4V   *(nH20+nHeI0)*Imp/100*ne/qe ;
                // cout << exp1I << "  " << exp1II << "  " << exp2III << "  " << exp1IV << "  " << exp4V << endl;
            }
        }
        // END bimpur

        // Update variables
        dnr.dne[im] = dne;
        dEr.dEe[im] = dEe;
        colrate.nue[im] = nue;
        dnr.dnH[im] = dnH;
        dEr.dEH[im] = dEH;
        colrate.nuH[im] = nuH;
        dnr.dnH2[im] = dnH2;
        dEr.dEH2[im] = dEH2;
        colrate.nuH2[im] = nuH2;
        dnr.dnHi[im] = dnHi;
        dEr.dEHi[im] = dEHi;
        colrate.nuHi[im] = nuHi;
        dnr.dnH2i[im] = dnH2i;
        dEr.dEH2i[im] = dEH2i;
        colrate.nuH2i[im] = nuH2i;
        dnr.dnH3i[im] = dnH3i;
        dEr.dEH3i[im] = dEH3i;
        colrate.nuH3i[im] = nuH3i;
        dnr.dnHeI[im] = dnHeI;
        dEr.dEHeI[im] = dEHeI;
        colrate.nuHeI[im] = nuHeI;
        dnr.dnHeII[im] = dnHeII;
        dEr.dEHeII[im] = dEHeII;
        colrate.nuHeII[im] = nuHeII;
        dnr.dnHeIII[im] = dnHeIII;
        dEr.dEHeIII[im] = dEHeIII;
        colrate.nuHeIII[im] = nuHeIII;
        dnr.dnCI[im] = dnCI;
        dnr.dnCII[im] = dnCII;
        dnr.dnCIII[im] = dnCIII;
        dnr.dnCIV[im] = dnCIV;
        dnr.dnCV[im] = dnCV;
        // cout << "   " << dnHeII << endl;
        // if (!bselfcol){
        colrateRF.nue[im] += colrate.nue[im]; // both are equal in case of bselfcol
        colrateRF.nuHi[im] += colrate.nuHi[im];
        colrateRF.nuH2i[im] += colrate.nuH2i[im];
        colrateRF.nuH3i[im] += colrate.nuH3i[im];
        colrateRF.nuHeII[im] += colrate.nuHeII[im];
        colrateRF.nuHeIII[im] += colrate.nuHeIII[im];
        //}
    }

    return;
}

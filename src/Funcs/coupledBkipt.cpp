#include "coupledpower.h"

void bkipt_func(double &alr) {
    double powall = 0.0;
    // part prop to ne cylindrical coordinates
    // double nelimit=1.e8;
    double lne = 0.0;
    // double nemean=0.0;
    // double eprof[NMESHP],inteprof=0.0;

    /*if (bkipt_decay)
    {
        for (int id=0; id<NMESHP; ++id)
        { // line integrated density
            eprof[id] = exp(-sqrt(pow(aR[id]-Rant,2.0)/pow(3.0*a/2.0,2.0)));
        }
        for (int id=0; id<NMESHP-1; ++id)
        { // line integrated density
            inteprof += 0.5*(eprof[id]*aR[id]+eprof[id+1]*aR[id+1])*(aR[id+1]-aR[id]);
        }
        for (int id=0; id<NMESHP; ++id)
        { // line integrated density
            eprof[id] = eprof[id]/inteprof;
        }
    }*/

    for (int id = 0; id < NMESHP - 1; ++id) { // line integrated density
        lne += 0.5 * (nr.ne[id] * aR[id] + nr.ne[id + 1] * aR[id + 1]) * (aR[id + 1] - aR[id]);
    }
    // nemean=lne/(pow(aR[NMESHP-1],2.0)-pow(aR[0],2.0));
#pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        /*if (nemean < nelimit)
        {
            if (bkipt_decay){  PRFe[id] = 6.24e18 * Prf/(2.0*pi*2.0*a)*eprof[id]; }
            else if (bkipt_propto){  PRFe[id] = 6.24e18 * Prf/Vpl*(pow(aR[NMESHP-1],2.0)-pow(aR[0],2.0)) * nr.ne[id] / lne; }
            else { PRFe[id] = 0.0; }
            PRFHi[id] = 0.0;
            PRFH2i[id] = 0.0;
            PRFH3i[id] = 0.0;
            PRFHeII[id] = 0.0;
            PRFHeIII[id] = 0.0;
        }
        else
        {
            if (bkipt_decay){  PRFe[id] = 6.24e18 * Prf/(2.0*pi*2.0*a)*eprof[id]   * nelimit/nemean; }
            else if (bkipt_propto){  PRFe[id] = 6.24e18 * Prf/Vpl*(pow(aR[NMESHP-1],2.0)-pow(aR[0],2.0)) * nr.ne[id] / lne   * nelimit/nemean; }
            else {*/
        PRFe_array[id] = 0.0; //}
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
        //}
    }

    if (1) //(nemean > nelimit)
    {
        double rpos[NPOINTS];               // mesh radial points (-a:a) [cm];
        int zsort[NSORT] = {1, 1, 1, 1, 2}; //,   1, 1, 1};     // ion charges in units of proton charge; (Hi H2i H3i    HeII HeIII    Di HDi D2i)
        int msort[NSORT] = {1, 2, 3, 4, 4}; //,   2, 3, 4};     // ion mass number;
        double density[NSORT][NPOINTS * 2]; // ion density distributions [cm-3];
        double magnf[NPOINTS];              // magnetic field radial profile [G];
        double cfreq[NSORT + 1][NPOINTS];   // collision frequencies [c-1];
        double powin = Prf;                 // input power [W];
        double powr[NSORT + 2][NPOINTS];    // RF power profile delivered to each ion sort, electrons, and full RF power including power due to penaltizing [W/cm3](output);
        double powr2[NSORT + 2][NPOINTS];   // RF power profile delivered to each ion sort, electrons, and full RF power including power due to penaltizing [W/cm3](output);
        int nsort = NSORT;                  // number of sorts of ions;
        int npoints = NPOINTS;              // number of mesh points;

#pragma omp parallel for
        for (int i = 0; i < NPOINTS; ++i) {
            rpos[i] = aR[i] - R;
            // cout << i << "/"<< npoints << "     rpos = " << rpos[i] << " cm and aR = " << aR[i] << endl;
            density[0][i * 2] = nr.nHi[i];  //* HtoHD;
            density[1][i * 2] = nr.nH2i[i]; //* HtoHD*HtoHD;         // H2i
            density[2][i * 2] = nr.nH3i[i];
            density[3][i * 2] = nr.nHeII[i];
            density[4][i * 2] = nr.nHeIII[i];
            // density[5][i*2]=nHi_RF[i]                    * (1-HtoHD);
            // density[6][i*2]=nH2i_RF[i]                   * 2.0*HtoHD*(1-HtoHD); // HDi
            // density[7][i*2]=nH2i_RF[i]                   * (1-HtoHD)*(1-HtoHD); // D2i
            density[0][i * 2 + 1] = nr.nHi[i];  //* HtoHD;
            density[1][i * 2 + 1] = nr.nH2i[i]; //* HtoHD*HtoHD;         // H2i
            density[2][i * 2 + 1] = nr.nH3i[i];
            density[3][i * 2 + 1] = nr.nHeII[i];
            density[4][i * 2 + 1] = nr.nHeIII[i];
            // density[5][i*2+1]=xnHi_RF[i]                 * (1-HtoHD);
            // density[6][i*2+1]=xnH2i_RF[i]                * 2.0*HtoHD*(1-HtoHD); // HDi
            // density[7][i*2+1]=xnH2i_RF[i]                * (1-HtoHD)*(1-HtoHD); // D2i
            cfreq[0][i] = colrateRF.nuHi[i];
            cfreq[1][i] = colrateRF.nuH2i[i]; // H2i
            cfreq[2][i] = colrateRF.nuH3i[i];
            cfreq[3][i] = colrateRF.nuHeII[i];
            cfreq[4][i] = colrateRF.nuHeIII[i];
            // cfreq[5][i]=colrateRF.nuHi[i]                    ;
            // cfreq[6][i]=colrateRF.nuH2i[i]                   ; // HDi
            // cfreq[7][i]=colrateRF.nuH2i[i]                   ; // D2i
            cfreq[5][i] = colrateRF.nue[i]; // e
            magnf[i] = Br[i] * 1.0e4;
            // PRFDi[i]=0.0;
            // PRFHDi[i]=0.0;
            // PRFD2i[i]=0.0;
        }
#pragma omp parallel for // smooth collisions to avoid too peaked absorption...
        for (int i = 1; i < NPOINTS - 1; ++i) {
            cfreq[0][i] = 0.25 * (colrateRF.nuHi[i - 1] + 2.0 * colrateRF.nuHi[i] + colrateRF.nuHi[i + 1]);
            cfreq[1][i] = 0.25 * (colrateRF.nuH2i[i - 1] + 2.0 * colrateRF.nuH2i[i] + colrateRF.nuH2i[i + 1]);
            cfreq[2][i] = 0.25 * (colrateRF.nuH3i[i - 1] + 2.0 * colrateRF.nuH3i[i] + colrateRF.nuH3i[i + 1]);
            cfreq[3][i] = 0.25 * (colrateRF.nuHeII[i - 1] + 2.0 * colrateRF.nuHeII[i] + colrateRF.nuHeII[i + 1]);
            cfreq[4][i] = 0.25 * (colrateRF.nuHeIII[i - 1] + 2.0 * colrateRF.nuHeIII[i] + colrateRF.nuHeIII[i + 1]);
            cfreq[5][i] = 0.25 * (colrateRF.nue[i - 1] + 2.0 * colrateRF.nue[i] + colrateRF.nue[i + 1]);
        }
        // cfreq[0][NPOINTS-1]=0.0                    ;
        // cfreq[1][NPOINTS-1]=0.0                    ;
        // cfreq[2][NPOINTS-1]=0.0                    ;
        // cfreq[3][NPOINTS-1]=0.0                    ;
        // cfreq[4][NPOINTS-1]=0.0                    ;
        // cfreq[5][NPOINTS-1]=0.0                    ;
#pragma omp parallel for collapse(2)
        for (int i = 0; i < NPOINTS; ++i) // This loops on the columns.
        {
            for (int l = 0; l < NSORT + 2; ++l) // This loops on the rows.
            {
                powr[l][i] = 0.0;
            }
        }
        // for (int i=0; i<NPOINTS;++i) //This loops on the columns.
        // {  cout << "   " << magnf[i] << endl; }
        int lastcall = 0; // false. if there will be more calls of the subroutine and .true. if the call is last
        /////////////// CALL_FORTRAN(rfpower)(mmin, mmax, mstep, freq, zsort, msort, nsort, rpos, npoints, density, magnf, cfreq, powin, alr, powr, lastcall);
        CALL_FORTRAN(rfpower)
        (freq, zsort, msort, nsort, rpos, npoints, density, magnf, cfreq, powin, alr, powr2, lastcall);

        double powe = 0.0;
        double powHi = 0.0;
        double powH2i = 0.0;
        double powH3i = 0.0;
        double powHeII = 0.0;
        double powHeIII = 0.0;
        // double powDi=0.0;
        // double powHDi=0.0;
        // double powD2i=0.0;
        // double powtot = 0.0;
        powall = 0.0;
// double powpatch = 0.0;
#pragma omp parallel for
        for (int id = 0; id < NPOINTS; ++id) { // check if power on a grid point is negative, if yes, zero it and remove it from total also
            for (int jd = 0; jd < NSORT + 1; ++jd) {
                if (powr2[jd][id] < 0.0) {
                    powr2[NSORT + 1][id] = powr2[NSORT + 1][id] - powr2[jd][id];
                    powr2[jd][id] = 0.0;
                }
            }
        }
#pragma omp parallel for
        for (int jd = 0; jd < NSORT + 2; ++jd) {
            powr[jd][0] = powr2[jd][0];
            powr[jd][NPOINTS - 1] = powr2[jd][NPOINTS - 1];
            for (int id = 1; id < NPOINTS - 1; ++id) {
                powr[jd][id] = 0.25 * (powr2[jd][id - 1] + 2.0 * powr2[jd][id] + powr2[jd][id + 1]);
            }
        }
        // #pragma omp parallel for
        for (int id = 0; id < NPOINTS; ++id) { // check if power on a grid point is negative, if yes, zero it and remove it from total also
            for (int jd = 0; jd < NSORT + 1; ++jd) {
                if (powr[jd][id] < 0.0) {
                    powr[NSORT + 1][id] = powr[NSORT + 1][id] - powr[jd][id];
                    powr[jd][id] = 0.0;
                }
            }
        }
        for (int jd = 0; jd < NSORT + 2; ++jd) // no power on lfs outer 2 grid points
        {
            powr[jd][NPOINTS - 1] = powr[jd][NPOINTS - 1] / 1.e10;
            powr[jd][NPOINTS - 2] = powr[jd][NPOINTS - 2] / 1.e10;
        }
        for (int id = 0; id < NMESHP - 1; ++id) { // integrate power to renormalise later to the total power
            powHi += pi * 0.5 * (powr[0][id] + powr[0][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
            powH2i += pi * 0.5 * (powr[1][id] + powr[1][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
            powH3i += pi * 0.5 * (powr[2][id] + powr[2][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
            powHeII += pi * 0.5 * (powr[3][id] + powr[3][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
            powHeIII += pi * 0.5 * (powr[4][id] + powr[4][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
            // powDi     += pi * 0.5*(powr[5][id]+powr[5][id+1])  *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*pi*300.0;
            // powHDi    += pi * 0.5*(powr[6][id]+powr[6][id+1])  *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*pi*300.0;
            // powD2i    += pi * 0.5*(powr[7][id]+powr[7][id+1])  *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*pi*300.0;
            powe += pi * 0.5 * (powr[5][id] + powr[5][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
            powall += pi * 0.5 * (powr[6][id] + powr[6][id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * pi * 300.0;
        }
        // total power to species
        // powtot = powe + powHi + powH2i + powH3i + powHeII + powHeIII; // + powDi+powHDi+powD2i;
        // total power to patch
        // powpatch = powall - powtot;

        //            /////// distribute patch power to all, keeping same profiles
        //            for (int id=0; id<NMESHP; ++id)
        //            {
        //                powr[0][id] = powr[0][id]+powr[0][id]/powHi    * powpatch * powHi   /(powe+powHi+powH2i+powH3i+powHeII+powHeIII);// +powDi+powHDi+powD2i);
        //                powr[1][id] = powr[1][id]+powr[1][id]/powH2i   * powpatch * powH2i  /(powe+powHi+powH2i+powH3i+powHeII+powHeIII);// +powDi+powHDi+powD2i);
        //                powr[2][id] = powr[2][id]+powr[2][id]/powH3i   * powpatch * powH3i  /(powe+powHi+powH2i+powH3i+powHeII+powHeIII);// +powDi+powHDi+powD2i);
        //                powr[3][id] = powr[3][id]+powr[3][id]/powHeII  * powpatch * powHeII /(powe+powHi+powH2i+powH3i+powHeII+powHeIII);// +powDi+powHDi+powD2i);
        //                powr[4][id] = powr[4][id]+powr[4][id]/powHeIII * powpatch * powHeIII/(powe+powHi+powH2i+powH3i+powHeII+powHeIII);// +powDi+powHDi+powD2i);
        //                // powr[5][id] = powr[5][id]+powr[5][id]/powDi    * powpatch * powDi   /(powe+powHi+powH2i+powH3i+powHeII+powHeIII +powDi+powHDi+powD2i);
        //                // powr[6][id] = powr[6][id]+powr[6][id]/powHDi   * powpatch * powHDi  /(powe+powHi+powH2i+powH3i+powHeII+powHeIII +powDi+powHDi+powD2i);
        //                // powr[7][id] = powr[7][id]+powr[7][id]/powD2i   * powpatch * powD2i  /(powe+powHi+powH2i+powH3i+powHeII+powHeIII +powDi+powHDi+powD2i);
        //                powr[5][id] = powr[5][id]+powr[5][id]/powe     * powpatch * powe    /(powe+powHi+powH2i+powH3i+powHeII+powHeIII);// +powDi+powHDi+powD2i);
        //            }

        /////// deposit patch power to electrons at patch location
        //            for (int id=0; id<NMESHP; ++id)
        //            {
        //                powr[8][id] = powr[9][id]-powr[0][id]-powr[1][id]-powr[2][id]-powr[3][id]-powr[4][id]-powr[5][id]-powr[6][id]-powr[7][id]; // all minus ions
        //            }

        if (1) {
            // distribute patch power evenly to all species at the patch location
            double powtotid = 0.0;
            double powpatchid = 0.0;
            // #pragma omp parallel for private(powtotid, powpatchid)
            for (int id = 0; id < NMESHP; ++id) {
                powtotid = powr[0][id] + powr[1][id] + powr[2][id] + powr[3][id] + powr[4][id] + powr[5][id]; // +powr[6][id]+powr[7][id]+powr[8][id];
                powpatchid = powr[6][id] - powtotid;
                powr[0][id] += powr[0][id] / powtotid * powpatchid;
                powr[1][id] += powr[1][id] / powtotid * powpatchid;
                powr[2][id] += powr[2][id] / powtotid * powpatchid;
                powr[3][id] += powr[3][id] / powtotid * powpatchid;
                powr[4][id] += powr[4][id] / powtotid * powpatchid;
                // powr[5][id] += powr[5][id]/powtotid * powpatchid;
                // powr[6][id] += powr[6][id]/powtotid * powpatchid;
                // powr[7][id] += powr[7][id]/powtotid * powpatchid;
                powr[5][id] += powr[5][id] / powtotid * powpatchid;
            }
        } else {
            // lhr power to electrons
            // #pragma omp parallel for
            for (int id = 0; id < NMESHP; ++id) {
                powr[5][id] += powr[6][id] - (powr[0][id] + powr[1][id] + powr[2][id] + powr[3][id] + powr[4][id] + powr[5][id]);
            }
        }

        // renormalise power such that powin = powout
        // #pragma omp parallel for
        for (int id = 0; id < NMESHP; ++id) {
            powr[0][id] = powr[0][id] / powall * powin;
            powr[1][id] = powr[1][id] / powall * powin;
            powr[2][id] = powr[2][id] / powall * powin;
            powr[3][id] = powr[3][id] / powall * powin;
            powr[4][id] = powr[4][id] / powall * powin;
            powr[5][id] = powr[5][id] / powall * powin;
            // powr[6][id] = powr[6][id]/powall    * powin;
            // powr[7][id] = powr[7][id]/powall    * powin;
            // powr[8][id] = powr[8][id]/powall    * powin;
        }

        // renormalise power such that powin = powout
#pragma omp parallel for
        for (int id = 0; id < NMESHP; ++id) {
            //                PRFHi[id]    += 6.24e18 * powr[0][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a); // (2.0*pi*300.0) / (2.0*a) because kipt module considers a different 'height'...
            //                PRFH2i[id]   += 6.24e18 * powr[1][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                PRFH3i[id]   += 6.24e18 * powr[2][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                PRFHeII[id]  += 6.24e18 * powr[3][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                PRFHeIII[id] += 6.24e18 * powr[4][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                // PRFDi[id]    += 6.24e18 * powr[5][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                // PRFHDi[id]   += 6.24e18 * powr[6][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                // PRFD2i[id]   += 6.24e18 * powr[7][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            //                PRFe[id]     += 6.24e18 * powr[8][id] * (1.0-nelimit/nemean)  * (2.0*pi*300.0) / (2.0*a);
            PRFHi_array[id] = 6.24e18 * powr[0][id] * (2.0 * pi * 300.0) / (2.0 * a); // (2.0*pi*300.0) / (2.0*a) because kipt module considers a different 'height'...
            PRFH2i_array[id] = 6.24e18 * powr[1][id] * (2.0 * pi * 300.0) / (2.0 * a);
            PRFH3i_array[id] = 6.24e18 * powr[2][id] * (2.0 * pi * 300.0) / (2.0 * a);
            PRFHeII_array[id] = 6.24e18 * powr[3][id] * (2.0 * pi * 300.0) / (2.0 * a);
            PRFHeIII_array[id] = 6.24e18 * powr[4][id] * (2.0 * pi * 300.0) / (2.0 * a);
            // PRFDi[id]    += 6.24e18 * powr[5][id] * (2.0*pi*300.0) / (2.0*a);
            // PRFHDi[id]   += 6.24e18 * powr[6][id] * (2.0*pi*300.0) / (2.0*a);
            // PRFD2i[id]   += 6.24e18 * powr[7][id] * (2.0*pi*300.0) / (2.0*a);
            PRFe_array[id] = 6.24e18 * powr[5][id] * (2.0 * pi * 300.0) / (2.0 * a);
        }

        // check
        powe = 0.0;
        powHi = 0.0;
        powH2i = 0.0;
        powH3i = 0.0;
        powHeII = 0.0;
        powHeIII = 0.0;
        // powDi=0.0;
        // powHDi=0.0;
        // powD2i=0.0;
        // powtot = 0.0;
        powall = 0.0;
        // powpatch=0.0;
        for (int id = 0; id < NMESHP - 1; ++id) { // line integrated density
            // cout <<"[ " << id << " " << powr[0][id] << " ],   " ;
            powHi += pi * 0.5 * (PRFHi_array[id] + PRFHi_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            powH2i += pi * 0.5 * (PRFH2i_array[id] + PRFH2i_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            powH3i += pi * 0.5 * (PRFH3i_array[id] + PRFH3i_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            powHeII += pi * 0.5 * (PRFHeII_array[id] + PRFHeII_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            powHeIII += pi * 0.5 * (PRFHeIII_array[id] + PRFHeIII_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            // powDi     += pi * 0.5*(powr[5][id]+powr[5][id+1])  *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*pi*300.0;
            // powHDi    += pi * 0.5*(powr[6][id]+powr[6][id+1])  *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*pi*300.0;
            // powD2i    += pi * 0.5*(powr[7][id]+powr[7][id+1])  *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*pi*300.0;
            powe += pi * 0.5 * (PRFe_array[id] + PRFe_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            powall += pi * 0.5 * (PRFHi_array[id] + PRFH2i_array[id] + PRFH3i_array[id] + PRFHeII_array[id] + PRFHeIII_array[id] + PRFe_array[id] /*+powr[6][id]  +powr[7][id]  +powr[8][id]*/) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
            powall += pi * 0.5 * (PRFHi_array[id + 1] + PRFH2i_array[id + 1] + PRFH3i_array[id + 1] + PRFHeII_array[id + 1] + PRFHeIII_array[id + 1] + PRFe_array[id + 1] /*+powr[6][id+1]+powr[7][id+1]+powr[8][id+1]*/) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0)) * 2.0 * a;
        }
        // cout << "  Power (Pj/Ptot) :          Pe    =  " << fixed << setprecision(3) << powe/powall << endl
        // << "                   m/q=1 :   PHi   =  " << fixed << setprecision(3) << powHi/powall << endl
        // << "                   m/q=2 :   PH2i  =  " << fixed << setprecision(3) << powH2i/powall  <<  "     PDi   =  " << fixed << setprecision(3) << powDi/powall  << "     PHeIII (kW) =  " << fixed << setprecision(3) << powHeIII/powall << endl
        // << "                   m/q=3 :   PH3i  =  " << fixed << setprecision(3) << powH3i/powall  <<  "     PHDi  =  " << fixed << setprecision(3) << powHDi/powall << endl
        // << "                   m/q=4 :   PD2i  =  " << fixed << setprecision(3) << powD2i/powall  <<  "     PHeII =  " << fixed << setprecision(3) << powHeII/powall << endl
        // << "                             Ptot  =  " << fixed << setprecision(3) << powall/1e3 << " kW   <--->    Ppatch    =  " << fixed << setprecision(3) << powpatch/powall << endl;
        // cout << "  Power (Pj/Ptot) :          Pe    =  " << fixed << setprecision(3) << powe/powall << "        alr      =  " << fixed << setprecision(4) << alr*2.0 << endl
        // << "                   m/q=1 :   PHi   =  " << fixed << setprecision(3) << powHi/powall << endl
        // << "                   m/q=2 :   PH2i  =  " << fixed << setprecision(3) << powH2i/powall  << "     PHeIII (kW) =  " << fixed << setprecision(3) << powHeIII/powall << endl
        // << "                   m/q=3 :   PH3i  =  " << fixed << setprecision(3) << powH3i/powall  << endl
        // << "                   m/q=4 :   PHeII =  " << fixed << setprecision(3) << powHeII/powall << endl
        // << "                             Ptot  =  " << fixed << setprecision(3) << powall/6.24e18/1e3 << " kW   <--->    of which " << fixed << setprecision(3) << powpatch/1e3 << " kW patched " << endl;
        // for (int id=0; id<NMESHP; ++id)
        // {
        //     PRFHi[id]    += PRFDi[id] ;
        //     PRFH2i[id]   += PRFHDi[id] + PRFD2i[id] ;
        // }
    }
    // powall=0.0;
    // for (int id=0; id<N; ++id)
    //{  // line integrated density
    //     powall += pi * (PRFe[id]+PRFHi[id]+PRFH2i[id]+PRFH3i[id]+PRFHeII[id]+PRFHeIII[id]) *(pow(aR[id+1],2.0)-pow(aR[id],2.0))  *2.0*a;
    // }
    //// cout << "    < check total power : " << powall/1000/6.24e18 << " kW, and mean ne : " << nemean << " > " << endl;
}
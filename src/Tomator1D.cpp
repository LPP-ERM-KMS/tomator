// Last update: 2021-10-05

// Libraries
#include <omp.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <filesystem>

// Function definitions
#include "Funcs/coupledpower.h"
#include "Funcs/extractor.h"
#include "Funcs/fileHandler.h"
#include "Funcs/functions.h"
#include "Funcs/timeStep.h"
#include "Funcs/collisions.h"
//#include "Vars/reactionrates.h"

// Global variables
#include "Vars/simparam.h"

using namespace std;

// function prototyping
void simulationLoop(double tstartloop, ofstream *outFile, int timeSteps);
string getCurrentDate();
void getSimParams(const char *json_file, int timeSteps);

// main
int main(int argc, char *argv[]) {
    #ifdef debug
	    cout << "Debug mode is enabled!" << endl;
    #endif

    if (argc < 2) {
        printf("Usage: %s <json_file_path>\n", argv[0]);
        return 1;
    }
    int timeSteps = -1; // Default value indicating full simulation

    // Check for additional command-line arguments
    for (int i = 2; i < argc; i++) {
        if (strcmp(argv[i], "-t") == 0 && i + 1 < argc) {
            // Parse the next argument as the number of timesteps
            timeSteps = atoi(argv[++i]);
        }
    }

    const char *json_file_path = argv[1];
    // Get Parameters from json file
    getSimParams(json_file_path, timeSteps); 

    double tstartloop;
    double tstartcalculation = omp_get_wtime();
    Eigen::initParallel();
    time_t now = time(0);
    char *timestamp = ctime(&now);

    printf("Start time: %s\n", timestamp);

    omp_set_dynamic(1); //dynamic thread adjustment
    {
        ofstream outFile;
        string fullfilename = sOutputfolder;

        cout << "Uploading data into \"" << fullfilename << "\"\n" << endl;
        outFile.open(fullfilename.c_str(), ios::out | ios::app); // , ios::app

        init_positions();
        solverInit();

        nHeI0 = computeN(pHe, Ta0);
        nH20 = computeN(pH2, Ta0);
        nCI0 = 0.001 * (pH2 + pHe) * 100.0 / (Ta0 * 11600.0 * kb) / 1.0e6;

        if (bfinput == false) {
            writeToOutFile(&outFile,timestamp);
            initializeSpecies(); //sets all initial densities
        } else {
            infile(sinputfile);
        }

        // nr_m1 = nr; HERE
        // Er_m1 = Er;
        nr1 = nr;
        Er1 = Er;
        tstartloop = omp_get_wtime();

        string basefolder = std::getenv("TOMATORSOURCE");
        string filename_ = basefolder + "/src/SimParams/Public/hydhel.tex";
        int n = filename_.length();
        char filename[n];
        std::strcpy(filename,filename_.c_str());

        if (timeSteps > 0) {
            initReactionDataMap(filename);
            cout << "Running simulation for " << timeSteps << " time steps." << endl;
        } else {
            initReactionDataMap(filename);
            cout << "Running full simulation." << endl;
        }

        simulationLoop(tstartloop, &outFile, timeSteps);

	cout << "Simulation loop ended." << endl;

        if (timeSteps <= 0) {
            writeLast(&outFile, tstartcalculation, omp_get_wtime());
        }
        

        outFile.close();
    }
    return 0;
}

void simulationLoop(double tstartloop, ofstream *outFile, int timeSteps) { // called in main()
    while (tmain < tmainend) {
        tstartloop = omp_get_wtime();

        #pragma omp parallel for
        for (int im = 0; im < NMESHP; ++im) {
            Tr.Te[im] = Er.Ee[im] / (ENERGY_FACTOR * nr.ne[im]);
            Tr.TH[im] = Er.EH[im] / (ENERGY_FACTOR * nr.nH[im]);
            Tr.TH2[im] = Er.EH2[im] / (ENERGY_FACTOR * nr.nH2[im]);
            Tr.THi[im] = Er.EHi[im] / (ENERGY_FACTOR * nr.nHi[im]);
            Tr.TH2i[im] = Er.EH2i[im] / (ENERGY_FACTOR * nr.nH2i[im]);
            Tr.TH3i[im] = Er.EH3i[im] / (ENERGY_FACTOR * nr.nH3i[im]);
            Tr.THeI[im] = Er.EHeI[im] / (ENERGY_FACTOR * nr.nHeI[im]);
            Tr.THeII[im] = Er.EHeII[im] / (ENERGY_FACTOR * nr.nHeII[im]);
            Tr.THeIII[im] = Er.EHeIII[im] / (ENERGY_FACTOR * nr.nHeIII[im]);
        }

        // Reset dnr, dEr and colrate
        #pragma omp parallel for
        for (int im = 0; im < NMESHP; ++im) {
            dnr.dne[im] = dnr.dnH[im] = dnr.dnH2[im] = dnr.dnHi[im] = dnr.dnH2i[im] = dnr.dnH3i[im] = 0.0;
            dnr.dnHeI[im] = dnr.dnHeII[im] = dnr.dnHeIII[im] = 0.0;
            dnr.dnCI[im] = dnr.dnCII[im] = dnr.dnCIII[im] = dnr.dnCIV[im] = dnr.dnCV[im] = 0.0;
            dnr.dxne[im] = dnr.dxnH[im] = dnr.dxnH2[im] = dnr.dxnHi[im] = dnr.dxnH2i[im] = dnr.dxnH3i[im] = 0.0;
            dnr.dxnHeI[im] = dnr.dxnHeII[im] = dnr.dxnHeIII[im] = 0.0;
            dnr.dxnCI[im] = dnr.dxnCII[im] = dnr.dxnCIII[im] = dnr.dxnCIV[im] = dnr.dxnCV[im] = 0.0;

            dEr.dEe[im] = dEr.dEH[im] = dEr.dEH2[im] = dEr.dEHi[im] = dEr.dEH2i[im] = dEr.dEH3i[im] = 0.0;
            dEr.dEHeI[im] = dEr.dEHeII[im] = dEr.dEHeIII[im] = 0.0;
            dEr.dxEe[im] = dEr.dxEH[im] = dEr.dxEH2[im] = dEr.dxEHi[im] = dEr.dxEH2i[im] = dEr.dxEH3i[im] = 0.0;
            dEr.dxEHeI[im] = dEr.dxEHeII[im] = dEr.dxEHeIII[im] = 0.0;

            colrate.nue[im] = colrate.nuH[im] = colrate.nuH2[im] = colrate.nuHi[im] = colrate.nuH2i[im] = colrate.nuH3i[im] = 0.0;
            colrate.nuHeI[im] = colrate.nuHeII[im] = colrate.nuHeIII[im] = 0.0;

            colrateRF.nue[im] = colrateRF.nuH[im] = colrateRF.nuH2[im] = colrateRF.nuHi[im] = colrateRF.nuH2i[im] = colrateRF.nuH3i[im] = 0.0;
            colrateRF.nuHeI[im] = colrateRF.nuHeII[im] = colrateRF.nuHeIII[im] = 0.0;
        }

        StepTimer[0] = omp_get_wtime();
        if (bcoll) {
            collisions(); // present in file collisions.cpp
        }

        #pragma omp parallel for
        for (int im = 0; im < NMESHP; ++im) {
            dEr.dEe[im] = ENERGY_FACTOR * dEr.dEe[im];
            dEr.dEH[im] = ENERGY_FACTOR * dEr.dEH[im];
            dEr.dEHi[im] = ENERGY_FACTOR * dEr.dEHi[im];
            dEr.dEH2[im] = ENERGY_FACTOR * dEr.dEH2[im];
            dEr.dEH2i[im] = ENERGY_FACTOR * dEr.dEH2i[im];
            dEr.dEH3i[im] = ENERGY_FACTOR * dEr.dEH3i[im];
            dEr.dEHeI[im] = ENERGY_FACTOR * dEr.dEHeI[im];
            dEr.dEHeII[im] = ENERGY_FACTOR * dEr.dEHeII[im];
            dEr.dEHeIII[im] = ENERGY_FACTOR * dEr.dEHeIII[im];
        }
        StepTimer[1] = omp_get_wtime();

        if (bedge) {
            limiters(); // present in file functions.cpp
        }

        if (vdrift) {
            vdrift_function(); // present in file functions.cpp
        }

        transpCoef(); // present in transport.cpp

        if (bpol) {
            bpol_function(); // present in file functions.cpp
        }

        StepTimer[5] = omp_get_wtime();

        double alr = 0.0;
        // pecabs=0.0;
        static int RFcounter = 0;
        double lne = 0.0;
        double lnn = 0.0;

        if (RFcounter == 0) {
            // tstart = omp_get_wtime();

            if (dtRFvar) {
                Rfdifmax = 0;
                Rfdif = 0;
                double ampfactor = 1.0;
                for (int im = 0; im < NMESHP; ++im) {
                    if ((aR[im] < lHFS) | (aR[im] > lLFS)) {
                        ampfactor = 2.0;
                    } else {
                        ampfactor = 1;
                    }
                    Rfdif = ampfactor * (nRF_save.ne[im] - nr.ne[im]) / nRF_save.ne[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nH[im] - nr.nH[im]) / nRF_save.nH[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nH2[im] - nr.nH2[im]) / nRF_save.nH2[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nHi[im] - nr.nHi[im]) / nRF_save.nHi[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nH2i[im] - nr.nH2i[im]) / nRF_save.nH2i[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nH3i[im] - nr.nH3i[im]) / nRF_save.nH3i[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nHeI[im] - nr.nHeI[im]) / nRF_save.nHeI[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nHeII[im] - nr.nHeII[im]) / nRF_save.nHeII[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    Rfdif = ampfactor * (nRF_save.nHeIII[im] - nr.nHeIII[im]) / nRF_save.nHeIII[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    // Rfdif =  	ampfactor*(nRF_save.nCI[im]				- nr.nCI[im])/nRF_save.nCI[im]; 	     	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif =  	ampfactor*(nRF_save.nCII[im]			- nr.nCII[im])/nRF_save.nCII[im];       if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif =  	ampfactor*(nRF_save.nCIII[im]			- nr.nCIII[im])/nRF_save.nCIII[im];     if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif =  	ampfactor*(nRF_save.nCIV[im]			- nr.nCIV[im])/nRF_save.nCIV[im];       if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif =  	ampfactor*(nRF_save.nCV[im]				- nr.nCV[im])/nRF_save.nCV[im];        	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    Rfdif = ampfactor * (ERF_save.Ee[im] - Er.Ee[im]) / ERF_save.Ee[im];
                    if (abs(Rfdifmax) < abs(Rfdif)) {
                        Rfdifmax = Rfdif;
                    }
                    // Rfdif = 	ampfactor*(ERF_save.EH[im]      	- Er.EH[im])/ERF_save.EH[im];        		if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EH2[im]    		- Er.EH2[im])/ERF_save.EH2[im];       	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EHi[im]     	- Er.EHi[im])/ERF_save.EHi[im];        	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EH2i[im]   	 	- Er.EH2i[im])/ERF_save.EH2i[im];      	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EH3i[im]    	- Er.EH3i[im])/ERF_save.EH3i[im];      	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EHeI[im]    	- Er.EHeI[im])/ERF_save.EHeI[im];      	if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EHeII[im]   	- Er.EHeII[im])/ERF_save.EHeII[im];     if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                    // Rfdif = 	ampfactor*(ERF_save.EHeIII[im]  	- Er.EHeIII[im])/ERF_save.EHeIII[im];   if(abs(Rfdifmax) < abs(Rfdif)) { Rfdifmax = Rfdif; }
                }
                Rfdifmax = abs(Rfdifmax);
                // cout << scientific << setprecision(5) << "Rfdifmax = " << Rfdifmax << endl;
                if (Rfdifmax < 1) {
                    dtRF *= min(1 / Rfdifmax, 1.1);
                } // cout << scientific << setprecision(5) << "dtRF increased to :" << dtRF << endl;}
                else {
                    dtRF *= min(1 / Rfdifmax, 0.9);
                } // cout << scientific << setprecision(5) << "dtRF decreased to :" << dtRF << endl;}
                if (minstopcrit < 2000 && dtRFconv) {
                    if (dtRFnunc == 0) {
                        dtRFnunc = 0.8 * dtRF;
                    }
                    dtRFtemp = pow(10, 10.0 / 9.0 * log10(dtRFmin) - log10(dtRFnunc) / 9.0 - minstopcrit / 100 * (log10(dtRFmin) / 18.0 - log10(dtRFnunc) / 18.0));
                    if (dtRFtemp < dtRF) {
                        dtRF = dtRFtemp;
                    }
                }

                if (dtRF > dtRFmax) {
                    dtRF = dtRFmax;
                }
                if (dtRF < dtRFmin) {
                    dtRF = dtRFmin;
                }
            }
            // END dtRFvar

            //  started coupled power" << endl;
            // coupledpower (TeRF,  nr.ne,  nr.nH2,  nr.nHi,  nr.nH2i,  nr.nH3i,  nr.nHeI,  nr.nHeII,  nr.nHeIII,
            //               nr.xne, nr.xnH2, nr.xnHi, nr.xnH2i, nr.xnH3i, nr.xnHeI, nr.xnHeII, nr.xnHeIII,
            //               colrateRF.nue,  colrateRF.nuHi,  colrateRF.nuH2i,  colrateRF.nuH3i,  colrateRF.nuHeII,  colrateRF.nuHeIII,
            // 							alr, pecabs,
            //               aR, Br, freq, Prf, a, R, Vpl, HtoHD,
            // 							PIerrorP, mytid, dtnew,
            //               PRFe_array, PRFHi_array, PRFH2i_array, PRFH3i_array, PRFHeII_array, PRFHeIII_array);

            coupledpower(freq,alr);

            // if (bkipt == true) {cout << " Wall clock time " << (omp_get_wtime()-tstart) << " s " << endl;}
            if (dtRFvar) {
                nRF_save = nr; //////////////////////////
                ERF_save = Er;
            }

            #pragma omp parallel for
            for (int im = 0; im < NMESHP; ++im) {
                if ((bkipt == true) && (bantlr == true)) // use antenna resistance
                {
                    PRFe_array_stat[im] = PRFe_array[im] * 2.0 * alr / (2.0 * alr + avlr);
                    PRFHi_array_stat[im] = PRFHi_array[im] * 2.0 * alr / (2.0 * alr + avlr);
                    PRFH2i_array_stat[im] = PRFH2i_array[im] * 2.0 * alr / (2.0 * alr + avlr);
                    PRFH3i_array_stat[im] = PRFH3i_array[im] * 2.0 * alr / (2.0 * alr + avlr);
                    PRFHeII_array_stat[im] = PRFHeII_array[im] * 2.0 * alr / (2.0 * alr + avlr);
                    PRFHeIII_array_stat[im] = PRFHeIII_array[im] * 2.0 * alr / (2.0 * alr + avlr);
                } else {
                    PRFe_array_stat[im] = PRFe_array[im];
                    PRFHi_array_stat[im] = PRFHi_array[im];
                    PRFH2i_array_stat[im] = PRFH2i_array[im];
                    PRFH3i_array_stat[im] = PRFH3i_array[im];
                    PRFHeII_array_stat[im] = PRFHeII_array[im];
                    PRFHeIII_array_stat[im] = PRFHeIII_array[im];
                }
            }
        }

        if (bupdateRFstep == true) { // update RF every N steps
            if (RFcounter < updateRF) {
                ++RFcounter;
            }
            else {
                RFcounter = 0;
            }
        }
        else { // update RF every dt
            if (tmain < tRF + dtRF) {
                ++RFcounter;
            }
            else {
                RFcounter = 0;
                tRF = tmain;
            }
        }

        /////////////////////////////
        if ((bkipt == true) && (bantlr == false)) { // check this condition !
            ///////////////////////////// // alphalaw
            lne = 0.0;
            lnn = 0.0;
            for (int id = 0; id < NMESHP - 1; ++id) {
                lne += 0.5 * (nr.ne[id] * aR[id] + nr.ne[id + 1] * aR[id + 1]) * (aR[id + 1] - aR[id]);
                lnn += 0.5 * ((nr.nH[id] + nr.nH2[id] + nr.nHeI[id] + nr.ne[id]) * aR[id] + (nr.nH[id + 1] + nr.nH2[id + 1] + nr.nHeI[id + 1] + nr.ne[id + 1]) * aR[id + 1]) * (aR[id + 1] - aR[id]);
                // mTe += 0.5*(Tr.Te[id]*aR[id]+Tr.Te[id+1]*aR[id+1])*(aR[id+1]-aR[id])/(pow(aR[NMESHP-1],2.0)-pow(aR[0],2.0)); // some average Te...
            }
            alphaval = lne / lnn;

            #pragma omp parallel for
            for (int im = 0; im < NMESHP; ++im) {
                // Shouldn't lne and lnn be set to 0 for every MESHPoint??, indeed!
                PRFe_id[im] = PRFe_array_stat[im] * (1.0 - exp(-10.0 * alphaval / alphalaw));
                PRFHi_id[im] = PRFHi_array_stat[im] * (1.0 - exp(-10.0 * alphaval / alphalaw));
                PRFH2i_id[im] = PRFH2i_array_stat[im] * (1.0 - exp(-10.0 * alphaval / alphalaw));
                PRFH3i_id[im] = PRFH3i_array_stat[im] * (1.0 - exp(-10.0 * alphaval / alphalaw));
                PRFHeII_id[im] = PRFHeII_array_stat[im] * (1.0 - exp(-10.0 * alphaval / alphalaw));
                PRFHeIII_id[im] = PRFHeIII_array_stat[im] * (1.0 - exp(-10.0 * alphaval / alphalaw));
            }
        }
        else {
            #pragma omp parallel for
            for (int im = 0; im < NMESHP; ++im) {
                PRFe_id[im] = PRFe_array_stat[im];
                PRFHi_id[im] = PRFHi_array_stat[im];
                PRFH2i_id[im] = PRFH2i_array_stat[im];
                PRFH3i_id[im] = PRFH3i_array_stat[im];
                PRFHeII_id[im] = PRFHeII_array_stat[im];
                PRFHeIII_id[im] = PRFHeIII_array_stat[im];
            }
        }

        #pragma omp parallel for
        for (int im = 0; im < NMESHP; ++im) {
            dEr.dEe[im] += PRFe_id[im];
            dEr.dEHi[im] += PRFHi_id[im];
            dEr.dEH2i[im] += PRFH2i_id[im];
            dEr.dEH3i[im] += PRFH3i_id[im];
            dEr.dEHeII[im] += PRFHeII_id[im];
            dEr.dEHeIII[im] += PRFHeIII_id[im];
        }

        StepTimer[6] = omp_get_wtime();

        if (((cnt_save == 0) & !(isnan(nr.ne[0]))) || (timeSteps > 0 && Nit >= timeSteps)) {
            writePhysicalStates(outFile);
            if (timeSteps > 0 && Nit >= timeSteps) {
                return;
            }
        }

        StepTimer[7] = omp_get_wtime();

        nr1bis = nr;
        Er1bis = Er;
        timeStep(); // solve a first time with the previously determined time step, then improve if needed
        StepTimer[8] = omp_get_wtime();

        if (Nit % Nlog == 0 && Nit != 0) {
            cout << scientific << setprecision(2)
                 << "N = " << Nit
                 << ",   t = " << tmain
                 << ",   last dt = " << dtnew
                 << ",   ne[cc] = " << nr.ne[cc]
                 << ",   Te[cc] = " << Tr.Te[cc];
            if (bgray | bram | bnefix ){
            cout << ",   Pecabs = " << pecabs;
            }
            cout << ",   D[cc] = " << Dion[cc]
                 << ",   V[cc] = " << Vion[cc];
                 // << ",   loop time = " << (omp_get_wtime() - tstartloop) << " s"
                 // << ",   Av dt = " << AvTimeStep
            if (btunedv){
            cout << ",   fD = " << Dfsave
                 << ",   fV = " << Vfsave ;
            }
            // if (bnefix){
            // cout << ",   PI_P = " << PIerrorP ;
            // }
            cout << endl;

        }

        nr1 = nr1bis;
        Er1 = Er1bis;

        if (btendvar || baccurvar) {

            tsum += oldtstep;
            if (cnt_save == 0) {

                tdifmax = 0.0;
                tdifne = 0.0;
                tdifTe = 0.0;
                for (int im = 0; im < NMESHP; ++im) {
                    tdifne += pow(abs(((n_save.ne[im] - nr.ne[im]) / n_save.ne[im]) / tsum), 2); //	if(abs(tdifmax) < abs(tdif)) { tdifmax = tdif; }
                    tdifTe += pow(abs(((E_save.Ee[im] - Er.Ee[im]) / E_save.Ee[im]) / tsum), 2); //	if(abs(tdifmax) < abs(tdif)) { tdifmax = tdif; }
                }
                tdifne = sqrt(tdifne);
                tdifTe = sqrt(tdifTe);
                tdifmax = max(tdifne, tdifTe);

                convsave = (int)(10 / minstopcrit / dtsave);
                if (convsave > convsavemax) {
                    convsave = convsavemax;
                }
                else if (convsave < convsavemin) {
                    convsave = convsavemin;
                }

                ++convsavecount;
                if (convsavecount < 4) {
                    minstopcrit = tdifmax;
                }
                else if (convsavecount < convsave) {
                    minstopcrit = (convsavecount * minstopcrit + tdifmax) / (convsavecount + 1);
                }
                else {
                    minstopcrit = (convsave * minstopcrit + tdifmax) / (convsave + 1);
                }

                if (baccurvar) {

                    if (minstopcrit < (accurcrit * accur) && accur > minaccur) {
                        accur *= 0.9;
                        cout << "accur diminished, accur = " << accur << endl;
                    }
                    // Certain treshold of convergence exceeded, increasing the accuracy (diminishing accur) becomes interesting.
                    if (accur < minaccur) {
                        accur = minaccur;
                    }
                }

                if (btendvar && minstopcrit < convcrit) { // Convergence criterium reached
                    cout << "Solution Converged, simulation ended prematurely" << endl;
                    break;
                }

                n_save = nr;
                E_save = Er;
                tsum = 0.0;
            }
            // if ((Nit%10000)==0){
            // 	cout << "Epsilon = " << minstopcrit << endl;
            // }
        }
        // Output to screen

        if ((cnt_save == 0) & !(isnan(nr.ne[0]))) {
            // if (bOutdt) {
            //     cout << scientific << setprecision(5)
            //          << "N = " << Nit
            //          << ",   t = " << tmain
            //          << ",   ne = " << nr.ne[cc]
            //          << ",   Te = " << Tr.Te[cc]
            //          << ",   loop time = " << (omp_get_wtime() - tstartloop) << " s"
            //          << ",   Av dt = " << AvTimeStep
            //          << ",   last dt = " << dtnew << endl;
            // } else {
            //     cout << scientific << setprecision(5)
            //          << "aR = " << aR[cc]
            //          << "  <dt> = " << dtout
            //          << ",   t = " << tmain
            //          << ",   ne = " << nr.ne[cc]
            //          << ",   Te = " << Tr.Te[cc]
            //          << ",   loop time = " << (omp_get_wtime() - tstartloop) << " s" << endl;
            // }
            cnt_save = -Nloopsave + 1;

            // #pragma omp barrier
            dtout = 0;
        }
        else {
            ++cnt_save;
            if (bOutdt) {
                cnt_save = -Nloopsave;
            }
        }
        if (bOutdt) {
            if (tmain > (double)cnt_loop * dtsave) {
                ++cnt_loop;
                cnt_save = 0;
            }
        }

        AverageTimer = (Nit * AverageTimer + (omp_get_wtime() - tstartloop)) / (Nit + 1);
        ColTimer = (Nit * ColTimer + (StepTimer[1] - StepTimer[0])) / (Nit + 1);
        ParTimer = (Nit * ParTimer + (StepTimer[5] - StepTimer[1])) / (Nit + 1);
        RFTimer = (Nit * RFTimer + (StepTimer[6] - StepTimer[5])) / (Nit + 1);
        TimeStepTimer = (Nit * TimeStepTimer + (StepTimer[8] - StepTimer[7])) / (Nit + 1);

        ++Nit;
    }
}

string getCurrentDate() {
    std::time_t t = std::time(nullptr);
    char buffer[100];
    if (std::strftime(buffer, sizeof(buffer), "Res_%Y%m%d_%H%M%S.csv", std::localtime(&t))) {
        return std::string(buffer);
    } else {
        return "Res_UnknownDateTime.csv";
    }
}

void getSimParams(const char *json_file, int timeSteps) {

    char *file_content = read_file(json_file);
    if (!file_content) {
        cout << "\033[1;31m";
        cout << "[Error]: " << json_file << " does not exist or is not readable" << endl;
        cout << "\033[0m";
        exit(EXIT_FAILURE);
    }

    extract_magnetic_field(file_content);
    extract_toroidal_machine_geometry(file_content);
    extract_neutral_pressure(file_content);
    extract_rf_power(file_content);
    extract_type(file_content);
    extract_general_ec(file_content);
    extract_necfix(file_content);
    extract_tomas(file_content);
    extract_general_ic(file_content);
    extract_other(file_content);
    extract_diffusion(file_content);
    extract_convection(file_content);
    extract_tune_transport(file_content);
    extract_physics_to_include(file_content);
    extract_initial_conditions(file_content);
    extract_edge_conditions(file_content);
    extract_simulation_grid(file_content);
    extract_input_file(file_content);
    extract_time_step(file_content);
    extract_time_step_for_rf_coupling(file_content);
    extract_advanced_time_step_settings(file_content);
    extract_output_parameters(file_content);
    extract_solver_parameters(file_content);

    if (NMESHP != nmeshp) { 
        // Error, print red message and exit
        cout << "\033[1;31m";
        cout << "[Error]: To change the number of meshpoints, modify line 6 in simparam.h and recompile" << endl;
        cout << "\033[0m";
        exit(EXIT_FAILURE);
    }

    if (bfinput) { // Get file name
        cout << "Input file: " << sinputfile << endl;
        extract_input_file_name(file_content);
        free(file_content);
        if (sinputfile.length() > 2) { // Check if the string has more than two characters
            sOutputfolder = sinputfile.substr(0, sinputfile.length());
        } else {
            // Handle the case where the string is too short to remove characters
            sOutputfolder = sinputfile; // Or however you want to handle this case
        }
        return;
        
    }

    if (bmanuel) {
        extract_bmanuel_input(file_content);
    }

    free(file_content);
    
    
    if (timeSteps > 0) {
        // Only used during test routine
        std::filesystem::path testFolderPath = "Data/Test";

        if (!std::filesystem::exists(testFolderPath)) {
            std::filesystem::create_directories(testFolderPath);
        }

        std::filesystem::path outputFilePath = testFolderPath / "output.csv";
        sOutputfolder = outputFilePath.string();

        std::ofstream outFile(outputFilePath, std::ios::out | std::ios::trunc);
        if (!outFile) {
            std::cerr << "Cannot open output file: " << outputFilePath << std::endl;
            return;
        }

        outFile.close();
        cout << "output: " << outputFilePath << endl;

        return;
    }

    std::string folderName = json_file;
    const std::string extension = ".json";
    size_t extPos = folderName.rfind(extension);
    if (extPos != std::string::npos) {
        folderName.erase(extPos, extension.length());
        folderName.erase(0, folderName.find_last_of("/\\") + 1);
    }

    // Use std::filesystem::path for safer path operations
    string basefolder = std::getenv("TOMATORRESULTS");
    std::filesystem::path basePath(basefolder);
    std::filesystem::path dataFolderPath = basePath / folderName;

    // Check if the directory exists and if not, create it
    if (!std::filesystem::exists(dataFolderPath)) {
        std::filesystem::create_directories(dataFolderPath);
    }

    // Update sOutputfolder (assuming it's a global or member variable)
    std::string date = getCurrentDate();
    std::filesystem::path outputFolderPath = dataFolderPath / date;
    sOutputfolder = basefolder + "/" + folderName + "/" + date;

    // TODO: copy input json file into the output folder

}

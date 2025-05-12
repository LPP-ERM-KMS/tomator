#include "fileHandler.h"

void infile(string name) {

    // Construct the input filename
    std::ifstream myInputFile(name);
    cout << "Reading from file: " << name << endl;
    if (!myInputFile.is_open()) {
        std::cerr << "Failed to open file: " << name << std::endl;
        return;
    }

    std::string line;
    std::vector<std::vector<double>> allData;
    bool isFirstLine = true;                 // Flag to determine if we are at the first line
    std::map<std::string, size_t> headerMap; // Map to associate headers with column indices
    size_t columnIndex = 0;

    while (getline(myInputFile, line)) {
        std::istringstream iss(line);
        std::string value;
        std::vector<double> lineData;

        if (isFirstLine) {
            // Process the header line
            while (getline(iss, value, ',')) {
                headerMap[value] = columnIndex++;
            }
            isFirstLine = false; // Reset the flag after reading the header
            continue;            // Skip the rest of the loop to avoid processing headers as data
        }


        columnIndex = 0; // Reset columnIndex for each new line of data
        while (getline(iss, value, ',')) {
            try {
                // Convert value to double and add to lineData
                lineData.push_back(std::stod(value));
            } catch (const std::invalid_argument &e) {
                cerr << "Invalid argument for stod, value: " << value << " at column index " << columnIndex << endl;
                // Handle the error, possibly by skipping the line or setting a default value
            }
            columnIndex++;
        }

        // Check if this line's tmain is the latest one
        if (!lineData.empty() && lineData[0] > tmain) {
            tmain = lineData[0];
            allData.clear(); // Clear previous data as we found a new latest time
        }

        if (lineData[0] == tmain) {
            allData.push_back(lineData); // Add data for the current latest tmain
        }
    }

    // Save the data to the global variables
    for (long unsigned int i = 0; i < allData.size(); i++) {
        int im = i;
        nr.ne[im] = allData[i][headerMap["ne"]];
        Er.Ee[im] = allData[i][headerMap["Ee"]];
        nr.xne[im] = allData[i][headerMap["xne"]];
        Er.xEe[im] = allData[i][headerMap["xEe"]];
        nr.nH[im] = allData[i][headerMap["nH"]];
        Er.EH[im] = allData[i][headerMap["EH"]];
        nr.xnH[im] = allData[i][headerMap["xnH"]];
        Er.xEH[im] = allData[i][headerMap["xEH"]];
        nr.nH2[im] = allData[i][headerMap["nH2"]];
        Er.EH2[im] = allData[i][headerMap["EH2"]];
        nr.xnH2[im] = allData[i][headerMap["xnH2"]];
        Er.xEH2[im] = allData[i][headerMap["xEH2"]];
        nr.nHi[im] = allData[i][headerMap["nHi"]];
        Er.EHi[im] = allData[i][headerMap["EHi"]];
        nr.xnHi[im] = allData[i][headerMap["xnHi"]];
        Er.xEHi[im] = allData[i][headerMap["xEHi"]];
        nr.nH2i[im] = allData[i][headerMap["nH2i"]];
        Er.EH2i[im] = allData[i][headerMap["EH2i"]];
        nr.xnH2i[im] = allData[i][headerMap["xnH2i"]];
        Er.xEH2i[im] = allData[i][headerMap["xEH2i"]];
        nr.nH3i[im] = allData[i][headerMap["nH3i"]];
        Er.EH3i[im] = allData[i][headerMap["EH3i"]];
        nr.xnH3i[im] = allData[i][headerMap["xnH3i"]];
        Er.xEH3i[im] = allData[i][headerMap["xEH3i"]];
        nr.nHeI[im] = allData[i][headerMap["nHeI"]];
        Er.EHeI[im] = allData[i][headerMap["EHeI"]];
        nr.xnHeI[im] = allData[i][headerMap["xnHeI"]];
        Er.xEHeI[im] = allData[i][headerMap["xEHeI"]];
        nr.nHeII[im] = allData[i][headerMap["nHeII"]];
        Er.EHeII[im] = allData[i][headerMap["EHeII"]];
        nr.xnHeII[im] = allData[i][headerMap["xnHeII"]];
        Er.xEHeII[im] = allData[i][headerMap["xEHeII"]];
        nr.nHeIII[im] = allData[i][headerMap["nHeIII"]];
        Er.EHeIII[im] = allData[i][headerMap["EHeIII"]];
        nr.xnHeIII[im] = allData[i][headerMap["xnHeIII"]];
        Er.xEHeIII[im] = allData[i][headerMap["xEHeIII"]];
        nr.nCI[im] = allData[i][headerMap["nCI"]];
        nr.xnCI[im] = allData[i][headerMap["xnCI"]];
        nr.nCII[im] = allData[i][headerMap["nCII"]];
        nr.xnCII[im] = allData[i][headerMap["xnCII"]];
        nr.nCIII[im] = allData[i][headerMap["nCIII"]];
        nr.xnCIII[im] = allData[i][headerMap["xnCIII"]];
        nr.nCIV[im] = allData[i][headerMap["nCIV"]];
        nr.xnCIV[im] = allData[i][headerMap["xnCIV"]];
        nr.nCV[im] = allData[i][headerMap["nCV"]];
        nr.xnCV[im] = allData[i][headerMap["xnCV"]];
    }

    tmain = allData[0][headerMap["tmain"]];

    dtnew = dtinit;

    myInputFile.close();
    cnt_loop = (int)ceil(tmain / dtsave) + 1;

    return;
}

void writeToOutFile(ofstream *outFile) {

    *outFile << "tmain,RadialPositions,"
             << "ne,Ee,dne,dEe,nue,"
             << "nH,EH,dnH,dEH,nuH,"
             << "nH2,EH2,dnH2,dEH2,nuH2,"
             << "nHi,EHi,dnHi,dEHi,nuHi,"
             << "nH2i,EH2i,dnH2i,dEH2i,nuH2i,"
             << "nH3i,EH3i,dnH3i,dEH3i,nuH3i,"
             << "nHeI,EHeI,dnHeI,dEHeI,nuHeI,"
             << "nHeII,EHeII,dnHeII,dEHeII,nuHeII,"
             << "nHeIII,EHeIII,dnHeIII,dEHeIII,nuHeIII,"
             << "nCI,xnCI,"
             << "nCII,xnCII,"
             << "nCIII,xnCIII,"
             << "nCIV,xnCIV,"
             << "nCV,xnCV,"
             << "PRFe,PRFHi,PRFH2i,PRFH3i,PRFHeII,PRFHeIII,"
             << "tnew" << endl;
             // << "WallTime" << endl;
    return;
}

void writeLast(ofstream *outFile, double tstartcalculation, double currentTime) {
    cout << scientific << setprecision(5) << "Average loop time = " << AverageTimer << "\n"
         << "Average collision time = " << ColTimer << "\n"
         << "Average parallel time = " << ParTimer << "\n"
         << "Average RF time = " << RFTimer << "\n"
         << "Average time step time = " << TimeStepTimer << endl;
    cout << scientific << setprecision(5) << "Average Time Step = " << AvTimeStep << endl;
    double calctime = currentTime - tstartcalculation;
    cout << "time to reach " << tmain << " was " << calctime << endl;

    ////////////////////////////////////////////////////
    //// Write last result
    ////////////////////////////////////////////////////
    cnt_save = 0;
    if ((cnt_save == 0) & !(isnan(nr.ne[0]))) {
        for (int im = 0; im < NMESHP; ++im) {
            outFile->fill(' ');
            outFile->width(3);
            *outFile << scientific << setprecision(6)
                     << tmain << "," << aR[im]
                     << "," << nr.ne[im] << "," << Er.Ee[im] << "," << nr.xne[im] << "," << Er.xEe[im] << "," << colrateRF.nue[im]
                     << "," << nr.nH[im] << "," << Er.EH[im] << "," << nr.xnH[im] << "," << Er.xEH[im] << "," << colrate.nuH[im]
                     << "," << nr.nH2[im] << "," << Er.EH2[im] << "," << nr.xnH2[im] << "," << Er.xEH2[im] << "," << colrate.nuH2[im]
                     << "," << nr.nHi[im] << "," << Er.EHi[im] << "," << nr.xnHi[im] << "," << Er.xEHi[im] << "," << colrateRF.nuHi[im]
                     << "," << nr.nH2i[im] << "," << Er.EH2i[im] << "," << nr.xnH2i[im] << "," << Er.xEH2i[im] << "," << colrateRF.nuH2i[im]
                     << "," << nr.nH3i[im] << "," << Er.EH3i[im] << "," << nr.xnH3i[im] << "," << Er.xEH3i[im] << "," << colrateRF.nuH3i[im]
                     << "," << nr.nHeI[im] << "," << Er.EHeI[im] << "," << nr.xnHeI[im] << "," << Er.xEHeI[im] << "," << colrate.nuHeI[im]
                     << "," << nr.nHeII[im] << "," << Er.EHeII[im] << "," << nr.xnHeII[im] << "," << Er.xEHeII[im] << "," << colrateRF.nuHeII[im]
                     << "," << nr.nHeIII[im] << "," << Er.EHeIII[im] << "," << nr.xnHeIII[im] << "," << Er.xEHeIII[im] << "," << colrateRF.nuHeIII[im]
                     << "," << nr.nCI[im] << "," << nr.xnCI[im]
                     << "," << nr.nCII[im] << "," << nr.xnCII[im]
                     << "," << nr.nCIII[im] << "," << nr.xnCIII[im]
                     << "," << nr.nCIV[im] << "," << nr.xnCIV[im]
                     << "," << nr.nCV[im] << "," << nr.xnCV[im]
                     << "," << PRFe_id[im]
                     << "," << PRFHi_id[im]
                     << "," << PRFH2i_id[im]
                     << "," << PRFH3i_id[im]
                     << "," << PRFHeII_id[im]
                     << "," << PRFHeIII_id[im];
            // if (bnefix) {
            //     *outFile << scientific << setprecision(6)
            //              << "       " << PIerrorV
            //              << "    " << PIerrorD
            //              << "    " << PIerrorP
            //              << "       " << Vfsave
            //              << "    " << Dfsave
            //              << "    " << pecabs;
            // }
            // if (bkipt) {
            //     *outFile << scientific << setprecision(6)
            //              << "    " << PRFe_id[im] << " " << PRFHi_id[im] << " " << PRFH2i_id[im] << " " << PRFH3i_id[im] << " " << PRFHeII_id[im] << " " << PRFHeIII_id[im];
            // }
            *outFile << scientific << setprecision(6)
                     // << "    " << dtnew << "    " << omp_get_wtime()-tstartcalculation << "    " << minstopcrit << "    " << dtRF << endl;
                     << "," << dtnew << endl;
        }
    }
}

void writePhysicalStates(ofstream *outFile) {
    for (int im = 0; im < NMESHP; ++im) {
        outFile->fill(' ');
        outFile->width(3);
        *outFile << scientific << setprecision(6)
                 << tmain << "," << aR[im]
                 << "," << nr.ne[im] << "," << Er.Ee[im] << "," << nr.xne[im] << "," << Er.xEe[im] << "," << colrateRF.nue[im]
                 << "," << nr.nH[im] << "," << Er.EH[im] << "," << nr.xnH[im] << "," << Er.xEH[im] << "," << colrate.nuH[im]
                 << "," << nr.nH2[im] << "," << Er.EH2[im] << "," << nr.xnH2[im] << "," << Er.xEH2[im] << "," << colrate.nuH2[im]
                 << "," << nr.nHi[im] << "," << Er.EHi[im] << "," << nr.xnHi[im] << "," << Er.xEHi[im] << "," << colrateRF.nuHi[im]
                 << "," << nr.nH2i[im] << "," << Er.EH2i[im] << "," << nr.xnH2i[im] << "," << Er.xEH2i[im] << "," << colrateRF.nuH2i[im]
                 << "," << nr.nH3i[im] << "," << Er.EH3i[im] << "," << nr.xnH3i[im] << "," << Er.xEH3i[im] << "," << colrateRF.nuH3i[im]
                 << "," << nr.nHeI[im] << "," << Er.EHeI[im] << "," << nr.xnHeI[im] << "," << Er.xEHeI[im] << "," << colrate.nuHeI[im]
                 << "," << nr.nHeII[im] << "," << Er.EHeII[im] << "," << nr.xnHeII[im] << "," << Er.xEHeII[im] << "," << colrateRF.nuHeII[im]
                 << "," << nr.nHeIII[im] << "," << Er.EHeIII[im] << "," << nr.xnHeIII[im] << "," << Er.xEHeIII[im] << "," << colrateRF.nuHeIII[im]
                 << "," << nr.nCI[im] << "," << nr.xnCI[im]
                 << "," << nr.nCII[im] << "," << nr.xnCII[im]
                 << "," << nr.nCIII[im] << "," << nr.xnCIII[im]
                 << "," << nr.nCIV[im] << "," << nr.xnCIV[im]
                 << "," << nr.nCV[im] << "," << nr.xnCV[im]
                 << "," << PRFe_id[im]
                 << "," << PRFHi_id[im]
                 << "," << PRFH2i_id[im]
                 << "," << PRFH3i_id[im]
                 << "," << PRFHeII_id[im]
                 << "," << PRFHeIII_id[im];
        // if (bnefix) {
        //     *outFile << scientific << setprecision(6)
        //              << "       " << PIerrorV
        //              << "    " << PIerrorD
        //              << "    " << PIerrorP
        //              << "       " << Vfsave
        //              << "    " << Dfsave
        //              << "    " << pecabs;
        // }
        // if (bkipt) {
        //     *outFile << scientific << setprecision(6)
        //              << "    " << PRFe_id[im] << " " << PRFHi_id[im] << " " << PRFH2i_id[im] << " " << PRFH3i_id[im] << " " << PRFHeII_id[im] << " " << PRFHeIII_id[im];
        // }
        *outFile << scientific << setprecision(6)
                 // << "    " << dtnew << "    " << omp_get_wtime()-tstartcalculation << "    " << minstopcrit << "    " << dtRF << endl;
                 << "," << dtnew << endl;
    }
}

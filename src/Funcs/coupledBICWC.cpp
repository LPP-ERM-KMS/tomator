#include "coupledpower.h"

double calculateTemperature(double energy, double density) {
    return energy / (1.5 * density);
}

void generateInput() {
    std::ofstream myOutputFile("Tomator_exe/data/TOMATOR/profiles_TOMATOR_WEST.txt");
    if (!myOutputFile.is_open()) {
        cout << "\033[1;31m"
             << "[Error]: Could not open file: "
             << "Qs.csv"
             << "\033[0m" << endl;
        exit(EXIT_FAILURE);
    }
    // Values to write: Rpos, n_e, n_Hi, n_H2i, n_H3i, T_e, T_Hi, T_H2i, T_H3i, nu_e, nu_Hi, nu_H2i, nu_H3i
    for (int im = 0; im < NMESHP; ++im) { // Assuming nH3i is always 0
        myOutputFile << aR[im] << "," << nr.ne[im] << "," << nr.nHi[im] << "," << nr.nH2i[im] << "," << nr.nH3i[im] << ","
                     << calculateTemperature(Er.Ee[im], nr.ne[im]) << "," << calculateTemperature(Er.EHi[im], nr.nHi[im]) << ","
                     << calculateTemperature(Er.EH2i[im], nr.nH2i[im]) << "," << calculateTemperature(Er.EH3i[im], nr.nH3i[im]) << ","
                     << colrateRF.nue[im] << "," << colrateRF.nuHi[im] << "," << colrateRF.nuH2i[im] << "," << colrateRF.nuH3i[im] << endl;
    }
    myOutputFile.close();

    return;
}

void handleOutput() {
    std::ifstream myInputFile("Tomator_exe/data/TOMATOR/Qs.csv");
    if (!myInputFile.is_open()) {
        cout << "\033[1;31m"
             << "[Error]: Could not open file: "
             << "Qs.csv"
             << "\033[0m" << endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    // Adjust the vector of vectors to hold 8 additional power deposition vectors
    std::vector<std::vector<double>> allData(11); // For Rpos, pd_e, pd_Hi, pd_H2i, pd_H3i, pd_Di, pd_HDi, pd_D2i, pd_D3i, pd_H2Di, pd_HD2i

    while (getline(myInputFile, line)) {
        std::istringstream iss(line);
        std::string value;
        size_t columnIndex = 0;

        while (getline(iss, value, ',')) {
            double numValue;
            try {
                numValue = std::stod(trim(value)); // Convert to double and trim if necessary
            } catch (const std::invalid_argument &e) {
                cerr << "Invalid argument for stod, value: " << value << " at column index " << columnIndex << endl;
                exit(EXIT_FAILURE);
            }

            // Assuming the predefined order matches the structure of allData
            allData[columnIndex].push_back(numValue);
            columnIndex++;
        }
    }
    myInputFile.close();

    // cout << "Data Dimensions:" << endl;
    // for (size_t i = 0; i < allData.size(); ++i) {
    //     cout << "Column " << i + 1 << ": " << allData[i].size() << " entries" << endl;
    // }

    // // Print the first 5 values of each vector, if available
    // cout << "\nFirst 5 values in each column:" << endl;
    // for (size_t i = 0; i < allData.size(); ++i) {
    //     cout << "Column " << i + 1 << ": ";
    //     for (size_t j = 0; j < std::min(allData[i].size(), size_t(5)); ++j) {
    //         cout << allData[i][j] << (j < std::min(allData[i].size(), size_t(5)) - 1 ? ", " : "");
    //     }
    //     cout << endl;
    // }

    Eigen::VectorXd originalRadialPositions(allData[0].size());
    for (size_t i = 0; i < allData[0].size(); ++i) {
        originalRadialPositions(i) = allData[0][i] * 1.00e2; // Convert to cm
    }
    Eigen::VectorXd newRadialPositions(NMESHP);

    // NewRadialPositions has the same values as aR
    for (int i = 0; i < NMESHP; ++i) {
        newRadialPositions(i) = aR[i];
    }

    std::vector<std::vector<double>> interpolatedAllData(11, std::vector<double>(NMESHP, 0.0));

    for (size_t j = 0; j < NMESHP; ++j) {
        interpolatedAllData[0][j] = newRadialPositions(j);
    }

    // Continue with interpolation for the rest of the data columns
    for (size_t i = 1; i < allData.size(); ++i) {
        Eigen::VectorXd originalData(allData[i].size());
        for (size_t j = 0; j < allData[i].size(); ++j) {
            originalData(j) = allData[i][j];
        }

        Eigen::VectorXd newData = interpolateData(newRadialPositions, originalRadialPositions, originalData);
        for (size_t j = 0; j < NMESHP; ++j) {
            interpolatedAllData[i][j] = newData(j);
        }
    }

    // // Print the dimensions of the interpolated data
    // cout << "\nInterpolated Data Dimensions:" << endl;
    // for (size_t i = 0; i < interpolatedAllData.size(); ++i) {
    //     cout << "Column " << i + 1 << ": " << interpolatedAllData[i].size() << " entries" << endl;
    // }

    // // Print the first 5 values of each vector, if available
    // cout << "\nFirst 5 values in each column after interpolation:" << endl;
    // for (size_t i = 0; i < interpolatedAllData.size(); ++i) {
    //     cout << "Column " << i + 1 << ": ";
    //     for (size_t j = 0; j < std::min(interpolatedAllData[i].size(), size_t(5)); ++j) {
    //         cout << interpolatedAllData[i][j] << (j < std::min(interpolatedAllData[i].size(), size_t(5)) - 1 ? ", " : "");
    //     }
    //     cout << endl;
    // }

    // We have pd_e, pd_Hi, pd_H2i, pd_H3i, pd_Di, pd_HDi, pd_D2i, pd_D3i, pd_H2Di, pd_HD2i
    // Put the interpolated data in the correct arrays
    for (int id = 0; id < NMESHP; id++) {
        PRFe_array[id] = interpolatedAllData[1][id];
        PRFHi_array[id] = interpolatedAllData[2][id] + interpolatedAllData[5][id]; // pd_Hi, pd_Di
        PRFH2i_array[id] = interpolatedAllData[3][id] + interpolatedAllData[6][id] + interpolatedAllData[7][id]; // pd_H2i, pd_HDi, pd_D2i
        PRFH3i_array[id] = interpolatedAllData[4][id] + interpolatedAllData[8][id] + interpolatedAllData[9][id] + interpolatedAllData[10][id]; // pd_H3i, pd_D3i, pd_H2Di, pd_HD2i
        PRFHeII_array[id] = 0;
        PRFHeIII_array[id] = 0;
        // PRFDi_array[id] = interpolatedAllData[5][id];
        // PRFHDi_array[id] = interpolatedAllData[6][id];
        // PRFD2i_array[id] = interpolatedAllData[7][id];
        // PRFD3i_array[id] = interpolatedAllData[8][id];
    }

    // Integrate
    double inteprof = 0.0;
    for (int id = 0; id < NMESHP - 1; ++id) {
        inteprof += 2.0 * b * pi * 0.5 * (PRFe_array[id] + PRFe_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (PRFHi_array[id] + PRFHi_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (PRFH2i_array[id] + PRFH2i_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (PRFH3i_array[id] + PRFH3i_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (PRFHeII_array[id] + PRFHeII_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (PRFHeIII_array[id] + PRFHeIII_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    }

    // Normalize so that the some of all integrals is Prf
    double norm = 6.24e18 * Prf * 1e3 / inteprof;
    for (int id = 0; id < NMESHP; id++) {
        PRFe_array[id] = PRFe_array[id] * norm;
        PRFHi_array[id] = PRFHi_array[id] * norm;
        PRFH2i_array[id] = PRFH2i_array[id] * norm;
        PRFH3i_array[id] = PRFH3i_array[id] * norm;
        PRFHeII_array[id] = PRFHeII_array[id] * norm;
        PRFHeIII_array[id] = PRFHeIII_array[id] * norm;
    }

    return;
}

void bICWC_func() {
    using std::chrono::duration_cast;
    using std::chrono::milliseconds;
    using std::chrono::steady_clock;

    auto start1 = steady_clock::now();

    generateInput();

    // Call the TOMATOR executable
    // TODO: doing cd and then singularity exec is not the best way to do this but the way the executable is done it is necessary to run it on its directory...
    auto start = steady_clock::now();                                                   // Start timing 
    int temp = system("cd Tomator_exe && singularity exec ubuntu_22.04.sif ./TOMATOR > /dev/null 2>&1");
    if (temp == -1) {
        cout << "\033[1;31m"
             << "[Error]: Could not execute TOMATOR"
             << "\033[0m" << endl;
        exit(EXIT_FAILURE);
    }
    auto end = steady_clock::now(); // End timing

    handleOutput();

    auto end1 = steady_clock::now(); // End timing
    auto executionTime = duration_cast<milliseconds>(end1 - start1);
    auto executableTime = duration_cast<milliseconds>(end - start);

    std::cout << "bICWC execution completed in " << executionTime.count() << " milliseconds.";
    std::cout << "\t|\tExecutable ran in " << executableTime.count() << " milliseconds.";

    // Calculate percentage of time spent in the executable
    double percentage = (executableTime.count() / static_cast<double>(executionTime.count())) * 100;
    printf("\t|\tPercentage of time spent in the executable: %.2f%\%\n", percentage);

    return;
}
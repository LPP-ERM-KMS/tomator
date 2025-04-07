#include "coupledpower.h"

void bmanuel_func(bool firstpass) {
    if (!firstpass) {
        return;
    }
    cout << "Running bmanuel_func" << endl;

    std::ifstream myInputFile(smanuel);
    if (!myInputFile.is_open()) {
        cout << "\033[1;31m"
             << "[Error]: Could not open file: " << smanuel << "\033[0m" << endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<double> radialPositions;
    std::vector<std::vector<double>> allData;
    bool isFirstLine = true;
    std::map<std::string, size_t> headerMap;
    size_t columnIndex = 0;

    while (getline(myInputFile, line)) {
        std::istringstream iss(line);
        std::string value;
        std::vector<double> lineData;

        if (isFirstLine) {
            // Process headers
            while (getline(iss, value, ',')) {
                std::string trimmedValue = trim(value); // Trim the header string
                headerMap[trimmedValue] = columnIndex++;
            }
            
            // Ensure necessary headers are present
            std::vector<std::string> requiredHeaders = {"RadialPositions", "PRFe", "PRFHi", "PRFH2i", "PRFH3i", "PRFHeII", "PRFHeIII"};
            for (const auto &header : requiredHeaders) {
                if (headerMap.find(header) == headerMap.end()) {
                    cout << "\033[1;31m"
                         << "[Error]: Could not find " << header << " in file: " << smanuel << "\033[0m" << endl;
                    exit(EXIT_FAILURE);
                }
            }
            isFirstLine = false;
            continue;
        }

        columnIndex = 0;
        while (getline(iss, value, ',')) {
            if (columnIndex == 0) { // Assuming RadialPositions is the first column
                try {
                    radialPositions.push_back(std::stod(value));
                } catch (const std::invalid_argument &e) {
                    cerr << "Invalid argument for stod, value: " << value << " at RadialPositions column" << endl;
                    exit(EXIT_FAILURE);
                }
            } else {
                try {
                    lineData.push_back(std::stod(value));
                } catch (const std::invalid_argument &e) {
                    cerr << "Invalid argument for stod, value: " << value << " at column index " << columnIndex << endl;
                    exit(EXIT_FAILURE);
                }
            }
            columnIndex++;
        }
        if (!lineData.empty()) {
            allData.push_back(lineData);
        }
    }
    myInputFile.close();

    size_t numberOfColumns = allData[0].size();
    std::vector<std::vector<double>> interpolatedAllData(NMESHP, std::vector<double>(numberOfColumns, 0.0));

    size_t dataPoints = radialPositions.size();
    if (dataPoints != NMESHP) {
        // Perform Interpolation
        Eigen::VectorXd originalRadialPositions(dataPoints);
        for (size_t i = 0; i < dataPoints; ++i) {
            originalRadialPositions(i) = radialPositions[i];
        }

        Eigen::VectorXd newRadialPositions = Eigen::VectorXd::LinSpaced(NMESHP, originalRadialPositions.minCoeff(), originalRadialPositions.maxCoeff());

        for (const auto &header : headerMap) {
            if (header.first == "RadialPositions")
                continue;

            Eigen::VectorXd originalDataVector(dataPoints);
            for (size_t i = 0; i < dataPoints; ++i) {
                originalDataVector(i) = allData[i][header.second - 1]; // -1 because RadialPositions is not in allData
            }

            Eigen::VectorXd interpolatedData = interpolateData(newRadialPositions, originalRadialPositions, originalDataVector);

            // Update allData with interpolatedData
            for (int i = 0; i < NMESHP; ++i) {
                interpolatedAllData[i][header.second - 1] = interpolatedData(i);
            }
        }

        // Update radialPositions with newRadialPositions
        radialPositions = std::vector<double>(newRadialPositions.data(), newRadialPositions.data() + newRadialPositions.size());
    }
    else {
        // No interpolation needed
        interpolatedAllData = allData;
    }

    // Integrate
    double inteprof = 0.0;

    // print interpolated data sizes
    // cout << "Interpolated data sizes: " << interpolatedAllData.size() << " " << interpolatedAllData[0].size() << endl;
    // for (int i = 0; i < interpolatedAllData.size(); ++i) {
    //     cout << aR[i] << " " << interpolatedAllData[i][headerMap["PRFe"]-1] << " " << interpolatedAllData[i][headerMap["PRFHi"]-1] << " " << interpolatedAllData[i][headerMap["PRFH2i"]-1] << " " << interpolatedAllData[i][headerMap["PRFH3i"]-1] << " " << interpolatedAllData[i][headerMap["PRFHeII"]-1] << " " << interpolatedAllData[i][headerMap["PRFHeIII"]-1] << endl;
    // }


    for (int id = 0; id < NMESHP - 1; ++id) {
        inteprof += 2.0 * b * pi * 0.5 * (interpolatedAllData[id][headerMap["PRFe"]-1] + interpolatedAllData[id + 1][headerMap["PRFe"]-1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (interpolatedAllData[id][headerMap["PRFHi"]-1] + interpolatedAllData[id + 1][headerMap["PRFHi"]-1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (interpolatedAllData[id][headerMap["PRFH2i"]-1] + interpolatedAllData[id + 1][headerMap["PRFH2i"]-1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (interpolatedAllData[id][headerMap["PRFH3i"]-1] + interpolatedAllData[id + 1][headerMap["PRFH3i"]-1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (interpolatedAllData[id][headerMap["PRFHeII"]-1] + interpolatedAllData[id + 1][headerMap["PRFHeII"]-1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof += 2.0 * b * pi * 0.5 * (interpolatedAllData[id][headerMap["PRFHeIII"]-1] + interpolatedAllData[id + 1][headerMap["PRFHeIII"]-1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    }

    // Normalize so that the some of all integrals is Prf
    double norm = 6.24e18 * (Prf * 1e3) / inteprof;
    for (int id = 0; id < NMESHP; id++) {
        PRFe_array[id] = interpolatedAllData[id][headerMap["PRFe"]-1] * norm;
        PRFHi_array[id] = interpolatedAllData[id][headerMap["PRFHi"]-1] * norm;
        PRFH2i_array[id] = interpolatedAllData[id][headerMap["PRFH2i"]-1] * norm;
        PRFH3i_array[id] = interpolatedAllData[id][headerMap["PRFH3i"]-1] * norm;
        PRFHeII_array[id] = interpolatedAllData[id][headerMap["PRFHeII"]-1] * norm;
        PRFHeIII_array[id] = interpolatedAllData[id][headerMap["PRFHeIII"]-1] * norm;
    }

    // Check if the total power is Prf
    // inteprof = 0.0;
    // for (int id = 0; id < NMESHP - 1; ++id) {
    //     inteprof += 2.0 * b * pi * 0.5 * (PRFe_array[id] + PRFe_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    //     inteprof += 2.0 * b * pi * 0.5 * (PRFHi_array[id] + PRFHi_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    //     inteprof += 2.0 * b * pi * 0.5 * (PRFH2i_array[id] + PRFH2i_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    //     inteprof += 2.0 * b * pi * 0.5 * (PRFH3i_array[id] + PRFH3i_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    //     inteprof += 2.0 * b * pi * 0.5 * (PRFHeII_array[id] + PRFHeII_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    //     inteprof += 2.0 * b * pi * 0.5 * (PRFHeIII_array[id] + PRFHeIII_array[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    // }
    // cout << "Difference between Prf and total power: " << Prf - inteprof << endl;

    return;
}

#include "coupledpower.h"

Eigen::VectorXd interpolateData(const Eigen::VectorXd &newPositions, const Eigen::VectorXd &originalPositions, const Eigen::VectorXd &originalData) {
    Eigen::VectorXd newData(newPositions.size());
    for (int i = 0; i < newPositions.size(); ++i) {
        // Find the closest points
        int k = 0;
        for (; k < originalPositions.size() - 1; ++k) {
            if (newPositions(i) < originalPositions(k + 1)) {
                break;
            }
        }
        double t = (newPositions(i) - originalPositions(k)) / (originalPositions(k + 1) - originalPositions(k));
        newData(i) = originalData(k) + t * (originalData(k + 1) - originalData(k)); // Linear interpolation
    }
    return newData;
}

std::string trim(const std::string &str) {
    std::string whitespace = " \t\n\r";
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

void coupledpower(const double &freq, double &alr) //, double PTEST[] )
{
    // static double alphaval, ionalpha;
    // alphaval=(nr.ne[centerval])/(nH2[centerval] + nHeI[centerval] + nr.ne[centerval]);
    // the RF module starts working only from about 1e9cm-3
    // we split the coupled power in 2 parts, if ne is below 1e9, then we put the power proportional to the electron density,
    // if it is above then we use the RF module and ensure for a smooth transition...
    // if (bkipt_decay | bkipt_propto)

    // cout << "bkipt = " << bkipt << endl;
    // cout << "bproptone = " << bproptone << endl;
    // cout << "bgray = " << bgray << endl;
    // cout << "bram = " << bram << endl;
    // cout << "bnefix = " << bnefix << endl;
    // cout << "bfixpowerfrac = " << bfixpowerfrac << endl;
    // cout << "bnopower = " << bnopower << endl;
    // cout << "bTOMAS = " << bTOMAS << endl;

    // cout << "bmanuel = " << bmanuel << endl;
    // cout << "bICWC = " << bICWC << endl;
    // exit(EXIT_SUCCESS);

    //if (bkipt) { // KIPT module
    //    bkipt_func(alr);
    //}

    /////// Prop-to-ne, cylindrical coordinates
    if (bproptone) {
        bproptone_func();
    }

    /////// Gray - central deposition
    if (bgray) {
        bgray_func();
    }

    ///// Tulchhi Ram IPR - central deposition
    if (bram) {
        bram_func();
    }

    ///// Tune coupled power for fixed density at ic
    if (bnefix) {
        bnefix_func();
    }

    ///// couple fixed fraction of launched power
    if (bfixpowerfrac) {
        bfixpowerfrac_func();
    }

    ///// couple fixed fraction of launched power
    if (blhr) {
        bicatlhr_func();
    }

    /////// no power
    if (bnopower) {
        bnopower_func();
    }

    //////////////////////////////////////////
    /////// TOMAS - deposition on ECR ////////
    //////////////////////////////////////////

    if (bTOMAS) {
        bTOMAS_func();
    }

    //////////////////////////////////////////
    /////// Manuel - deposition on ECR ///////
    //////////////////////////////////////////

    if (bmanuel) {
        bmanuel_func(Nit == 0);
    }

    if (bICWC) {
        bICWC_func();
    }

    if ((tmain < (1.02 * dtpramp)) & !(bmanuel | bICWC)) {
        for (int id = 0; id < NMESHP; ++id) {
            PRFe_array[id] = PRFe_array[id] * ((tmain + 0.02 * dtpramp) / (1.02 * dtpramp));
            PRFHi_array[id] = PRFHi_array[id] * ((tmain + 0.02 * dtpramp) / (1.02 * dtpramp));
            PRFH2i_array[id] = PRFH2i_array[id] * ((tmain + 0.02 * dtpramp) / (1.02 * dtpramp));
            PRFH3i_array[id] = PRFH3i_array[id] * ((tmain + 0.02 * dtpramp) / (1.02 * dtpramp));
            PRFHeII_array[id] = PRFHeII_array[id] * ((tmain + 0.02 * dtpramp) / (1.02 * dtpramp));
            PRFHeIII_array[id] = PRFHeIII_array[id] * ((tmain + 0.02 * dtpramp) / (1.02 * dtpramp));
        }
    }
}

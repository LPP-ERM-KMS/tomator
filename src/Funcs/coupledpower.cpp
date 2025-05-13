
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

void coupledpower() 
{

    /////// Prop-to-ne, cylindrical coordinates
    if (bproptone) {
        bproptone_func();
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

    if ((tmain < (1.02 * dtpramp)) & !(bmanuel)) {
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

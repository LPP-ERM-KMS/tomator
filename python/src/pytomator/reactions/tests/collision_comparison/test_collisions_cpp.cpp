/**
 * C++ Collision Sources Test Driver
 * 
 * This program calls the actual collisions() function from Tomator1D 
 * with controlled test conditions and outputs the source terms (dn, dE)
 * for direct comparison with Python's compute_collision_sources().
 * 
 * Build: See build_test_collisions.sh
 * 
 * Usage: ./test_collisions_cpp <hydhel_path> [Te] [ne]
 * 
 * Output: JSON-formatted collision source terms for all species
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

// Include Tomator headers (relative paths from collision_comparison directory)
#include "../../../../../src/Vars/datastructures.h"
#include "../../../../../src/Vars/simparam.h"
#include "../../../../../src/Vars/globalVariables.h"
#include "../../../../../src/Vars/reactionrates.h"
#include "../../../../../src/Funcs/collisions.h"

// Function to initialize HYDHEL data
extern void initReactionDataMap(char *fileName);

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <hydhel.tex path> [Te=10] [ne_cgs=1e12] [bion=1] [bcx=1] [belas=1] [bcoulomb=1]" << std::endl;
        return 1;
    }
    
    char* hydhelPath = argv[1];
    double Te_test = (argc > 2) ? std::stod(argv[2]) : 10.0;
    double ne_test = (argc > 3) ? std::stod(argv[3]) : 1e12;
    bool arg_bion = (argc > 4) ? (std::stoi(argv[4]) != 0) : true;
    bool arg_bcx = (argc > 5) ? (std::stoi(argv[5]) != 0) : true;
    bool arg_belas = (argc > 6) ? (std::stoi(argv[6]) != 0) : true;
    bool arg_bcoulomb = (argc > 7) ? (std::stoi(argv[7]) != 0) : true;
    bool arg_bH = (argc > 8) ? (std::stoi(argv[8]) != 0) : true;
    bool arg_bH2 = (argc > 9) ? (std::stoi(argv[9]) != 0) : true;
    bool arg_bHe = (argc > 10) ? (std::stoi(argv[10]) != 0) : true;
    
    // Initialize HYDHEL reaction data
    std::cerr << "Initializing HYDHEL data from: " << hydhelPath << std::endl;
    initReactionDataMap(hydhelPath);
    
    // Set up test conditions (same as Python comparison)
    double TH = 0.5;
    double THi = Te_test * 0.8;
    double TH2 = TH;
    double TH2i = THi;
    double TH3i = THi;
    double THeI = TH;
    double THeII = THi;
    double THeIII = THi;
    
    double nH_cgs = 1e10;
    double nHi_cgs = ne_test * 0.9;
    double nH2_cgs = 1e11;
    double nH2i_cgs = 1e9;
    double nH3i_cgs = 1e8;
    double nHeI_cgs = 1e12;
    double nHeII_cgs = ne_test * 0.1;
    double nHeIII_cgs = 1e6;
    
    // Set simulation flags (match Python comparison settings)
    bH = arg_bH;
    bH2 = arg_bH2;
    bHe = arg_bHe;
    bADAS = true;  // Use ADAS cooling to match Python
    bion = arg_bion;
    bcx = arg_bcx;
    belas = arg_belas;
    bcoulomb = arg_bcoulomb;
    
    // Initialize all mesh points to zero first
    for (int im = 0; im < NMESHP; ++im) {
        // Densities (CGS: cm^-3)
        nr.ne[im] = 0.0;
        nr.nH[im] = 0.0;
        nr.nH2[im] = 0.0;
        nr.nHi[im] = 0.0;
        nr.nH2i[im] = 0.0;
        nr.nH3i[im] = 0.0;
        nr.nHeI[im] = 0.0;
        nr.nHeII[im] = 0.0;
        nr.nHeIII[im] = 0.0;
        nr.nCI[im] = 0.0;
        nr.nCII[im] = 0.0;
        nr.nCIII[im] = 0.0;
        nr.nCIV[im] = 0.0;
        nr.nCV[im] = 0.0;
        
        // Temperatures (eV)
        Tr.Te[im] = 0.1;
        Tr.TH[im] = 0.1;
        Tr.TH2[im] = 0.1;
        Tr.THi[im] = 0.1;
        Tr.TH2i[im] = 0.1;
        Tr.TH3i[im] = 0.1;
        Tr.THeI[im] = 0.1;
        Tr.THeII[im] = 0.1;
        Tr.THeIII[im] = 0.1;
    }
    
    // Set test conditions at mesh point 0
    int im = 0;
    nr.ne[im] = ne_test;
    nr.nH[im] = nH_cgs;
    nr.nHi[im] = nHi_cgs;
    nr.nH2[im] = nH2_cgs;
    nr.nH2i[im] = nH2i_cgs;
    nr.nH3i[im] = nH3i_cgs;
    nr.nHeI[im] = nHeI_cgs;
    nr.nHeII[im] = nHeII_cgs;
    nr.nHeIII[im] = nHeIII_cgs;
    
    Tr.Te[im] = Te_test;
    Tr.TH[im] = TH;
    Tr.THi[im] = THi;
    Tr.TH2[im] = TH2;
    Tr.TH2i[im] = TH2i;
    Tr.TH3i[im] = TH3i;
    Tr.THeI[im] = THeI;
    Tr.THeII[im] = THeII;
    Tr.THeIII[im] = THeIII;
    
    std::cerr << "Running collisions() with:" << std::endl;
    std::cerr << "  Te = " << Te_test << " eV" << std::endl;
    std::cerr << "  ne = " << ne_test << " cm^-3" << std::endl;
    std::cerr << "  bH = " << bH << ", bH2 = " << bH2 << ", bHe = " << bHe << ", bADAS = " << bADAS << std::endl;
    std::cerr << "  bion = " << bion << ", bcx = " << bcx << ", belas = " << belas << ", bcoulomb = " << bcoulomb << std::endl;
    
    // Call the actual C++ collisions function
    collisions();
    
    // Output results as JSON to stdout
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "{" << std::endl;
    std::cout << "  \"test_conditions\": {" << std::endl;
    std::cout << "    \"Te\": " << Te_test << "," << std::endl;
    std::cout << "    \"ne_cgs\": " << ne_test << "," << std::endl;
    std::cout << "    \"nH_cgs\": " << nH_cgs << "," << std::endl;
    std::cout << "    \"nHi_cgs\": " << nHi_cgs << "," << std::endl;
    std::cout << "    \"nH2_cgs\": " << nH2_cgs << "," << std::endl;
    std::cout << "    \"nH2i_cgs\": " << nH2i_cgs << "," << std::endl;
    std::cout << "    \"nH3i_cgs\": " << nH3i_cgs << "," << std::endl;
    std::cout << "    \"nHeI_cgs\": " << nHeI_cgs << "," << std::endl;
    std::cout << "    \"nHeII_cgs\": " << nHeII_cgs << "," << std::endl;
    std::cout << "    \"nHeIII_cgs\": " << nHeIII_cgs << "," << std::endl;
    std::cout << "    \"TH\": " << TH << "," << std::endl;
    std::cout << "    \"THi\": " << THi << "," << std::endl;
    std::cout << "    \"TH2\": " << TH2 << "," << std::endl;
    std::cout << "    \"TH2i\": " << TH2i << "," << std::endl;
    std::cout << "    \"TH3i\": " << TH3i << "," << std::endl;
    std::cout << "    \"THeI\": " << THeI << "," << std::endl;
    std::cout << "    \"THeII\": " << THeII << "," << std::endl;
    std::cout << "    \"THeIII\": " << THeIII << "," << std::endl;
    std::cout << "    \"bH2\": " << bH2 << "," << std::endl;
    std::cout << "    \"bHe\": " << bHe << "," << std::endl;
    std::cout << "    \"bADAS\": " << bADAS << "," << std::endl;
    std::cout << "    \"bion\": " << bion << "," << std::endl;
    std::cout << "    \"bcx\": " << bcx << std::endl;
    std::cout << "  }," << std::endl;
    
    std::cout << "  \"dn_cgs\": {" << std::endl;
    std::cout << "    \"e\": " << dnr.dne[im] << "," << std::endl;
    std::cout << "    \"H\": " << dnr.dnH[im] << "," << std::endl;
    std::cout << "    \"Hi\": " << dnr.dnHi[im] << "," << std::endl;
    std::cout << "    \"H2\": " << dnr.dnH2[im] << "," << std::endl;
    std::cout << "    \"H2i\": " << dnr.dnH2i[im] << "," << std::endl;
    std::cout << "    \"H3i\": " << dnr.dnH3i[im] << "," << std::endl;
    std::cout << "    \"HeI\": " << dnr.dnHeI[im] << "," << std::endl;
    std::cout << "    \"HeII\": " << dnr.dnHeII[im] << "," << std::endl;
    std::cout << "    \"HeIII\": " << dnr.dnHeIII[im] << std::endl;
    std::cout << "  }," << std::endl;
    
    std::cout << "  \"dE_cgs\": {" << std::endl;
    std::cout << "    \"e\": " << dEr.dEe[im] << "," << std::endl;
    std::cout << "    \"H\": " << dEr.dEH[im] << "," << std::endl;
    std::cout << "    \"Hi\": " << dEr.dEHi[im] << "," << std::endl;
    std::cout << "    \"H2\": " << dEr.dEH2[im] << "," << std::endl;
    std::cout << "    \"H2i\": " << dEr.dEH2i[im] << "," << std::endl;
    std::cout << "    \"H3i\": " << dEr.dEH3i[im] << "," << std::endl;
    std::cout << "    \"HeI\": " << dEr.dEHeI[im] << "," << std::endl;
    std::cout << "    \"HeII\": " << dEr.dEHeII[im] << "," << std::endl;
    std::cout << "    \"HeIII\": " << dEr.dEHeIII[im] << std::endl;
    std::cout << "  }" << std::endl;
    std::cout << "}" << std::endl;
    
    return 0;
}

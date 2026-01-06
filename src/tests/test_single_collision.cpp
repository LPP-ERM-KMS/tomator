/**
 * test_single_collision.cpp
 * 
 * Standalone test to run collisions() with a single flag enabled.
 * Compiles against the Tomator source code.
 * 
 * Compile from tomator/src/build:
 *   g++ -O2 -std=c++17 -fopenmp -I.. -I../Eigen \
 *       ../tests/test_single_collision.cpp \
 *       ../Funcs/collisions.cpp ../Funcs/functions.cpp \
 *       ../SimParams/simparam.cpp ../Vars/globalVariables.cpp \
 *       ../Funcs/fileHandler.cpp ../Funcs/meanfilter.cpp \
 *       -o test_single_collision
 * 
 * Usage: ./test_single_collision <flag_name>
 */

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>

// Include Tomator headers
#include "../Vars/globalVariables.h"
#include "../Vars/simparam.h"
#include "../Vars/reactionrates.h"
#include "../Funcs/collisions.h"

using namespace std;

// Test conditions
const double TEST_Te = 10.0;       // eV
const double TEST_ne = 1e12;       // cm^-3 (C++ uses CGS internally)
const double TEST_TH = 0.5;        // eV
const double TEST_nH = 1e10;       // cm^-3
const double TEST_THi = 8.0;       // eV
const double TEST_nHi = 9e11;      // cm^-3
const double TEST_TH2 = 0.3;       // eV
const double TEST_nH2 = 5e9;       // cm^-3
const double TEST_TH2i = 2.0;      // eV
const double TEST_nH2i = 1e9;      // cm^-3
const double TEST_TH3i = 1.5;      // eV
const double TEST_nH3i = 1e8;      // cm^-3
const double TEST_THeI = 0.5;      // eV
const double TEST_nHeI = 1e10;     // cm^-3
const double TEST_THeII = 5.0;     // eV
const double TEST_nHeII = 5e9;     // cm^-3
const double TEST_THeIII = 10.0;   // eV
const double TEST_nHeIII = 1e8;    // cm^-3

void set_all_flags_false() {
    // Main section flags - keep true so individual flags control reactions
    bH = true;
    bH2 = true;
    bHe = true;
    bcx = true;
    bion = true;
    belas = true;
    bcoulomb = false;  // Disable for cleaner comparison
    bADAS = true;      // Use ADAS tables to match Python default
    
    // Individual bH flags
    bH_exc = false;
    bH_ion = false;
    bH_3body = false;
    bH_rec = false;
    
    // Individual bH2 flags
    bH2_elas = false;
    bH2_exc = false;
    bH2_dis = false;
    bH2_ion = false;
    bH2i_rec = false;
    bH2_dision = false;
    bH2i_dis = false;
    bH2i_disexc = false;
    bH2i_disrec = false;
    bH3i_disrec = false;
    bH3i_dis = false;
    
    // Individual bHe flags
    bHeI_ion = false;
    bHeII_ion = false;
    bHe_cooling = false;
    bHeII_rec = false;
    bHeIII_rec = false;
    
    // Individual bcx flags
    bcx_HiH = false;
    bcx_HiH2 = false;
    bcx_H2iH2 = false;
    bcx_HeIIH = false;
    bcx_HeIIHeI = false;
    bcx_HeIIIH = false;
    bcx_HeIIIHeI = false;
    
    // Individual bion flags
    bion_HiH_exca = false;
    bion_HiH_excb = false;
    bion_HiH_ion = false;
    bion_HiH2_exca = false;
    bion_HiH2_excb = false;
    bion_HiH2_325 = false;
    bion_HiH2i_326 = false;
    bion_H2iH2_H3i = false;
    bion_HiHeI_ion = false;
    bion_HeIIH2_cxdis = false;
    
    // Individual belas flags
    belas_HiH = false;
    belas_HH = false;
    belas_HH2 = false;
    belas_HiH2 = false;
    belas_H2iH = false;
    belas_H2iH2 = false;
    belas_H3iH = false;
    belas_H3iH2 = false;
    belas_H2H2 = false;
    belas_HHeI = false;
    belas_HeIIH = false;
    belas_HiHeI = false;
    belas_HeIIHeI = false;
    belas_HeIHeI = false;
    belas_HeIH2 = false;
    belas_HeIIH2 = false;
}

bool set_flag_by_name(const string& flag_name) {
    // bH flags
    if (flag_name == "bH_exc") { bH_exc = true; return true; }
    if (flag_name == "bH_ion") { bH_ion = true; return true; }
    if (flag_name == "bH_3body") { bH_3body = true; return true; }
    if (flag_name == "bH_rec") { bH_rec = true; return true; }
    
    // bH2 flags
    if (flag_name == "bH2_elas") { bH2_elas = true; return true; }
    if (flag_name == "bH2_exc") { bH2_exc = true; return true; }
    if (flag_name == "bH2_dis") { bH2_dis = true; return true; }
    if (flag_name == "bH2_ion") { bH2_ion = true; return true; }
    if (flag_name == "bH2i_rec") { bH2i_rec = true; return true; }
    if (flag_name == "bH2_dision") { bH2_dision = true; return true; }
    if (flag_name == "bH2i_dis") { bH2i_dis = true; return true; }
    if (flag_name == "bH2i_disexc") { bH2i_disexc = true; return true; }
    if (flag_name == "bH2i_disrec") { bH2i_disrec = true; return true; }
    if (flag_name == "bH3i_disrec") { bH3i_disrec = true; return true; }
    if (flag_name == "bH3i_dis") { bH3i_dis = true; return true; }
    
    // bHe flags
    if (flag_name == "bHeI_ion") { bHeI_ion = true; return true; }
    if (flag_name == "bHeII_ion") { bHeII_ion = true; return true; }
    if (flag_name == "bHe_cooling") { bHe_cooling = true; return true; }
    if (flag_name == "bHeII_rec") { bHeII_rec = true; return true; }
    if (flag_name == "bHeIII_rec") { bHeIII_rec = true; return true; }
    
    // bcx flags
    if (flag_name == "bcx_HiH") { bcx_HiH = true; return true; }
    if (flag_name == "bcx_HiH2") { bcx_HiH2 = true; return true; }
    if (flag_name == "bcx_H2iH2") { bcx_H2iH2 = true; return true; }
    if (flag_name == "bcx_HeIIH") { bcx_HeIIH = true; return true; }
    if (flag_name == "bcx_HeIIHeI") { bcx_HeIIHeI = true; return true; }
    if (flag_name == "bcx_HeIIIH") { bcx_HeIIIH = true; return true; }
    if (flag_name == "bcx_HeIIIHeI") { bcx_HeIIIHeI = true; return true; }
    
    // bion flags
    if (flag_name == "bion_HiH_exca") { bion_HiH_exca = true; return true; }
    if (flag_name == "bion_HiH_excb") { bion_HiH_excb = true; return true; }
    if (flag_name == "bion_HiH_ion") { bion_HiH_ion = true; return true; }
    if (flag_name == "bion_HiH2_exca") { bion_HiH2_exca = true; return true; }
    if (flag_name == "bion_HiH2_excb") { bion_HiH2_excb = true; return true; }
    if (flag_name == "bion_HiH2_325") { bion_HiH2_325 = true; return true; }
    if (flag_name == "bion_HiH2i_326") { bion_HiH2i_326 = true; return true; }
    if (flag_name == "bion_H2iH2_H3i") { bion_H2iH2_H3i = true; return true; }
    if (flag_name == "bion_HiHeI_ion") { bion_HiHeI_ion = true; return true; }
    if (flag_name == "bion_HeIIH2_cxdis") { bion_HeIIH2_cxdis = true; return true; }
    
    // belas flags
    if (flag_name == "belas_HiH") { belas_HiH = true; return true; }
    if (flag_name == "belas_HH") { belas_HH = true; return true; }
    if (flag_name == "belas_HH2") { belas_HH2 = true; return true; }
    if (flag_name == "belas_HiH2") { belas_HiH2 = true; return true; }
    if (flag_name == "belas_H2iH") { belas_H2iH = true; return true; }
    if (flag_name == "belas_H2iH2") { belas_H2iH2 = true; return true; }
    if (flag_name == "belas_H3iH") { belas_H3iH = true; return true; }
    if (flag_name == "belas_H3iH2") { belas_H3iH2 = true; return true; }
    if (flag_name == "belas_H2H2") { belas_H2H2 = true; return true; }
    if (flag_name == "belas_HHeI") { belas_HHeI = true; return true; }
    if (flag_name == "belas_HeIIH") { belas_HeIIH = true; return true; }
    if (flag_name == "belas_HiHeI") { belas_HiHeI = true; return true; }
    if (flag_name == "belas_HeIIHeI") { belas_HeIIHeI = true; return true; }
    if (flag_name == "belas_HeIHeI") { belas_HeIHeI = true; return true; }
    if (flag_name == "belas_HeIH2") { belas_HeIH2 = true; return true; }
    if (flag_name == "belas_HeIIH2") { belas_HeIIH2 = true; return true; }
    
    return false;  // Unknown flag
}

void setup_test_conditions() {
    // Set simulation parameters
    nmeshp = 1;  // Single grid point for testing
    dt = 1e-6;
    accur = 0.1;
    
    // Initialize density and temperature arrays at grid point 0
    // C++ uses CGS (cm^-3)
    nr.ne[0] = TEST_ne;
    nr.nH[0] = TEST_nH;
    nr.nHi[0] = TEST_nHi;
    nr.nH2[0] = TEST_nH2;
    nr.nH2i[0] = TEST_nH2i;
    nr.nH3i[0] = TEST_nH3i;
    nr.nHeI[0] = TEST_nHeI;
    nr.nHeII[0] = TEST_nHeII;
    nr.nHeIII[0] = TEST_nHeIII;
    
    Tr.Te[0] = TEST_Te;
    Tr.TH[0] = TEST_TH;
    Tr.THi[0] = TEST_THi;
    Tr.TH2[0] = TEST_TH2;
    Tr.TH2i[0] = TEST_TH2i;
    Tr.TH3i[0] = TEST_TH3i;
    Tr.THeI[0] = TEST_THeI;
    Tr.THeII[0] = TEST_THeII;
    Tr.THeIII[0] = TEST_THeIII;
}

void output_json() {
    // Output results in JSON format (CGS units: cm^-3/s for dn, eV*cm^-3/s for dE)
    cout << fixed << setprecision(10);
    cout << "{" << endl;
    
    // dn (density changes)
    cout << "  \"dn\": {" << endl;
    cout << "    \"e\": " << dnr.dne[0] << "," << endl;
    cout << "    \"H\": " << dnr.dnH[0] << "," << endl;
    cout << "    \"Hi\": " << dnr.dnHi[0] << "," << endl;
    cout << "    \"H2\": " << dnr.dnH2[0] << "," << endl;
    cout << "    \"H2i\": " << dnr.dnH2i[0] << "," << endl;
    cout << "    \"H3i\": " << dnr.dnH3i[0] << "," << endl;
    cout << "    \"HeI\": " << dnr.dnHeI[0] << "," << endl;
    cout << "    \"HeII\": " << dnr.dnHeII[0] << "," << endl;
    cout << "    \"HeIII\": " << dnr.dnHeIII[0] << endl;
    cout << "  }," << endl;
    
    // dE (energy changes)
    cout << "  \"dE\": {" << endl;
    cout << "    \"e\": " << dEr.dEe[0] << "," << endl;
    cout << "    \"H\": " << dEr.dEH[0] << "," << endl;
    cout << "    \"Hi\": " << dEr.dEHi[0] << "," << endl;
    cout << "    \"H2\": " << dEr.dEH2[0] << "," << endl;
    cout << "    \"H2i\": " << dEr.dEH2i[0] << "," << endl;
    cout << "    \"H3i\": " << dEr.dEH3i[0] << "," << endl;
    cout << "    \"HeI\": " << dEr.dEHeI[0] << "," << endl;
    cout << "    \"HeII\": " << dEr.dEHeII[0] << "," << endl;
    cout << "    \"HeIII\": " << dEr.dEHeIII[0] << endl;
    cout << "  }," << endl;
    
    // nu (collision frequencies)
    cout << "  \"nu\": {" << endl;
    cout << "    \"e\": " << colrateRF.nue[0] << "," << endl;
    cout << "    \"Hi\": " << colrateRF.nuHi[0] << "," << endl;
    cout << "    \"H2i\": " << colrateRF.nuH2i[0] << "," << endl;
    cout << "    \"H3i\": " << colrateRF.nuH3i[0] << "," << endl;
    cout << "    \"HeII\": " << colrateRF.nuHeII[0] << "," << endl;
    cout << "    \"HeIII\": " << colrateRF.nuHeIII[0] << endl;
    cout << "  }" << endl;
    
    cout << "}" << endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <flag_name>" << endl;
        cerr << "Example: " << argv[0] << " bH_ion" << endl;
        return 1;
    }
    
    string flag_name = argv[1];
    
    // Initialize reaction rate data from hydhel.tex file
    // Try multiple possible paths
    char hydhelPath[512];
    const char* possiblePaths[] = {
        "../SimParams/Public/hydhel.tex",
        "../../SimParams/Public/hydhel.tex",
        "../../../src/SimParams/Public/hydhel.tex",
        "/home/ITER/wautert/Documents/Tomator-RMA/tomator/src/SimParams/Public/hydhel.tex"
    };
    
    bool found = false;
    for (const char* path : possiblePaths) {
        FILE* f = fopen(path, "r");
        if (f) {
            fclose(f);
            strcpy(hydhelPath, path);
            found = true;
            break;
        }
    }
    
    if (!found) {
        cerr << "Error: Could not find hydhel.tex file" << endl;
        return 1;
    }
    
    initReactionDataMap(hydhelPath);
    
    // Set all individual flags to false
    set_all_flags_false();
    
    // Enable the requested flag
    if (!set_flag_by_name(flag_name)) {
        cerr << "Unknown flag: " << flag_name << endl;
        return 1;
    }
    
    // Setup test conditions
    setup_test_conditions();
    
    // Run collisions
    collisions();
    
    // Output results
    output_json();
    
    return 0;
}

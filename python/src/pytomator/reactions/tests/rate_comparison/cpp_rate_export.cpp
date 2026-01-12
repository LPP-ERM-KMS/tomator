/**
 * C++ Rate Export Driver
 * 
 * This program exports reaction rates calculated using the exact C++ implementation
 * to CSV format for comparison with Python implementations.
 * 
 * Build: See build_cpp_rate_export.sh
 * 
 * Usage: ./cpp_rate_export <hydhel_path> <output_csv>
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include "reactionrates.h"

// External variables from other source files
extern bool bADAS;       // from simparam.cpp
extern const double pi;  // from constants.c

// Function declarations from reactionrates.cpp
void initReactionDataMap(char *fileName);
double RR(int reac, double E, double T);
double RRH2(int reac, double ne, double Te);
double RRHe(int reac, double ne, double Te);
double RRCX(int reac, double T1, double T2);  // uppercase
double RRion(int reac, double T1, double T2);
double RRel(int reac, double T1, double T2);

// Analytical formula for REAC218a from collisions.cpp line 220
// This is what's actually used in the simulation (not the 2D polynomial in RR())
double REAC218a_analytical(double Te) {
    if (Te >= 1000.0) return 0.0;
    double k = 7.982e-11 / (sqrt(Te / 2.713e-4) * 
                            pow((1 + sqrt(Te / 2.713e-4)), (1.0 - 0.7480)) * 
                            pow((1 + sqrt(Te / 60.631)), (1.0 + 0.7480)));
    return k;
}

// Generate logarithmically spaced array
std::vector<double> logspace(double start, double end, int n) {
    std::vector<double> result(n);
    double logStart = std::log10(start);
    double logEnd = std::log10(end);
    double step = (logEnd - logStart) / (n - 1);
    for (int i = 0; i < n; ++i) {
        result[i] = std::pow(10.0, logStart + i * step);
    }
    return result;
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <hydhel.tex path> <output.csv>" << std::endl;
        return 1;
    }
    
    char* hydhelPath = argv[1];
    std::string outputPath = argv[2];
    
    // Initialize HYDHEL reaction data
    std::cout << "Initializing HYDHEL data from: " << hydhelPath << std::endl;
    initReactionDataMap(hydhelPath);
    
    // Temperature and density grids
    auto Te_grid = logspace(0.1, 1000.0, 100);  // 100 points from 0.1 to 1000 eV
    std::vector<double> ne_slices = {1e8, 1e11, 1e14};  // cm^-3
    std::vector<double> T2_slices = {1.0, 10.0, 100.0}; // eV
    
    // Open output file
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Cannot open output file: " << outputPath << std::endl;
        return 1;
    }
    
    // Set precision for output
    outFile << std::scientific << std::setprecision(10);
    
    // Write header
    outFile << "cpp_func,cpp_id,param1_name,param1_value,param2_name,param2_value,rate_cm3_s" << std::endl;
    
    std::cout << "Exporting reaction rates to: " << outputPath << std::endl;
    
    // =====================================================================
    // RR Reactions (HYDHEL polynomial fits) - Te only
    // =====================================================================
    std::cout << "  Exporting RR reactions..." << std::endl;
    
    // List of RR reactions (1D: Te only, E=Te for most)
    std::vector<std::pair<int, std::string>> rr_reactions_1d = {
        {REAC211, "REAC211"},
        {REAC212, "REAC212"},
        {REAC213, "REAC213"},
        {REAC214a, "REAC214a"},
        {REAC214b, "REAC214b"},
        {REAC215, "REAC215"},
        // {REAC216, "REAC216"},  // Excluded - H- not simulated
        {REAC217, "REAC217"},
        {REAC221a, "REAC221a"},
        {REAC221b, "REAC221b"},
        {REAC222, "REAC222"},
        {REAC223, "REAC223"},
        {REAC224, "REAC224"},
        {REAC225, "REAC225"},
        {REAC226, "REAC226"},
        {REAC227, "REAC227"},
        {REAC228, "REAC228"},
        {REAC229, "REAC229"},
        {REAC2210, "REAC2210"},
        {REAC2211, "REAC2211"},
        {REAC2212, "REAC2212"},
        {REAC2213, "REAC2213"},
        {REAC2214, "REAC2214"},
        {REAC2215, "REAC2215"},
        {REAC2216, "REAC2216"},
        // {REAC2217, "REAC2217"},  // Excluded - not called in collisions.cpp
    };
    
    for (const auto& [reac_id, reac_name] : rr_reactions_1d) {
        for (double Te : Te_grid) {
            double rate = RR(reac_id, Te, Te);  // E=Te for electron-impact
            outFile << "RR," << reac_name << ",Te," << Te << ",none,0," << rate << std::endl;
        }
    }
    
    // REAC218a: Use analytical formula from collisions.cpp (not the unused RR() 2D polynomial)
    // REAC218b: Excluded - not actually used in collisions.cpp
    std::cout << "  Exporting REAC218a (analytical formula from collisions.cpp)..." << std::endl;
    for (double Te : Te_grid) {
        double rate = REAC218a_analytical(Te);
        outFile << "RR," << "REAC218a" << ",Te," << Te << ",none,0," << rate << std::endl;
    }
    
    // =====================================================================
    // RRH2 Reactions (Dirk W. tables)
    // =====================================================================
    std::cout << "  Exporting RRH2 reactions..." << std::endl;
    
    // 2D reactions: DISS, IONI
    std::vector<std::pair<int, std::string>> rrh2_reactions_2d = {
        {DISS, "DISS"},
        {IONI, "IONI"},
    };
    
    for (const auto& [reac_id, reac_name] : rrh2_reactions_2d) {
        for (double ne : ne_slices) {
            for (double Te : Te_grid) {
                if (Te <= 100.0) {  // These tables only go up to 100 eV
                    double rate = RRH2(reac_id, ne, Te);
                    outFile << "RRH2," << reac_name << ",Te," << Te << ",ne," << ne << "," << rate << std::endl;
                }
            }
        }
    }
    
    // 1D reactions: RECO, ELAS
    std::vector<std::pair<int, std::string>> rrh2_reactions_1d = {
        {RECO, "RECO"},
        {ELAS, "ELAS"},
    };
    
    for (const auto& [reac_id, reac_name] : rrh2_reactions_1d) {
        for (double Te : Te_grid) {
            double rate = RRH2(reac_id, 1e11, Te);  // ne doesn't matter for these
            outFile << "RRH2," << reac_name << ",Te," << Te << ",none,0," << rate << std::endl;
        }
    }
    
    // =====================================================================
    // RRHe Reactions (ADAS tables)
    // =====================================================================
    std::cout << "  Exporting RRHe reactions..." << std::endl;
    
    bADAS = true;  // Use ADAS data
    
    std::vector<std::pair<int, std::string>> rrhe_reactions = {
        {IHE1, "IHE1"},
        {IHE2, "IHE2"},
        {RHE2, "RHE2"},
        {RHE3, "RHE3"},
    };
    
    for (const auto& [reac_id, reac_name] : rrhe_reactions) {
        for (double ne : ne_slices) {
            for (double Te : Te_grid) {
                double rate = RRHe(reac_id, ne, Te);
                outFile << "RRHe," << reac_name << ",Te," << Te << ",ne," << ne << "," << rate << std::endl;
            }
        }
    }
    
    // =====================================================================
    // RRCX Reactions (Charge Exchange)
    // =====================================================================
    std::cout << "  Exporting RRCX reactions..." << std::endl;
    
    std::vector<std::pair<int, std::string>> rrcx_reactions = {
        {CXHe3H, "CXHe3H"},
        {CXHe2H, "CXHe2H"},
        {CXHe3He1, "CXHe3He1"},
    };
    
    for (const auto& [reac_id, reac_name] : rrcx_reactions) {
        for (double T2 : T2_slices) {
            for (double T1 : Te_grid) {
                double rate = RRCX(reac_id, T1, T2);  // uppercase function name
                outFile << "RRCX," << reac_name << ",T1," << T1 << ",T2," << T2 << "," << rate << std::endl;
            }
        }
    }
    
    // =====================================================================
    // RRion Reactions (Ion-Neutral)
    // =====================================================================
    std::cout << "  Exporting RRion reactions..." << std::endl;
    
    std::vector<std::pair<int, std::string>> rrion_reactions = {
        {REAC311, "REAC311"},
        {REAC312, "REAC312"},
        {REAC316, "REAC316"},
        {REAC321, "REAC321"},
        {REAC322, "REAC322"},
        {REAC323, "REAC323"},
        {REAC325, "REAC325"},
        {REAC326, "REAC326"},
        {REAC332, "REAC332"},
        {REAC431, "REAC431"},
        {REAC433, "REAC433"},
        // {REAC441, "REAC441"},  // Excluded - commented out in collisions.cpp
        {REAC523, "REAC523"},
        {REAC531, "REAC531"},
        {REAC631, "REAC631"},
    };
    
    for (const auto& [reac_id, reac_name] : rrion_reactions) {
        for (double T2 : T2_slices) {
            for (double T1 : Te_grid) {
                double rate = RRion(reac_id, T1, T2);
                outFile << "RRion," << reac_name << ",T1," << T1 << ",T2," << T2 << "," << rate << std::endl;
            }
        }
    }
    
    // =====================================================================
    // RRel Reactions (Elastic Collisions)
    // =====================================================================
    std::cout << "  Exporting RRel reactions..." << std::endl;
    
    std::vector<std::pair<int, std::string>> rrel_reactions = {
        {HH2, "HH2"},
        {HiH, "HiH"},
        {H2iH, "H2iH"},
        {H3iH, "H3iH"},
        {HeIIH, "HeIIH"},
        {HiH2, "HiH2"},
        {H2iH2, "H2iH2"},
        {H3iH2, "H3iH2"},
        {HeIIH2, "HeIIH2"},
        {HHe, "HHe"},
        {HiHe, "HiHe"},
        {HeIIHe, "HeIIHe"},
        {HeIHeI, "HeIHeI"},
        {H2H2, "H2H2"},
        {HH, "HH"},
    };
    
    for (const auto& [reac_id, reac_name] : rrel_reactions) {
        for (double T2 : T2_slices) {
            for (double T1 : Te_grid) {
                double rate = RRel(reac_id, T1, T2);
                outFile << "RRel," << reac_name << ",T1," << T1 << ",T2," << T2 << "," << rate << std::endl;
            }
        }
    }
    
    outFile.close();
    std::cout << "Export complete!" << std::endl;
    
    return 0;
}

#include "coupledpower.h"
#include "../Vars/simparam.h"
#include "transport.h"
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

void bTOMASIC_func(const double &freq) {
    //# Write to csv file
    //## open and initialize in /tmp
    std::ofstream csvwritefile;
    csvwritefile.open("/tmp/DensAndTemp.csv"); 
    csvwritefile << "Ra,Ne,Te,nHi,THi,nH,nH2,";
    csvwritefile << freq;
    csvwritefile << ",";
    csvwritefile << Prf;
    csvwritefile << "\n";
    //## write arrays
    for (int i=0; i<NMESHP; i++)
        {
            csvwritefile << aR[i];
            csvwritefile << ",";
            // e and H^+
            csvwritefile << nr.ne[i];
            csvwritefile << ",";
            csvwritefile << Tr.Te[i];
            csvwritefile << ",";
            csvwritefile << nr.nHi[i];
            csvwritefile << ",";
            csvwritefile << Tr.THi[i];
            csvwritefile << ",";
            // H and H2
            csvwritefile << nr.nH[i];
            csvwritefile << ",";
            csvwritefile << nr.nH2[i];
            csvwritefile << "\n";
        }
    csvwritefile.close(); 
    //# Execute ngsolve sim and read result
    //## Form command
    std::string BaseFolder = std::getenv("TOMATORSOURCE");
    std::string ScriptPosAdd = "/src/External/1DTOMAS.py";
    std::string ScriptPos = BaseFolder + ScriptPosAdd;
    std::string command = "python ";
    command += ScriptPos;
    //## Execution (also waits for the process to complete)
    cout << "Executing python" << endl;
    system(command.c_str());
    //## read result
    //### Create a vector of <string, double vector> pairs to store the result
    std::vector<std::pair<std::string, std::vector<double>>> result;

    // Create an input filestream
    std::ifstream myFile("/tmp/PowerDeposition.csv");

    // Make sure the file is open
    if(!myFile.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    double val;

    // Read the column names
    if(myFile.good())
    {
        // Extract the first line in the file
        std::getline(myFile, line);

        // Create a stringstream from line
        std::stringstream ss(line);

        // Extract each column name
        while(std::getline(ss, colname, ',')){
            // Initialize and add <colname, int vector> pairs to result
            result.push_back({colname, std::vector<double> {}});
        }
    }

    // Read data, line by line
    while(std::getline(myFile, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);

        // Keep track of the current column index
        int colIdx = 0;

        // Extract each double
        while(ss >> val){

            // Add the current integer to the 'colIdx' column's values vector
            result.at(colIdx).second.push_back(val);

            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();

            // Increment the column index
            colIdx++;
        }
    }
    double ScaleFactor[NMESHP];
    double CutCircleHeight; // toroidal area of element
    double b = 26.0;
    #pragma omp parallel for
    for (int i = 0; i < NMESHP-1; ++i) {
	double FlatArea = pi*(pow(aR[i+1],2.0)-pow(aR[i],2.0));
        if ((78.0-24.0<aR[i]) && (aR[i]<78.0+24.0)) {
            CutCircleHeight = 2.0*b*sqrt(1.0-pow((aR[i]-78.0)/b,2.0));
        }
        else{
            CutCircleHeight = 2.0*b*sqrt(1.0-pow((26.0-1.0)/b,2.0));
        }
        ScaleFactor[i] = 6.241509074461e18/(CutCircleHeight*FlatArea); //W per meshpoint to eV/cm^3s
    }
    for (int id = 0; id < NMESHP; ++id) {
        PRFe_array[id] = ScaleFactor[id]*result.at(0).second[id];  //(eV/cm^3s)
        PRFHi_array[id] = ScaleFactor[id]*result.at(1).second[id]; //(eV/cm^3s)
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
    //# Write to csv file, to make sure it checks out
    //## open and initialize in /tmp
    std::ofstream csvwritedebugfile;
    csvwritedebugfile.open("/tmp/TomatorArrays.csv"); 
    csvwritedebugfile << "id,aR[id],PRFe,PRFHi,nH,nH_2\n";
    //## write arrays
    for (int i=0; i<NMESHP; i++)
        {
            csvwritedebugfile << i;
            csvwritedebugfile << ",";
            csvwritedebugfile << aR[i];
            csvwritedebugfile << ",";
            // e and H^+
            csvwritedebugfile << PRFe_array[i];
            csvwritedebugfile << ",";
            csvwritedebugfile << PRFHi_array[i];
            csvwritedebugfile << ",";
            // H and H2
            csvwritedebugfile << nr.nH[i];
            csvwritedebugfile << ",";
            csvwritedebugfile << nr.nH2[i];
            csvwritedebugfile << "\n";
        }
    csvwritedebugfile.close(); 

}

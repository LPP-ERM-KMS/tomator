#ifndef FILEHANDLER_H
#define FILEHANDLER_H

// Libraries
#include <fstream>

// Variable definitions
#include "../Vars/globalVariables.h"
#include "../Vars/matrixh.h"
#include "../Vars/positions.h"
#include "../Vars/simparam.h"

void infile(string);
void writeHeader(ofstream *outFile, const string *sfilenamechar);
void writeToOutFile(ofstream *outFile);
void writeLast(ofstream *outFile, double tstartcalculation, double currentTime);
void writePhysicalStates(ofstream *outFile);

#endif // FILEHANDLER_H

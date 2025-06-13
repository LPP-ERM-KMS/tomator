#ifndef COUPLEDPOWER_H
#define COUPLEDPOWER_H

// Libraries
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <chrono>

// Variable definitions
#include "../Vars/globalVariables.h"
#include "../Vars/graydat.h"
#include "../Vars/matrixh.h"
#include "../Vars/positions.h"
#include "../Vars/reactionrates.h"
#include "../Vars/simparam.h"

#define CALL_FORTRAN(function) function##_

#define NSORT 5        //+3          // number of sorts of ions;
#define NPOINTS NMESHP // number of mesh points;

extern "C" void CALL_FORTRAN(rfpower)(const double &, int[], int[], int &, double[], int &, double[][NPOINTS * 2], double[], double[][NPOINTS], double &, double &, double[][NPOINTS], int &);

Eigen::VectorXd interpolateData(const Eigen::VectorXd &, const Eigen::VectorXd &, const Eigen::VectorXd &);
std::string trim(const std::string &);
void coupledpower(const double &, double &);

void bproptone_func();
void bgray_func();
void bram_func();
void bnefix_func();
void bfixpowerfrac_func();
void bnopower_func();
void bTOMAS_func();
void bicatlhr_func();
void bmanuel_func(bool);
void bICWC_func();


#endif // COUPLEDPOWER_H

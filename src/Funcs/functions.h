#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// Libraries
#include <cmath>
#include <omp.h>

// Variable definitions
#include "../Vars/globalVariables.h"
#include "../Vars/matrixh.h"
#include "../Vars/positions.h"
#include "../Vars/reactionrates.h"


void init_positions();
void initAuxiliarPos();
double computeN(double, double);
void initializeSpecies();
double ndamp(double);
void limiters();
void vdrift_function();
void bpol_function();

#endif // FUNCTIONS_H

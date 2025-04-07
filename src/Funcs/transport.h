#ifndef TRANSPORT_H
#define TRANSPORT_H

// Libraries
#include <iostream>

using namespace std;

// Variable definitions
#include "../Vars/globalVariables.h"
#include "../Vars/positions.h"
#include "../Vars/solverVariables.h"

void transpCoef();
void transportions(double);
void transpH(double);
void transpHe(double);
void transpC(double);


#endif // TRANSPORT_H

#pragma once

#include "Element.h"

using namespace std;

void StressSPR(double* stress_SPR, double* stressG, double* PrePositions, double* PositionG,
	unsigned int* Ele_NodeNumber, unsigned int NUME, unsigned int NUMNP);

#ifndef _INTERIOR_ENERGY__
#define _INTERIOR_ENERGY__

#include <stdlib.h>
#include "basics.h"

using namespace std;

double InteriorLoopEnergy(int leftSize, int rightSize, int bp_i, int bp_j, int bp_before_i, int bp_before_j, int bp_i_plus1, int bp_j_minus1, int bp_before_i_minus1, int bp_before_j_plus1);

#endif   // _INTERIOR_ENERGY_

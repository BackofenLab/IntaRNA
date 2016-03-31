
#include "stacking_energy.h"

double StackingEnergy(int bp_i, int bp_j, int bp_before_i, int bp_before_j)
{
   double energy;
   
   if ((PairTest(bp_before_i, bp_before_j) == true) && (PairTest(bp_i, bp_j) == true))
      energy = stacking_energies[64*bp_i+16*bp_before_i+4*bp_j+bp_before_j];
   else
      energy = MAX_DOUBLE;

   return energy;
}



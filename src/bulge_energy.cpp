#include "bulge_energy.h"

double BulgeEnergy(int size, int bp_i, int bp_j, int bp_before_i, int bp_before_j)
{
   double energy = 0.0;
   
   // valid base pair test:
   //****************************
   if ((PairTest(bp_before_i, bp_before_j) == false) || (PairTest(bp_i, bp_j) == false))
      return MAX_DOUBLE;
   else
   {
      int bp_before = BP2int(bp_before_i, bp_before_j);

      // no bulges of size 0:
      //******************************************
      if (size == 0)
      {
         cerr << "No Bulge!" << endl;
         exit(1);
      }
      
      // loop_destabilizing_energies:
      //******************************
      if (size <= 30)
        energy += loop_destabilizing_energies[3*size-2];
      else
      {
         energy += loop_destabilizing_energies[3*30-2];
         // energy += 1.75*RT*log((double)(size/30.0)); replaced by
         energy += loop_extrapolation*log((double)(size/30.0));
      }

      // additional bulge-energy (size 1: stacking, >1: nothing, but term. A-U-penalty:
      //**********************************************************************************
      if (size == 1)
         energy += stacking_energies[64*bp_i+16*bp_before_i+4*bp_j+bp_before_j];
      else
      {
         if ((bp_before == 0) || (bp_before == 3) || (bp_before == 4) || (bp_before == 5))
            energy += terminal_AU;
         if ((BP2int(bp_i, bp_j) == 0) || (BP2int(bp_i,bp_j) == 3) || (BP2int(bp_i,bp_j) == 4) || (BP2int(bp_i,bp_j) == 5))
            energy += terminal_AU;
      }
      return energy;
   }
}

#include "dangling_energy.h"

/*****************************************
* 5' jd 3'    5'  i 3'
* 3' i  5'    3' dj 5'
*****************************************/
double ed3(int bp_i, int bp_j, int dangle3)
{
   double energy = MAX_DOUBLE;

   if (PairTest(bp_i, bp_j) == true)
      energy = single_base_stacking_energy[16*bp_j+4*bp_i+dangle3];

   return energy;
}

/*****************************************
* 5' j   3'    5' di 3'
* 3' id  5'    3'  j 5'
*****************************************/
double ed5(int bp_i, int bp_j, int dangle5)
{
   double energy = MAX_DOUBLE;

   if (PairTest(bp_i, bp_j) == true)
      energy = single_base_stacking_energy[64+16*bp_j+4*bp_i+dangle5];

   return energy;
}


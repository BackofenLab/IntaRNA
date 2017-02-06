

#include "interior_energy.h"

/************************************************************

finds the energy of an interior loop, if the assignments of
both closing BPs and all free bases are known

(only energy fragments that are important for the difference 
of two assignments of the loop)

************************************************************/

double InteriorLoopEnergy(int leftSize, int rightSize, int bp_i, int bp_j, int bp_before_i, int bp_before_j, int bp_i_plus1, int bp_j_minus1, int bp_before_i_minus1, int bp_before_j_plus1)
{
   /*        bp_i_plus1 ........ bp_before_i_minus1
   ------bp_i                                      bp_before_i------------
   
   ------bp_j                                      bp_before_j------------
             bp_j_minus1..........bp_before_j_plus1
             
   bp_i_plus1 and bp_before_i_minus1 might be the same nucleotide (if leftSize == 1)
   */
              
   // valid base pair test:
   //****************************
   if ((PairTest(bp_before_i, bp_before_j) == false) || (PairTest(bp_i, bp_j) == false))
      return MAX_DOUBLE;
   else
   {
      double energy = 0.0;
      // double asym = 0.0;
      int x1, x2, y1, y2, i1, i2, i3, i4;
      int size = leftSize + rightSize;

      int bp_new = BP2int(bp_i, bp_j);
      int bp_before = BP2int(bp_before_i, bp_before_j);

      // no bulges and loops of size 0:
      //--------------------------------------------
      if ((leftSize == 0) || (rightSize == 0))
      {
         printf("No interiorLoop!\n");
         exit(1);
      }

      // Special cases:
      //------------------------------------------------------
      if ((leftSize == 1) && (rightSize == 1))
      {
         x1 = bp_i_plus1;
         y1 = bp_j_minus1;
         energy += interior_loop_1_1_energy[96*bp_new+24*x1+4*bp_before+y1];
      }
      else if ((leftSize == 1) && (rightSize == 2)) /*x(i-1)-x(i)=2 and y(i)-y(i-1)=3*/
      {
         x1 = bp_i_plus1;
         y1 = bp_before_j_plus1;
         y2 = bp_j_minus1;
         energy += interior_loop_1_2_energy[384*bp_new+96*y1+24*x1+4*bp_before+y2];
      }
      else if ((leftSize == 2) && (rightSize == 1)) /*x(i-1)-x(i)=3 and y(i)-y(i-1)=2*/
      {
         x1 = bp_i_plus1;
         x2 = bp_before_i_minus1;
         y1 = bp_j_minus1;
         energy += interior_loop_1_2_energy[384*bp_before+96*x1+24*y1+4*bp_new+x2];
      }
      else if ((leftSize == 2) && (rightSize == 2))
      {
         x1 = bp_i_plus1;
         x2 = bp_before_i_minus1;
         y1 = bp_before_j_plus1;
         y2 = bp_j_minus1;
         energy += interior_loop_2_2_energy[1536*bp_new+256*bp_before+64*x1+16*y2+4*x2+y1];
      }
      else
      {
         if (size <= 30)
            energy += loop_destabilizing_energies[3*(size-1)];
         else
         {
            energy += loop_destabilizing_energies[3*(30-1)];
            // energy += 1.75*RT*log((double)(size/30.0)); replaced by
            energy += loop_extrapolation*log((double)(size/30.0));
         }

         // terminal mismatches on both closings:
         //----------------------------------------
         // Attention: there are dependencies, if the IL has size 1 at one side of the loop,
         //            then this base is involved in the terminal mismatches at both closings
         if (leftSize == 1)      //        bp_before_i - bp_before_j
         {                       //  i1                                      i2
                                 //                                          i3
                                 //               bp_i - bp_j
            i1 = bp_i_plus1;
            i2 = bp_before_j_plus1;
            i3 = bp_j_minus1;
            energy += mismatch_energies_interior[64*bp_before_j+16*i2+4*bp_before_i+i1] + mismatch_energies_interior[64*bp_i+16*i1+4*bp_j+i3];
         }
         else if (rightSize == 1)//        bp_before_i - bp_before_j
         {                       //  i1                                      i3
                                 //  i2
                                 //               bp_i - bp_j
            i1 = bp_before_i_minus1;
            i2 = bp_i_plus1;
            i3 = bp_j_minus1;
            energy += mismatch_energies_interior[64*bp_before_j+16*i3+4*bp_before_i+i1] + mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i3];
         }
         else                    //        bp_before_i - bp_before_j
         {                       //  i1                                      i3
                                 //  i2                                      i4
                                 //               bp_i - bp_j
            i1 = bp_before_i_minus1;
            i2 = bp_i_plus1;
            i3 = bp_before_j_plus1;
            i4 = bp_j_minus1;
            energy += mismatch_energies_interior[64*bp_before_j+16*i3+4*bp_before_i+i1] + mismatch_energies_interior[64*bp_i+16*i2+4*bp_j+i4];
         }

         // asymmetry-penalty:
         //------------------------
         // (for 1x2 IL no extra penalty)
         if (leftSize != rightSize)
         {
            // asym = 0.5 * abs(leftSize-rightSize); if (asym > 3.0) asym = 3.0; energy+=asym; replaced by
            energy += Minimum( max_ninio_correction, (ninio_correction * abs(leftSize-rightSize)) );
         }
      }
      return energy;
   }
}



